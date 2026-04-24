# -*- coding: utf-8 -*-
"""
无传感器控制 Demo — Active Flux 估计 vs 开环积分对比
====================================================
基于 tutorials_ep6_svpwm.py 的 FOC 仿真框架，实现：

1. **Active Flux Observer (with PI correction)**
   - 电压模型积分 + PI 校正（防止直流漂移）
   - 从 active flux 提取转子位置和转速

2. **Open-loop pure integrator (no correction)**
   - 纯电压模型积分（无修正），用于展示直流漂移问题

3. **速度估计方法对比**
   - 方法 (i):  1 阶低通滤波 (LPF) — Δθ/Δt + LPF
   - 方法 (ii): 4 阶 ESO — 把 θ_AF 当 "测量"，速度和负载转矩作为观测器内部状态

4. **对比图**
   - 角度估计 vs 编码器真值
   - 转速估计对比（LPF vs 4th-order ESO）
   - 磁链轨迹（αβ 平面）
   - ESO 负载转矩估计
   - 角度误差时域波形

原理说明
--------
定子磁链电压模型：
    dψ_s/dt = u_s - R·i_s

Active flux（有功磁链）：
    ψ_AF = ψ_s - Lq·i_s

转子位置由 ψ_AF 的方向确定：
    θ_est = arctan2(ψ_AF_β, ψ_AF_α)

速度估计 — 两种路径：
  路径 1: ω = LPF(Δθ/Δt)                   — 1 阶，等效低通滤波
  路径 2: 4th-order ESO(θ_AF) → xS[1] = ω  — 4 阶，内建速度和负载估计
"""

# %%
############################################# PACKAGES
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for headless execution
from pylab import np, plt, mpl
import copy
import time as _time

# Import the core simulation framework
from tutorials_ep6_svpwm import (
    The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental,
)
from tuner import (
    get_coeffs_dc_motor_current_regulator,
    get_coeffs_dc_motor_SPEED_regulator,
)


# ======================================================================
# Helper: angle wrapping to [-π, π]
# ======================================================================
def wrap_angle(theta):
    """Wrap angle (or array of angles) to [-π, π]."""
    return (theta + np.pi) % (2 * np.pi) - np.pi


def angle_error(a, b):
    """Element-wise smallest signed angle difference a - b, in [-π, π]."""
    d = a - b
    return (d + np.pi) % (2 * np.pi) - np.pi


def angle_diff_scalar(a, b):
    """Scalar version of smallest signed angle difference a - b, in [-π, π].
    Both a, b should be in [0, 2π] or any range (will be wrapped)."""
    a = a % (2 * np.pi)
    b = b % (2 * np.pi)
    d1 = a - b
    if d1 > 0:
        d2 = a - (b + 2 * np.pi)
    else:
        d2 = (2 * np.pi + a) - b
    if abs(d1) < abs(d2):
        return d1
    else:
        return d2


# ======================================================================
# 4th-order ESO dynamics (pure Python, no Numba dependency)
# ======================================================================
class SpeedESO:
    """
    4th-order Extended State Observer for speed estimation from angle.

    状态定义：
        x[0] = θ_est    (estimated angle, [elec. rad])
        x[1] = ω_est    (estimated elec. angular velocity, [rad/s])
        x[2] = d_est    (estimated total disturbance / npp*Js, [rad/s^2])
        x[3] = dd_est   (estimated disturbance rate, [rad/s^3])

    观测器动力学 (continuous-time):
        dx[0]/dt = ell1 * (θ_meas - x[0]) + x[1]
        dx[1]/dt = ell2 * (θ_meas - x[0]) + (Tem + x_load) * npp/Js
        dx[2]/dt = ell3 * (θ_meas - x[0]) + x[3]
        dx[3]/dt = ell4 * (θ_meas - x[0])

    其中 (Tem + x_load) * npp/Js 是加速度项，在这里我们不使用
    电磁转矩前馈（因为我们是 sensorless，不一定有准确的 Tem），
    所以设 Tem=0，让 ESO 的扰动状态 x[2] 去估计全部的加速度/扰动。

    Gains: ell1..ell4 via pole placement at -omega_ob (4 重极点):
        ell1 = 4 * ω_ob
        ell2 = 6 * ω_ob²
        ell3 = 4 * ω_ob³ * Js / npp
        ell4 = ω_ob⁴
    """

    def __init__(self, omega_ob, npp, Js, dt):
        self.npp = npp
        self.Js = Js
        self.dt = dt
        self.omega_ob = omega_ob

        # 4th-order observer gains (4 repeated poles at -omega_ob)
        # Note: ell3 and ell4 must include the Js/npp factor to account for
        # the npp/Js scaling in the plant model (fx[1] = ... + x[2]*npp/Js).
        # This matches the corrected formula in demo_4th_order_ESO.py.
        self.ell1 = 4 * omega_ob
        self.ell2 = 6 * omega_ob**2
        self.ell3 = 4 * omega_ob**3 * Js / npp
        self.ell4 =     omega_ob**4 * Js / npp

        # State: [θ_est, ω_est, d_est, dd_est]
        self.x = np.zeros(4)

    def reset(self, theta0=0.0):
        self.x = np.zeros(4)
        self.x[0] = theta0

    def _dynamics(self, x, theta_meas, Tem_ff=0.0):
        """Compute dx/dt given current state and measurement.

        Parameters
        ----------
        Tem_ff : float
            Electromagnetic torque feedforward [N·m].
            With Tem_ff=0: x[2] estimates net torque (Tem - TLoad) ≈ 0 in SS.
            With Tem_ff=CTRL.Tem: x[2] estimates -TLoad (the unknown disturbance).
        """
        # Wrap x[0] to [-π, π] before computing angle difference
        # This is critical for RK4 intermediate stages where x[0] may drift.
        x0_wrapped = (x[0] + np.pi) % (2 * np.pi) - np.pi
        output_error = angle_diff_scalar(theta_meas, x0_wrapped)

        fx = np.zeros(4)
        fx[0] = self.ell1 * output_error + x[1]
        fx[1] = self.ell2 * output_error + (Tem_ff + x[2]) * self.npp / self.Js
        fx[2] = self.ell3 * output_error + x[3]
        fx[3] = self.ell4 * output_error
        return fx

    def step(self, theta_meas, Tem_ff=0.0):
        """Advance one step using RK4 integration.

        Parameters
        ----------
        theta_meas : float
            Measured (or estimated) electrical angle [rad].
        Tem_ff : float
            Electromagnetic torque feedforward [N·m]. If provided,
            x[2] will estimate -TLoad. If 0, x[2] estimates Tem-TLoad.
        """
        hs = self.dt
        x = self.x.copy()  # must copy to avoid modifying state during RK4

        k1 = self._dynamics(x, theta_meas, Tem_ff) * hs
        k2 = self._dynamics(x + k1 * 0.5, theta_meas, Tem_ff) * hs
        k3 = self._dynamics(x + k2 * 0.5, theta_meas, Tem_ff) * hs
        k4 = self._dynamics(x + k3, theta_meas, Tem_ff) * hs

        self.x = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

        # Wrap angle to [-π, π]
        self.x[0] = (self.x[0] + np.pi) % (2 * np.pi) - np.pi

    @property
    def theta_est(self):
        return self.x[0]

    @property
    def omega_est(self):
        """Estimated electrical angular velocity [rad/s]."""
        return self.x[1]

    @property
    def omega_elec(self):
        """Estimated electrical angular velocity [rad/s]."""
        return self.x[1]

    @property
    def omega_rpm(self):
        """Estimated mechanical speed [rpm]."""
        return self.x[1] / (2 * np.pi * self.npp) * 60

    @property
    def disturbance_est(self):
        """Estimated total disturbance (proportional to load torque) [N·m equivalent]."""
        return self.x[2]


class NaturalSpeedObserver:
    """3rd-order Luenberger equivalent of the PID Natural Speed Observer.
    Estimates speed and load torque directly from measured currents and commanded voltage.
    """
    def __init__(self, omega_ob=200.0, npp=4, Js=0.001, R=1.0, Ld=0.001, Lq=0.001, KE=0.1, dt=1e-4):
        self.omega_ob = omega_ob
        self.npp = npp
        self.Js = Js
        self.R = R
        self.Ld = Ld
        self.Lq = Lq
        self.KE = KE
        self.dt = dt

        # State: [iq_est, omega_r_elec_est, TL_est]
        self.x = np.zeros(3)

    def reset(self, omega0=0.0):
        self.x = np.zeros(3)
        self.x[1] = omega0

    def _dynamics(self, x, iq_meas, uq_cmd, id_meas):
        iq_est, omega_est, TL_est = x[0], x[1], x[2]
        
        KActive = (self.Ld - self.Lq) * id_meas + self.KE
        Phi = self.Lq * id_meas + KActive
        
        # Prevent division by close-to-zero Phi
        Phi_safe = Phi if abs(Phi) > 1e-5 else 1e-5 * np.sign(Phi) if Phi != 0 else 1e-5

        kt = 1.5 * (self.npp ** 2) * KActive

        # Luenberger gains from characteristic polynomial (s + ω_ob)³:
        #   l1 = 3·ω_ob - R/Lq
        #   l2 = kt/Js - 3·ω_ob²·Lq/Φ
        #   l3 = ω_ob³·Lq·Js / (npp·Φ)   [must be POSITIVE for stability]
        l1 = 3 * self.omega_ob - self.R / self.Lq
        l2 = kt / self.Js - 3 * self.Lq * (self.omega_ob ** 2) / Phi_safe
        l3 = self.Lq * self.Js * (self.omega_ob ** 3) / (self.npp * Phi_safe)

        err = iq_meas - iq_est

        fx = np.zeros(3)
        fx[0] = -self.R / self.Lq * iq_est - Phi / self.Lq * omega_est + uq_cmd / self.Lq + l1 * err
        fx[1] = kt / self.Js * iq_est - self.npp / self.Js * TL_est + l2 * err
        fx[2] = l3 * err

        return fx

    def step(self, iq_meas, uq_cmd, id_meas):
        """Advance one step using RK4."""
        dt = self.dt
        k1 = self._dynamics(self.x, iq_meas, uq_cmd, id_meas)
        k2 = self._dynamics(self.x + 0.5 * dt * k1, iq_meas, uq_cmd, id_meas)
        k3 = self._dynamics(self.x + 0.5 * dt * k2, iq_meas, uq_cmd, id_meas)
        k4 = self._dynamics(self.x + dt * k3, iq_meas, uq_cmd, id_meas)
        self.x += (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

    @property
    def omega_est(self):
        return self.x[1]

    @property
    def omega_elec(self):
        return self.x[1]

    @property
    def omega_rpm(self):
        return self.x[1] / (2 * np.pi * self.npp) * 60

    @property
    def disturbance_est(self):
        """Estimated total disturbance (opposite sign of TL for compatible ESO FF injection)."""
        return -self.x[2]



# ======================================================================
# Sensorless simulation (Active Flux observer + open-loop integrator)
# ======================================================================
def run_sensorless_demo(d,
                        zeta=15, CLBW_Hz=1000,
                        af_Kp=500.0, af_Ki=5000.0,
                        R_mismatch_factor=1.5,
                        voltage_offset_alpha=0.02,
                        voltage_offset_beta=-0.015,
                        eso_omega_ob=200.0,
                        speed_observer='eso',
                        use_sensorless_speed=False,
                        use_sensorless_torque_ff=False,
                        pure_p_current=False,
                        verbose=True):
    """
    Run a sensorless control demo simulation.

    Two flux observers run in parallel alongside the true encoder:
      (A) Active Flux observer with PI correction
      (B) Pure open-loop integrator (no correction)

    Two speed estimation methods run on observer (A)'s angle output:
      (i)  1st-order: Δθ/Δt + LPF (simple differentiation)
      (ii) 4th-order ESO: feeds θ_AF as measurement, extracts speed and
           load torque as internal observer states

    Sensorless closed-loop switches
    --------------------------------
    use_sensorless_speed : bool
        If True, the FOC speed loop uses the ESO-estimated speed ω_est
        instead of the encoder speed. The Park transformation angle also
        switches to θ_AF (Active Flux estimate).
    use_sensorless_torque_ff : bool
        If True, the ESO's load torque estimate is fed forward into the
        iq command via CTRL.total_disrubance_feedforward.

    Parameters
    ----------
    d : dict
        Motor / controller parameter dictionary.
    eso_omega_ob : float
        Observer bandwidth [rad/s] for the 4th-order speed ESO.
    """
    dd = copy.deepcopy(d)

    # ------- PI tuning -------
    R_true = dd['init_R']
    R_obs  = R_true * R_mismatch_factor
    L      = dd['init_Lq']
    J_s    = dd['init_Js']
    n_pp   = dd['init_npp']
    KE     = dd['init_KE']
    KA     = KE
    CL_TS  = dd['CL_TS']
    VL_TS  = dd['CL_TS'] * dd['VL_EXE_PER_CL_EXE']

    currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R_true, L, CLBW_Hz)
    currentBandwidth_radPerSec = currentKp / L
    speedKp, speedKi = get_coeffs_dc_motor_SPEED_regulator(
        J_s, n_pp, KA, zeta, currentBandwidth_radPerSec
    )

    dd['CL_SERIES_KP'] = currentKp
    dd['CL_SERIES_KI'] = currentKi
    dd['VL_SERIES_KP'] = speedKp
    dd['VL_SERIES_KI'] = speedKi

    # Determine if we need per-step injection
    sensorless_closed_loop = use_sensorless_speed or use_sensorless_torque_ff

    if verbose:
        print(f'  R_true = {R_true:.4f} Ω, R_obs = {R_obs:.4f} Ω (factor={R_mismatch_factor})')
        print(f'  Voltage offsets: α={voltage_offset_alpha:.4f} V, β={voltage_offset_beta:.4f} V')
        print(f'  KE = {KE} Wb, L = {L*1e3:.4f} mH')
        print(f'  ESO bandwidth: ω_ob = {eso_omega_ob} rad/s')
        if sensorless_closed_loop:
            print(f'  *** Sensorless closed-loop mode ***')
            print(f'      use_sensorless_speed    = {use_sensorless_speed}')
            print(f'      use_sensorless_torque_ff = {use_sensorless_torque_ff}')

    # ------- Build simulation objects -------
    CTRL = The_Motor_Controller(
        CL_TS=dd['CL_TS'],
        VL_TS=VL_TS,
        init_npp=dd['init_npp'],
        init_IN=dd['init_IN'],
        init_R=dd['init_R'],
        init_Ld=dd['init_Ld'],
        init_Lq=dd['init_Lq'],
        init_KE=dd['init_KE'],
        init_Rreq=dd['init_Rreq'],
        init_Js=dd['init_Js'],
        DC_BUS_VOLTAGE=dd['DC_BUS_VOLTAGE'],
    )
    CTRL.bool_apply_decoupling_voltages_to_current_regulation = dd['CTRL.bool_apply_decoupling_voltages_to_current_regulation']
    CTRL.bool_apply_sweeping_frequency_excitation = False
    CTRL.bool_overwrite_speed_commands = True
    CTRL.bool_zero_id_control = dd['CTRL.bool_zero_id_control']
    CTRL.bool_apply_speed_closed_loop_control = True

    # ---- Sensorless injection setup ----
    if sensorless_closed_loop:
        # Enable the built-in speed observer path so that
        # CTRL.omega_r_elec reads from CTRL.xS[1], and
        # CTRL.total_disrubance_feedforward reads from CTRL.xS[2].
        CTRL.index_separate_speed_estimation = 1

        # Zero out the built-in observer gains — we will write xS externally.
        CTRL.ell1 = 0.0
        CTRL.ell2 = 0.0
        CTRL.ell3 = 0.0
        CTRL.ell4 = 0.0

        # Enable sensorless Park angle (θ from vartheta_d instead of encoder)
        if use_sensorless_speed:
            CTRL.bool_use_sensorless_theta = 1

        # Enable disturbance feedforward
        if use_sensorless_torque_ff:
            CTRL.use_disturbance_feedforward_rejection = 1
        else:
            CTRL.use_disturbance_feedforward_rejection = 0

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=dd['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])

    # Current & speed regulators
    if dd.get('CTRL.bool_apply_decoupling_voltages_to_current_regulation', False):
        local_Ki_factor = 1.0
    else:
        local_Ki_factor = dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False', 10)
        
    if pure_p_current:
        local_Ki_factor = 0.0
        print("  [Notice] Current regulator is set to PURE PROPORTIONAL (Ki=0).")

    reg_id = The_PID_Regulator(currentKp, currentKp * currentKi * local_Ki_factor, 0.0, 0.0,
                                dd['DC_BUS_VOLTAGE'] / 1.732, dd['DC_BUS_VOLTAGE'] / 1.732, CL_TS)
    reg_iq = The_PID_Regulator(currentKp, currentKp * currentKi * local_Ki_factor, 0.0, 0.0,
                                dd['DC_BUS_VOLTAGE'] / 1.732, dd['DC_BUS_VOLTAGE'] / 1.732, CL_TS)
    reg_speed = The_PID_Regulator(speedKp, speedKp * speedKi, 0.0, 0.0,
                                   dd['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * dd['init_IN'],
                                   dd['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * dd['init_IN'],
                                   VL_TS)

    MACHINE_TS = CL_TS / ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD

    # ============================================================
    #  Observer states
    # ============================================================
    # (A) Active Flux Observer with PI correction
    af_psi_s = np.zeros(2)         # stator flux [α, β]
    af_psi_s[0] = KE               # init to approximate PM flux
    af_corr_int = np.zeros(2)      # PI correction integrator

    # (B) Open-loop pure integrator
    ol_psi_s = np.zeros(2)
    ol_psi_s[0] = KE

    # ============================================================
    #  Speed estimation: two methods
    # ============================================================
    # Method (i): 1st-order LPF on dθ/dt
    af_omega_lpf = 0.0
    af_theta_prev = 0.0
    lpf_tau   = 0.002  # LPF time constant [s]
    lpf_alpha = CL_TS / (lpf_tau + CL_TS)

    # Method (ii): Dynamics Observer (ESO or NSO)
    if speed_observer == 'eso':
        dyn_obs = SpeedESO(omega_ob=eso_omega_ob, npp=n_pp, Js=J_s, dt=CL_TS)
    else:
        # NSO initialized with the same omega_ob for fair tracking bandwidth comparison
        dyn_obs = NaturalSpeedObserver(omega_ob=eso_omega_ob, npp=n_pp, Js=J_s, 
                                       R=R_obs, Ld=dd['init_Ld'], Lq=dd['init_Lq'], KE=KE, dt=CL_TS)

    # For open-loop integrator speed (still uses LPF, for reference)
    ol_omega_est = 0.0
    ol_theta_prev = 0.0

    # ============================================================
    #  Simulation parameters
    # ============================================================
    total_time = dd['TIME_SLICE'] * dd['NUMBER_OF_SLICES']
    controller_down_sampling = int(CL_TS / MACHINE_TS)

    # For sensorless closed-loop, run one CL_TS per slice
    if sensorless_closed_loop:
        sim_TIME_SLICE = CL_TS
        sim_NUMBER_OF_SLICES = int(total_time / CL_TS)
    else:
        sim_TIME_SLICE = dd['TIME_SLICE']
        sim_NUMBER_OF_SLICES = dd['NUMBER_OF_SLICES']

    # Storage arrays (at CL_TS rate)
    n_ctrl = int(total_time / CL_TS) + 100
    t_trace       = np.zeros(n_ctrl)
    theta_true    = np.zeros(n_ctrl)
    omega_true    = np.zeros(n_ctrl)
    theta_af      = np.zeros(n_ctrl)
    omega_af_lpf  = np.zeros(n_ctrl)  # speed via LPF
    omega_af_eso  = np.zeros(n_ctrl)  # speed via 4th-order ESO
    theta_af_eso  = np.zeros(n_ctrl)  # ESO's own angle estimate
    load_eso      = np.zeros(n_ctrl)  # ESO disturbance estimate
    load_true     = np.zeros(n_ctrl)  # true load torque
    psi_af_ab     = np.zeros((n_ctrl, 2))
    theta_ol      = np.zeros(n_ctrl)
    omega_ol      = np.zeros(n_ctrl)
    psi_ol_ab     = np.zeros((n_ctrl, 2))
    cmd_rpm_trace = np.zeros(n_ctrl)

    # ============================================================
    #  Main simulation loop
    # ============================================================
    if verbose:
        print('Starting sensorless demo simulation...')
    t_start_sim = _time.time()

    CTRL.cmd_rpm = 0.0
    ctrl_idx = 0

    for slice_idx in range(sim_NUMBER_OF_SLICES):
        t0_slice = slice_idx * sim_TIME_SLICE
        t_mid = t0_slice + sim_TIME_SLICE * 0.5

        # ---- Speed command & load profile ----
        if t_mid < 0.3:
            CTRL.cmd_rpm = 200 * min(t_mid / 0.3, 1.0)
            ACM.TLoad = 0.0
        elif t_mid < 0.6:
            CTRL.cmd_rpm = 200
            ACM.TLoad = 0.0
        elif t_mid < 0.8:
            CTRL.cmd_rpm = 200
            ACM.TLoad = 0.15    # step load applied
        elif t_mid < 1.0:
            CTRL.cmd_rpm = 200
            ACM.TLoad = 0.0     # load removed
        elif t_mid < 1.3:
            CTRL.cmd_rpm = -200
            ACM.TLoad = 0.0
        else:
            CTRL.cmd_rpm = -200
            ACM.TLoad = 0.0

        # ---- Inject sensorless observer outputs BEFORE the next control step ----
        if sensorless_closed_loop and ctrl_idx > 0:
            # Inject ESO-estimated angle → CTRL.vartheta_d (used for Park if enabled)
            CTRL.vartheta_d = theta_af[ctrl_idx - 1]  # latest AF angle estimate

            # Inject Observer-estimated speed → CTRL.xS[1] (used as omega_r_elec)
            if use_sensorless_speed:
                CTRL.xS[1] = dyn_obs.omega_elec  # Observer's ω in elec. rad/s

            # Inject Observer disturbance → CTRL.xS[2] (used as feedforward)
            if use_sensorless_torque_ff:
                CTRL.xS[2] = dyn_obs.disturbance_est

        # ---- Run one simulation slice ----
        machine_times, watch_data = ACMSimPyIncremental(
            t0=t0_slice, TIME=sim_TIME_SLICE,
            ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed
        )

        # Down-sample to CL_TS for observer stepping
        for k in range(0, len(machine_times), controller_down_sampling):
            if ctrl_idx >= n_ctrl:
                break

            t_now = machine_times[k]
            t_trace[ctrl_idx] = t_now

            # True values
            theta_true[ctrl_idx] = watch_data[0][k]
            omega_true[ctrl_idx] = watch_data[1][k]  # rpm
            cmd_rpm_trace[ctrl_idx] = watch_data[12][k]
            load_true[ctrl_idx] = ACM.TLoad

            # Measured currents and commanded voltages
            i_alpha = watch_data[6][k]
            i_beta  = watch_data[7][k]
            u_alpha = watch_data[28][k] + voltage_offset_alpha
            u_beta  = watch_data[29][k] + voltage_offset_beta

            # ============================================================
            # (A) Active Flux Observer with PI amplitude correction
            # ============================================================
            af_active_alpha = af_psi_s[0] - L * i_alpha
            af_active_beta  = af_psi_s[1] - L * i_beta
            af_active_amp   = np.sqrt(af_active_alpha**2 + af_active_beta**2)

            psi_af_error = KE - af_active_amp
            af_corr_int[0] += af_Ki * psi_af_error * CL_TS
            af_corr_int[1] += af_Ki * psi_af_error * CL_TS

            if af_active_amp > 1e-10:
                af_unit_alpha = af_active_alpha / af_active_amp
                af_unit_beta  = af_active_beta  / af_active_amp
            else:
                af_unit_alpha = 1.0
                af_unit_beta  = 0.0

            corr_alpha = (af_Kp * psi_af_error + af_corr_int[0]) * af_unit_alpha
            corr_beta  = (af_Kp * psi_af_error + af_corr_int[1]) * af_unit_beta

            af_psi_s[0] += CL_TS * (u_alpha - R_obs * i_alpha + corr_alpha)
            af_psi_s[1] += CL_TS * (u_beta  - R_obs * i_beta  + corr_beta)

            af_active_alpha = af_psi_s[0] - L * i_alpha
            af_active_beta  = af_psi_s[1] - L * i_beta
            theta_af_now = np.arctan2(af_active_beta, af_active_alpha)

            # ---- Speed method (i): LPF on dθ/dt ----
            d_theta_af = angle_error(np.array([theta_af_now]), np.array([af_theta_prev]))[0]
            raw_omega_af = d_theta_af / CL_TS
            af_omega_lpf_val = af_omega_lpf * (1 - lpf_alpha) + raw_omega_af * lpf_alpha
            af_omega_lpf = af_omega_lpf_val
            af_theta_prev = theta_af_now

            # ---- Speed method (ii): Dynamics Observer (4th-order ESO or NSO) ----
            if speed_observer == 'eso':
                Tem_ctrl = watch_data[27][k]  # CTRL.Tem (controller-computed torque)
                dyn_obs.step(theta_af_now, Tem_ff=Tem_ctrl)
                theta_af_eso_val = dyn_obs.theta_est
            else:
                # Use numpy properties from watch_data to avoid massive Numba PyObject boxing overhead
                # i_alpha = watch_data[6][k], i_beta = watch_data[7][k]
                i_a = watch_data[6][k]
                i_b = watch_data[7][k]
                u_a = watch_data[28][k] + voltage_offset_alpha  # CTRL.cmd_uab[0]
                u_b = watch_data[29][k] + voltage_offset_beta   # CTRL.cmd_uab[1]
                
                # Manual Park transform relying on the AF angle
                cos_th = np.cos(theta_af_now)
                sin_th = np.sin(theta_af_now)
                
                id_meas = i_a * cos_th + i_b * sin_th
                iq_meas = -i_a * sin_th + i_b * cos_th
                uq_cmd  = -u_a * sin_th + u_b * cos_th
                
                dyn_obs.step(iq_meas, uq_cmd, id_meas)
                theta_af_eso_val = theta_af_now  # NSO doesn't have an internal angle state, just log the AF angle

            # Store results
            theta_af[ctrl_idx]     = theta_af_now
            omega_af_lpf[ctrl_idx] = af_omega_lpf_val / (2 * np.pi * n_pp) * 60  # rpm
            omega_af_eso[ctrl_idx] = dyn_obs.omega_rpm
            theta_af_eso[ctrl_idx] = theta_af_eso_val
            load_eso[ctrl_idx]     = dyn_obs.disturbance_est  # disturbance state
            psi_af_ab[ctrl_idx]    = [af_active_alpha, af_active_beta]

            # ============================================================
            # (B) Open-loop pure integrator (no correction)
            # ============================================================
            ol_psi_s[0] += CL_TS * (u_alpha - R_obs * i_alpha)
            ol_psi_s[1] += CL_TS * (u_beta  - R_obs * i_beta)

            ol_active_alpha = ol_psi_s[0] - L * i_alpha
            ol_active_beta  = ol_psi_s[1] - L * i_beta
            theta_ol_now = np.arctan2(ol_active_beta, ol_active_alpha)

            d_theta_ol = angle_error(np.array([theta_ol_now]), np.array([ol_theta_prev]))[0]
            raw_omega_ol = d_theta_ol / CL_TS
            ol_omega_est = ol_omega_est * (1 - lpf_alpha) + raw_omega_ol * lpf_alpha
            ol_theta_prev = theta_ol_now

            theta_ol[ctrl_idx]  = theta_ol_now
            omega_ol[ctrl_idx]  = ol_omega_est / (2 * np.pi * n_pp) * 60
            psi_ol_ab[ctrl_idx] = [ol_active_alpha, ol_active_beta]

            ctrl_idx += 1

    elapsed = _time.time() - t_start_sim
    if verbose:
        print(f'  Simulation completed in {elapsed:.2f} s, {ctrl_idx} control steps.')

    # Trim arrays
    n = ctrl_idx
    results = {
        't':           t_trace[:n],
        'theta_true':  theta_true[:n],
        'omega_true':  omega_true[:n],
        'cmd_rpm':     cmd_rpm_trace[:n],
        'theta_af':    theta_af[:n],
        'omega_af_lpf': omega_af_lpf[:n],
        'omega_af_eso': omega_af_eso[:n],
        'theta_af_eso': theta_af_eso[:n],
        'load_eso':    load_eso[:n],
        'load_true':   load_true[:n],
        'psi_af_ab':   psi_af_ab[:n],
        'theta_ol':    theta_ol[:n],
        'omega_ol':    omega_ol[:n],
        'psi_ol_ab':   psi_ol_ab[:n],
        'n_pp':        n_pp,
        'KE':          KE,
        'Js':          J_s,
        'R_true':      R_true,
        'R_obs':       R_obs,
        'v_offset_a':  voltage_offset_alpha,
        'v_offset_b':  voltage_offset_beta,
        'eso_omega_ob': eso_omega_ob,
        'speed_observer': speed_observer,
        'use_sensorless_speed': use_sensorless_speed,
        'use_sensorless_torque_ff': use_sensorless_torque_ff,
    }
    return results


# ======================================================================
# Plotting
# ======================================================================
def plot_sensorless_results(results, save_path='fig_sensorless_demo'):
    """Generate comprehensive comparison plots."""

    plt.style.use('bmh')
    mpl.rc('font', family='Times New Roman', size=10.0)
    mpl.rc('legend', fontsize=8)
    mpl.rcParams['lines.linewidth'] = 1.0
    mpl.rcParams['mathtext.fontset'] = 'stix'

    t    = results['t']
    n_pp = results['n_pp']
    KE   = results['KE']

    err_af = np.degrees(angle_error(results['theta_af'], results['theta_true']))
    err_ol = np.degrees(angle_error(results['theta_ol'], results['theta_true']))

    n_ss = len(t) // 5  # skip first 20% for RMS

    # ================================================================
    # Figure 1: Angle estimation & errors (4 subplots)
    # ================================================================
    fig1, axes = plt.subplots(4, 1, dpi=150, facecolor='w', figsize=(13, 14), sharex=True)

    # (a) Speed command and true speed
    ax = axes[0]
    ax.plot(t, results['cmd_rpm'], 'k--', alpha=0.4, linewidth=0.8, label='cmd rpm')
    ax.plot(t, results['omega_true'], '#1f77b4', linewidth=1.2, label='True speed (encoder)')
    ax.set_ylabel('Speed [rpm]')
    ax.set_title(
        'Sensorless Control Demo: Active Flux Observer vs Open-Loop Integrator\n'
        f'(R mismatch: ×{results["R_obs"]/results["R_true"]:.1f}, '
        f'voltage offsets: α={results["v_offset_a"]:.3f}V, β={results["v_offset_b"]:.3f}V)',
        fontsize=11, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # (b) Angle: Active Flux vs encoder
    ax = axes[1]
    theta_true_w = np.degrees(wrap_angle(results['theta_true']))
    theta_af_w   = np.degrees(wrap_angle(results['theta_af']))
    ax.plot(t, theta_true_w, '#1f77b4', linewidth=0.6, alpha=0.5, label=r'$\theta_{true}$ (encoder)')
    ax.plot(t, theta_af_w,   '#d62728', linewidth=0.6, label=r'$\theta_{AF}$ (Active Flux)')
    ax.set_ylabel('Angle [deg]')
    ax.set_title('Angle Estimation — Active Flux Observer (with PI correction)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # (c) Angle: Open-loop vs encoder
    ax = axes[2]
    theta_ol_w = np.degrees(wrap_angle(results['theta_ol']))
    ax.plot(t, theta_true_w, '#1f77b4', linewidth=0.6, alpha=0.5, label=r'$\theta_{true}$ (encoder)')
    ax.plot(t, theta_ol_w,   '#2ca02c', linewidth=0.6, label=r'$\theta_{OL}$ (Open-loop)')
    ax.set_ylabel('Angle [deg]')
    ax.set_title('Angle Estimation — Open-Loop Integrator (DC drift visible!)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # (d) Angle errors
    ax = axes[3]
    rms_af = np.sqrt(np.mean(err_af[n_ss:]**2))
    rms_ol = np.sqrt(np.mean(err_ol[n_ss:]**2))
    ax.plot(t, err_af, '#d62728', linewidth=0.8, label=f'AF error (RMS={rms_af:.1f}°)')
    ax.plot(t, err_ol, '#2ca02c', linewidth=0.8, alpha=0.7, label=f'OL error (RMS={rms_ol:.1f}°)')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.set_ylabel('Angle Error [deg]')
    ax.set_xlabel('Time [s]')
    ax.set_title('Angle Estimation Error Comparison')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)

    fig1.tight_layout()
    fig1.savefig(f'{save_path}_angles.png', dpi=200, bbox_inches='tight')
    print(f'  Saved {save_path}_angles.png')

    # ================================================================
    # Figure 2: Speed estimation — LPF vs 4th-order ESO (4 subplots)
    # ================================================================
    fig2, axes2 = plt.subplots(4, 1, dpi=150, facecolor='w', figsize=(13, 14), sharex=True)

    omega_ob = results['eso_omega_ob']
    obs_name = "4th-order ESO" if results.get('speed_observer', 'eso') == 'eso' else "Natural Speed Observer"

    # (a) Speed comparison: true, LPF, ESO/NSO
    ax = axes2[0]
    ax.plot(t, results['omega_true'], '#1f77b4', linewidth=1.0, alpha=0.6, label='True speed (encoder)')
    ax.plot(t, results['omega_af_lpf'], '#d62728', linewidth=0.9,
            label=r'AF + LPF ($\tau$=2ms, 1st-order)')
    ax.plot(t, results['omega_af_eso'], '#9467bd', linewidth=1.0,
            label=f'AF + {obs_name} ($\\omega_{{ob}}$={omega_ob:.0f} rad/s)')
    ax.set_ylabel('Speed [rpm]')
    ax.set_title(
        f'Speed Estimation Comparison: LPF (1st-order) vs {obs_name}\n'
        '(both methods use the same Active Flux angle estimate as input)',
        fontsize=11, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # (b) Speed error
    ax = axes2[1]
    spd_err_lpf = results['omega_af_lpf'] - results['omega_true']
    spd_err_eso = results['omega_af_eso'] - results['omega_true']
    rms_spd_lpf = np.sqrt(np.mean(spd_err_lpf[n_ss:]**2))
    rms_spd_eso = np.sqrt(np.mean(spd_err_eso[n_ss:]**2))
    ax.plot(t, spd_err_lpf, '#d62728', linewidth=0.7, alpha=0.7,
            label=f'LPF error (RMS={rms_spd_lpf:.2f} rpm)')
    ax.plot(t, spd_err_eso, '#9467bd', linewidth=0.7,
            label=f'{obs_name} error (RMS={rms_spd_eso:.2f} rpm)')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.set_ylabel('Speed Error [rpm]')
    ax.set_title('Speed Estimation Error')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)

    # (c) ESO load-torque estimation
    ax = axes2[2]
    # ESO disturbance state x[2] represents acceleration perturbation
    # Convert: load_torque ≈ -x[2] (opposite sign, since load opposes motion)
    ax.plot(t, results['load_true'], '#1f77b4', linewidth=1.5, alpha=0.7, label='True load torque [N·m]')
    # x[2] * npp / Js gives acceleration; multiply back to get torque equivalent
    n_pp_val = results['n_pp']
    Js_val   = results['Js']
    load_est_Nm = -results['load_eso']  # opposite sign convention
    ax.plot(t, load_est_Nm, '#ff7f0e', linewidth=1.0,
            label=f'{obs_name} disturbance est.')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.set_ylabel('Load Torque [N·m]')
    ax.set_title(
        f'{obs_name} Load Torque Estimation (ω_ob={omega_ob:.0f} rad/s)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # (d) Observer Output Angle Error (only meaningful for ESO, but plot for both)
    ax = axes2[3]
    eso_out_err = wrap_angle(results['theta_af'] - results['theta_af_eso'])
    ax.plot(t, np.degrees(eso_out_err), '#2ca02c', linewidth=1.0, label='Observer Output Angle Error (deg)')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.set_ylabel('Output Error [deg]')
    ax.set_xlabel('Time [s]')
    ax.set_title(r'Observer Angle Output Error ($\theta_{meas} - \theta_{est}$) — 0 for NSO')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    obs_suffix = 'speed_nso' if results.get('speed_observer', 'eso') == 'nso' else 'speed_eso'
    fig2.tight_layout()
    fig2.savefig(f'{save_path}_{obs_suffix}.png', dpi=200, bbox_inches='tight')
    print(f'  Saved {save_path}_{obs_suffix}.png')

    # ================================================================
    # Figure 3: Flux trajectories (α-β plane)
    # ================================================================
    fig3, (ax1, ax2) = plt.subplots(1, 2, dpi=150, facecolor='w', figsize=(14, 6))

    theta_circle = np.linspace(0, 2 * np.pi, 200)
    ideal_x = KE * np.cos(theta_circle)
    ideal_y = KE * np.sin(theta_circle)
    n_start = len(t) // 5

    ax1.plot(ideal_x, ideal_y, '#1f77b4', linewidth=1.5, linestyle='--', alpha=0.6, label=f'Ideal |ψ|={KE} Wb')
    ax1.plot(results['psi_af_ab'][n_start:, 0], results['psi_af_ab'][n_start:, 1],
             '#d62728', linewidth=0.3, alpha=0.7, label='AF observer')
    ax1.set_xlabel(r'$\psi_{AF,\alpha}$ [Wb]')
    ax1.set_ylabel(r'$\psi_{AF,\beta}$ [Wb]')
    ax1.set_title('Active Flux Observer\n(circle maintained by PI correction)', fontsize=10)
    ax1.set_aspect('equal')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.3)

    ax2.plot(ideal_x, ideal_y, '#1f77b4', linewidth=1.5, linestyle='--', alpha=0.6, label=f'Ideal |ψ|={KE} Wb')
    ax2.plot(results['psi_ol_ab'][n_start:, 0], results['psi_ol_ab'][n_start:, 1],
             '#2ca02c', linewidth=0.3, alpha=0.7, label='Open-loop int.')
    ax2.set_xlabel(r'$\psi_{AF,\alpha}$ [Wb]')
    ax2.set_ylabel(r'$\psi_{AF,\beta}$ [Wb]')
    ax2.set_title('Open-Loop Integrator\n(DC drift causes center offset)', fontsize=10)
    ax2.set_aspect('equal')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig3.suptitle(r'Active Flux Trajectory in $\alpha$-$\beta$ Plane', fontsize=13, fontweight='bold', y=1.02)
    fig3.tight_layout()
    fig3.savefig(f'{save_path}_flux_trajectory.png', dpi=200, bbox_inches='tight')
    print(f'  Saved {save_path}_flux_trajectory.png')

    # ================================================================
    # Print summary
    # ================================================================
    max_af = np.max(np.abs(err_af[n_ss:]))
    max_ol = np.max(np.abs(err_ol[n_ss:]))

    print('\n' + '=' * 70)
    print('  Sensorless Estimation Accuracy Comparison')
    print('=' * 70)
    print(f'  --- Angle Estimation ---')
    print(f'  {"Metric":<35s} {"Active Flux":>13s} {"Open-Loop":>13s}')
    print(f'  {"-"*35:<35s} {"-"*13:>13s} {"-"*13:>13s}')
    print(f'  {"RMS angle error [deg]":<35s} {rms_af:>13.2f} {rms_ol:>13.2f}')
    print(f'  {"Max angle error [deg]":<35s} {max_af:>13.2f} {max_ol:>13.2f}')
    print()
    obs_col = 'NSO (3rd)' if results.get('speed_observer', 'eso') == 'nso' else 'ESO (4th)'
    print(f'  --- Speed Estimation (from Active Flux angle) ---')
    print(f'  {"Metric":<35s} {"LPF (1st)":>13s} {obs_col:>13s}')
    print(f'  {"-"*35:<35s} {"-"*13:>13s} {"-"*13:>13s}')
    print(f'  {"RMS speed error [rpm]":<35s} {rms_spd_lpf:>13.2f} {rms_spd_eso:>13.2f}')
    max_spd_lpf = np.max(np.abs(spd_err_lpf[n_ss:]))
    max_spd_eso = np.max(np.abs(spd_err_eso[n_ss:]))
    print(f'  {"Max speed error [rpm]":<35s} {max_spd_lpf:>13.2f} {max_spd_eso:>13.2f}')
    print('=' * 70)

    return fig1, fig2, fig3


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':
    print('=' * 70)
    print('  无传感器控制 Demo — Active Flux 估计 + ESO 速度观测')
    print('  Sensorless Demo — Active Flux + 4th-order ESO Speed Observer')
    print('=' * 70)

    # ---------- 电机参数（小电感电机，与 ep6 一致） ----------
    d = {
        'CL_TS': 1e-4,
        'VL_EXE_PER_CL_EXE': 5,
        'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
        'TIME_SLICE': 0.1,
        'NUMBER_OF_SLICES': 20,  # 2.0 s total
        'init_npp': 22,
        'init_IN': 1.3 * 6 / 1.414,
        'init_R': 0.035,
        'init_Ld': 1 * 0.036e-3,
        'init_Lq': 1 * 0.036e-3,
        'init_KE': 0.0125,
        'init_Rreq': 0.0,
        'init_Js': 0.44e-4,
        'DC_BUS_VOLTAGE': 5,
        'CTRL.bool_apply_speed_closed_loop_control': True,
        'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
        'CTRL.bool_apply_sweeping_frequency_excitation': False,
        'CTRL.bool_overwrite_speed_commands': True,
        'CTRL.bool_zero_id_control': True,
        'FOC_delta': 15,
        'FOC_desired_VLBW_HZ': 120,
        'FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False': 10,
        'CL_SERIES_KP': None,
        'CL_SERIES_KI': None,
        'VL_SERIES_KP': None,
        'VL_SERIES_KI': None,
        'VL_LIMIT_OVERLOAD_FACTOR': 3.0,
        'disp.Kp': 0.0,
        'disp.Ki': 0.0,
        'disp.Kd': 0.0,
        'disp.tau': 0.0,
        'disp.OutLimit': 0.0,
        'disp.IntLimit': 0.0,
    }

    # ====================================================
    # Parse command-line switches
    # ====================================================
    import argparse
    parser = argparse.ArgumentParser(description='Sensorless Active Flux + ESO Demo')
    parser.add_argument('--use-sensorless-speed', action='store_true', default=False,
                        help='Use ESO-estimated speed for FOC speed loop (and AF angle for Park)')
    parser.add_argument('--use-sensorless-torque-ff', action='store_true', default=False,
                        help='Use ESO load torque estimate for feedforward into iq command')
    parser.add_argument('--eso-omega-ob', type=float, default=200.0,
                        help='Dynamics observer bandwidth [rad/s] (default: 200)')
    parser.add_argument('--speed-observer', type=str, choices=['eso', 'nso'], default='eso',
                        help='Choose intermediate dynamics observer: eso (4th-order) or nso (Natural Speed Observer)')
    parser.add_argument('--clbw', type=float, default=1000.0,
                        help='Current loop bandwidth [Hz]')
    parser.add_argument('--af-kp', type=float, default=500.0,
                        help='Active flux error correction Kp')
    parser.add_argument('--af-ki', type=float, default=5000.0,
                        help='Active flux error correction Ki')
    parser.add_argument('--zeta', type=float, default=15.0,
                        help='Speed loop damping ratio')
    parser.add_argument('--pure-p-current', action='store_true', default=False,
                        help='Use pure P controller for current loops instead of PI')
    args = parser.parse_args()

    # ====================================================
    # Run sensorless demo
    # ====================================================
    results = run_sensorless_demo(
        d,
        zeta=args.zeta,
        CLBW_Hz=args.clbw,
        af_Kp=args.af_kp,
        af_Ki=args.af_ki,
        R_mismatch_factor=1.5,
        voltage_offset_alpha=0.02,
        voltage_offset_beta=-0.015,
        eso_omega_ob=args.eso_omega_ob,
        speed_observer=args.speed_observer,
        use_sensorless_speed=args.use_sensorless_speed,
        use_sensorless_torque_ff=args.use_sensorless_torque_ff,
        pure_p_current=args.pure_p_current,
        verbose=True,
    )

    # ====================================================
    # Plot results
    # ====================================================
    suffix = ""
    if args.use_sensorless_speed:
        suffix += "_spd"
    if args.use_sensorless_torque_ff:
        suffix += "_tqff"
        
    base_name = f'fig_sensorless_demo{suffix}'
    fig1, fig2, fig3 = plot_sensorless_results(results, save_path=base_name)

    plt.close('all')
    print('\n--- 所有图已保存 ---')
    print(f'  {base_name}_angles.png          (角度估计对比)')
    print(f'  {base_name}_speed_eso.png       (转速: LPF vs 4th-order ESO)')
    print(f'  {base_name}_flux_trajectory.png  (αβ 磁链轨迹)')
