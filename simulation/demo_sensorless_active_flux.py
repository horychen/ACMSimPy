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

3. **对比图**
   - 角度估计 vs 编码器真值
   - 转速估计 vs 真实转速
   - 磁链轨迹（αβ 平面）
   - 角度误差时域波形

原理说明
--------
定子磁链电压模型：
    dψ_s/dt = u_s - R·i_s

Active flux（有功磁链）：
    ψ_AF = ψ_s - Lq·i_s

对于 SPMSM (Ld = Lq)：
    ψ_AF ≈ ψ_PM（永磁磁链向量）

转子位置由 ψ_AF 的方向确定：
    θ_est = arctan2(ψ_AF_β, ψ_AF_α)

现实中的问题
------------
开环积分器存在 DC 漂移问题，主要来源：
  - 电阻 R 的估计误差（温度变化等）
  - 电压测量偏置（逆变器死区、AD 偏移等）
  - 积分初值误差
本 demo 故意引入 R 估计偏差和电压偏置来模拟真实情况。
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


# ======================================================================
# Sensorless simulation (Active Flux observer + open-loop integrator)
# ======================================================================
def run_sensorless_demo(d,
                        zeta=15, CLBW_Hz=1000,
                        af_Kp=500.0, af_Ki=5000.0,
                        R_mismatch_factor=1.5,
                        voltage_offset_alpha=0.02,
                        voltage_offset_beta=-0.015,
                        verbose=True):
    """
    Run a sensorless control demo simulation.

    Two flux observers run in parallel alongside the true encoder:
      (A) Active Flux observer with PI correction
      (B) Pure open-loop integrator (no correction)

    The motor is still controlled using the *true* encoder (sensored FOC)
    so that we can objectively evaluate the observation accuracy of both
    methods without instability from feeding back a bad estimate.

    Parameters
    ----------
    d : dict
        Motor / controller parameter dictionary.
    zeta, CLBW_Hz : float
        Speed-loop damping and current-loop bandwidth for PI tuning.
    af_Kp, af_Ki : float
        PI gains for the Active Flux observer correction term.
    R_mismatch_factor : float
        Observer uses R_obs = R_true * factor. >1 means overestimated R.
        This deliberately introduces a common real-world error source.
    voltage_offset_alpha, voltage_offset_beta : float
        Simulated voltage measurement offsets [V] injected into the
        observers' voltage input. Represents inverter deadzone, ADC
        offset, etc.
    verbose : bool
        Print progress.

    Returns
    -------
    results : dict
        Time traces and observer states for plotting.
    """
    dd = copy.deepcopy(d)

    # ------- PI tuning -------
    R_true = dd['init_R']            # true resistance (used by FOC)
    R_obs  = R_true * R_mismatch_factor  # observer's (wrong) estimate
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

    if verbose:
        print(f'  R_true = {R_true:.4f} Ω, R_obs = {R_obs:.4f} Ω (factor={R_mismatch_factor})')
        print(f'  Voltage offsets: α={voltage_offset_alpha:.4f} V, β={voltage_offset_beta:.4f} V')
        print(f'  KE = {KE} Wb, L = {L*1e3:.4f} mH')

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

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=dd['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])

    # Current & speed regulators
    if dd.get('CTRL.bool_apply_decoupling_voltages_to_current_regulation', False):
        local_Ki_factor = 1.0
    else:
        local_Ki_factor = dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False', 10)
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
    #  Simulation parameters
    # ============================================================
    total_time = dd['TIME_SLICE'] * dd['NUMBER_OF_SLICES']
    controller_down_sampling = int(CL_TS / MACHINE_TS)

    # Storage arrays (at CL_TS rate to save memory)
    n_ctrl = int(total_time / CL_TS) + 100  # small overflow buffer
    t_trace       = np.zeros(n_ctrl)
    theta_true    = np.zeros(n_ctrl)
    omega_true    = np.zeros(n_ctrl)
    theta_af      = np.zeros(n_ctrl)
    omega_af      = np.zeros(n_ctrl)
    psi_af_ab     = np.zeros((n_ctrl, 2))
    theta_ol      = np.zeros(n_ctrl)
    omega_ol      = np.zeros(n_ctrl)
    psi_ol_ab     = np.zeros((n_ctrl, 2))
    cmd_rpm_trace = np.zeros(n_ctrl)

    # Speed estimation via differentiation of angle (with LPF)
    af_omega_est  = 0.0
    ol_omega_est  = 0.0
    af_theta_prev = 0.0
    ol_theta_prev = 0.0
    lpf_tau   = 0.002  # LPF time constant for speed estimation [s]
    lpf_alpha = CL_TS / (lpf_tau + CL_TS)

    # ============================================================
    #  Main simulation loop
    # ============================================================
    if verbose:
        print('Starting sensorless demo simulation...')
    t_start_sim = _time.time()

    CTRL.cmd_rpm = 0.0
    ctrl_idx = 0
    n_slices = dd['NUMBER_OF_SLICES']

    for slice_idx in range(n_slices):
        t0_slice = slice_idx * dd['TIME_SLICE']
        t_mid = t0_slice + dd['TIME_SLICE'] * 0.5

        # ---- Speed command profile ----
        # 0.0 – 0.3 s : ramp to +200 rpm
        # 0.3 – 0.6 s : hold  +200 rpm
        # 0.6 – 0.8 s : step load applied
        # 0.8 – 1.0 s : hold  +200 rpm, load removed
        # 1.0 – 1.3 s : speed reversal to -200 rpm
        # 1.3 – 2.0 s : hold  -200 rpm
        if t_mid < 0.3:
            CTRL.cmd_rpm = 200 * min(t_mid / 0.3, 1.0)
        elif t_mid < 0.6:
            CTRL.cmd_rpm = 200
            ACM.TLoad = 0.0
        elif t_mid < 0.8:
            CTRL.cmd_rpm = 200
            ACM.TLoad = 0.15
        elif t_mid < 1.0:
            CTRL.cmd_rpm = 200
            ACM.TLoad = 0.0
        elif t_mid < 1.3:
            CTRL.cmd_rpm = -200
        else:
            CTRL.cmd_rpm = -200

        machine_times, watch_data = ACMSimPyIncremental(
            t0=t0_slice, TIME=dd['TIME_SLICE'],
            ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed
        )

        # watch_data layout (from ep6):
        #  0: ACM.theta_d (mod 2π)
        #  1: ACM.omega_r_mech in rpm
        # 12: CTRL.cmd_rpm
        # 28: CTRL.cmd_uab[0]
        # 29: CTRL.cmd_uab[1]
        #  6: CTRL.iab[0]
        #  7: CTRL.iab[1]

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

            # Measured currents and commanded voltages
            i_alpha = watch_data[6][k]
            i_beta  = watch_data[7][k]
            # Add deliberate voltage offsets (simulating real-world imperfections)
            u_alpha = watch_data[28][k] + voltage_offset_alpha
            u_beta  = watch_data[29][k] + voltage_offset_beta

            # ============================================================
            # (A) Active Flux Observer with PI amplitude correction
            # ============================================================
            # Step 1: Compute current active flux estimate
            af_active_alpha = af_psi_s[0] - L * i_alpha
            af_active_beta  = af_psi_s[1] - L * i_beta
            af_active_amp   = np.sqrt(af_active_alpha**2 + af_active_beta**2)

            # Step 2: PI correction on amplitude error
            #   error = |ψ_AF_cmd| - |ψ_AF_est|
            #   The correction is applied in the direction of the flux vector
            psi_af_error = KE - af_active_amp
            af_corr_int[0] += af_Ki * psi_af_error * CL_TS
            af_corr_int[1] += af_Ki * psi_af_error * CL_TS

            # Correction direction: along ψ_AF unit vector
            if af_active_amp > 1e-10:
                af_unit_alpha = af_active_alpha / af_active_amp
                af_unit_beta  = af_active_beta  / af_active_amp
            else:
                af_unit_alpha = 1.0
                af_unit_beta  = 0.0

            corr_alpha = (af_Kp * psi_af_error + af_corr_int[0]) * af_unit_alpha
            corr_beta  = (af_Kp * psi_af_error + af_corr_int[1]) * af_unit_beta

            # Step 3: Euler integration of stator flux with PI correction
            #   Note: uses R_obs (which may differ from R_true)
            af_psi_s[0] += CL_TS * (u_alpha - R_obs * i_alpha + corr_alpha)
            af_psi_s[1] += CL_TS * (u_beta  - R_obs * i_beta  + corr_beta)

            # Step 4: Extract angle from active flux
            af_active_alpha = af_psi_s[0] - L * i_alpha
            af_active_beta  = af_psi_s[1] - L * i_beta
            theta_af_now = np.arctan2(af_active_beta, af_active_alpha)

            # Step 5: Speed estimation via d(theta)/dt with LPF
            d_theta_af = angle_error(np.array([theta_af_now]), np.array([af_theta_prev]))[0]
            raw_omega_af = d_theta_af / CL_TS  # elec. rad/s
            af_omega_est = af_omega_est * (1 - lpf_alpha) + raw_omega_af * lpf_alpha
            af_theta_prev = theta_af_now

            theta_af[ctrl_idx]  = theta_af_now
            omega_af[ctrl_idx]  = af_omega_est / (2 * np.pi * n_pp) * 60  # to rpm
            psi_af_ab[ctrl_idx] = [af_active_alpha, af_active_beta]

            # ============================================================
            # (B) Open-loop pure integrator (no correction)
            # ============================================================
            # Same voltage model, same wrong R_obs, no correction at all
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
            omega_ol[ctrl_idx]  = ol_omega_est / (2 * np.pi * n_pp) * 60  # to rpm
            psi_ol_ab[ctrl_idx] = [ol_active_alpha, ol_active_beta]

            ctrl_idx += 1

    elapsed = _time.time() - t_start_sim
    if verbose:
        print(f'  Simulation completed in {elapsed:.2f} s, {ctrl_idx} control steps.')

    # Trim arrays
    n = ctrl_idx
    results = {
        't':          t_trace[:n],
        'theta_true': theta_true[:n],
        'omega_true': omega_true[:n],
        'cmd_rpm':    cmd_rpm_trace[:n],
        'theta_af':   theta_af[:n],
        'omega_af':   omega_af[:n],
        'psi_af_ab':  psi_af_ab[:n],
        'theta_ol':   theta_ol[:n],
        'omega_ol':   omega_ol[:n],
        'psi_ol_ab':  psi_ol_ab[:n],
        'n_pp':       n_pp,
        'KE':         KE,
        'R_true':     R_true,
        'R_obs':      R_obs,
        'v_offset_a': voltage_offset_alpha,
        'v_offset_b': voltage_offset_beta,
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

    # Pre-compute angle errors
    err_af = np.degrees(angle_error(results['theta_af'], results['theta_true']))
    err_ol = np.degrees(angle_error(results['theta_ol'], results['theta_true']))

    n_ss = len(t) // 5  # skip initial 20% for RMS calculation

    # ================================================================
    # Figure 1: Time-domain comparison (5 subplots)
    # ================================================================
    fig1, axes = plt.subplots(5, 1, dpi=150, facecolor='w', figsize=(13, 16), sharex=True)

    # --- (a) Speed command and true speed ---
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

    # --- (b) Speed estimation comparison ---
    ax = axes[1]
    ax.plot(t, results['omega_true'], '#1f77b4', linewidth=0.8, alpha=0.5, label='True speed')
    ax.plot(t, results['omega_af'], '#d62728', linewidth=1.0, label='Active Flux est.')
    ax.plot(t, results['omega_ol'], '#2ca02c', linewidth=0.8, alpha=0.7, label='Open-loop est.')
    ax.set_ylabel('Estimated Speed [rpm]')
    ax.set_title('Speed Estimation Comparison')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # --- (c) Angle estimation — Active Flux ---
    ax = axes[2]
    theta_true_wrapped = np.degrees(wrap_angle(results['theta_true']))
    theta_af_wrapped   = np.degrees(wrap_angle(results['theta_af']))
    ax.plot(t, theta_true_wrapped, '#1f77b4', linewidth=0.6, alpha=0.5, label=r'$\theta_{true}$ (encoder)')
    ax.plot(t, theta_af_wrapped,   '#d62728', linewidth=0.6, label=r'$\theta_{AF}$ (Active Flux)')
    ax.set_ylabel('Angle [deg]')
    ax.set_title('Angle Estimation — Active Flux Observer (with PI correction)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # --- (d) Angle estimation — Open-loop ---
    ax = axes[3]
    theta_ol_wrapped = np.degrees(wrap_angle(results['theta_ol']))
    ax.plot(t, theta_true_wrapped, '#1f77b4', linewidth=0.6, alpha=0.5, label=r'$\theta_{true}$ (encoder)')
    ax.plot(t, theta_ol_wrapped,   '#2ca02c', linewidth=0.6, label=r'$\theta_{OL}$ (Open-loop)')
    ax.set_ylabel('Angle [deg]')
    ax.set_title('Angle Estimation — Open-Loop Integrator (DC drift visible!)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # --- (e) Angle errors ---
    ax = axes[4]
    rms_af = np.sqrt(np.mean(err_af[n_ss:]**2))
    rms_ol = np.sqrt(np.mean(err_ol[n_ss:]**2))
    ax.plot(t, err_af, '#d62728', linewidth=0.8,
            label=f'AF error (RMS={rms_af:.1f}°)')
    ax.plot(t, err_ol, '#2ca02c', linewidth=0.8, alpha=0.7,
            label=f'OL error (RMS={rms_ol:.1f}°)')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.set_ylabel('Angle Error [deg]')
    ax.set_xlabel('Time [s]')
    ax.set_title('Angle Estimation Error Comparison')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)

    fig1.tight_layout()
    fig1.savefig(f'{save_path}_angles_speed.png', dpi=200, bbox_inches='tight')
    print(f'  Saved {save_path}_angles_speed.png')

    # ================================================================
    # Figure 2: Flux trajectories (α-β plane)
    # ================================================================
    fig2, (ax1, ax2) = plt.subplots(1, 2, dpi=150, facecolor='w', figsize=(14, 6))

    theta_circle = np.linspace(0, 2 * np.pi, 200)
    ideal_x = KE * np.cos(theta_circle)
    ideal_y = KE * np.sin(theta_circle)

    n_start = len(t) // 5

    # -- Active Flux observer --
    ax1.plot(ideal_x, ideal_y, '#1f77b4', linewidth=1.5, linestyle='--', alpha=0.6, label=f'Ideal |ψ|={KE} Wb')
    ax1.plot(results['psi_af_ab'][n_start:, 0], results['psi_af_ab'][n_start:, 1],
             '#d62728', linewidth=0.3, alpha=0.7, label='AF observer')
    ax1.set_xlabel(r'$\psi_{AF,\alpha}$ [Wb]')
    ax1.set_ylabel(r'$\psi_{AF,\beta}$ [Wb]')
    ax1.set_title('Active Flux Observer\n(circle maintained by PI correction)', fontsize=10)
    ax1.set_aspect('equal')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.grid(True, alpha=0.3)

    # -- Open-loop integrator --
    ax2.plot(ideal_x, ideal_y, '#1f77b4', linewidth=1.5, linestyle='--', alpha=0.6, label=f'Ideal |ψ|={KE} Wb')
    ax2.plot(results['psi_ol_ab'][n_start:, 0], results['psi_ol_ab'][n_start:, 1],
             '#2ca02c', linewidth=0.3, alpha=0.7, label='Open-loop int.')
    ax2.set_xlabel(r'$\psi_{AF,\alpha}$ [Wb]')
    ax2.set_ylabel(r'$\psi_{AF,\beta}$ [Wb]')
    ax2.set_title('Open-Loop Integrator\n(DC drift causes center offset)', fontsize=10)
    ax2.set_aspect('equal')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig2.suptitle(r'Active Flux Trajectory in $\alpha$-$\beta$ Plane', fontsize=13, fontweight='bold', y=1.02)
    fig2.tight_layout()
    fig2.savefig(f'{save_path}_flux_trajectory.png', dpi=200, bbox_inches='tight')
    print(f'  Saved {save_path}_flux_trajectory.png')

    # ================================================================
    # Print summary
    # ================================================================
    max_af = np.max(np.abs(err_af[n_ss:]))
    max_ol = np.max(np.abs(err_ol[n_ss:]))

    print('\n' + '=' * 65)
    print('  Sensorless Estimation Accuracy Comparison')
    print('=' * 65)
    print(f'  {"Metric":<35s} {"Active Flux":>13s} {"Open-Loop":>13s}')
    print(f'  {"-"*35:<35s} {"-"*13:>13s} {"-"*13:>13s}')
    print(f'  {"RMS angle error [deg]":<35s} {rms_af:>13.2f} {rms_ol:>13.2f}')
    print(f'  {"Max angle error [deg]":<35s} {max_af:>13.2f} {max_ol:>13.2f}')
    print(f'  {"R mismatch factor":<35s} {"—":>13s} {"—":>13s}')
    print(f'  {"  (observer R / true R)":<35s} {results["R_obs"]/results["R_true"]:>13.2f} {results["R_obs"]/results["R_true"]:>13.2f}')
    print(f'  {"Voltage offset α [V]":<35s} {results["v_offset_a"]:>13.4f} {results["v_offset_a"]:>13.4f}')
    print(f'  {"Voltage offset β [V]":<35s} {results["v_offset_b"]:>13.4f} {results["v_offset_b"]:>13.4f}')
    print('=' * 65)

    return fig1, fig2


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':
    print('=' * 65)
    print('  无传感器控制 Demo — Active Flux 估计 vs 开环积分')
    print('  Sensorless Demo — Active Flux Observer vs Open-Loop Integrator')
    print('=' * 65)

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
    # Run sensorless demo
    # ====================================================
    results = run_sensorless_demo(
        d,
        zeta=15,
        CLBW_Hz=1000,
        af_Kp=500.0,   # Active Flux observer proportional gain
        af_Ki=5000.0,  # Active Flux observer integral gain
        # --- Deliberately introduced imperfections ---
        R_mismatch_factor=1.5,        # observer overestimates R by 50%
        voltage_offset_alpha=0.02,    # 20 mV offset on α channel
        voltage_offset_beta=-0.015,   # -15 mV offset on β channel
        verbose=True,
    )

    # ====================================================
    # Plot results
    # ====================================================
    fig1, fig2 = plot_sensorless_results(results, save_path='fig_sensorless_demo')

    plt.close('all')
    print('\n--- 所有图已保存 ---')
    print('  fig_sensorless_demo_angles_speed.png')
    print('  fig_sensorless_demo_flux_trajectory.png')
