#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
=======================================================================
  Break NSO-AF Coupling: αβ-frame Speed Estimation for servo
=======================================================================
Key finding from previous diagnostic:
  - S2c (encoder speed + AF angle) = 0.77° ★ PERFECT
  - S3 (NSO speed + AF angle) = 104° ❌ FAIL
  - Root cause: NSO uses AF angle for Park → positive feedback loop

Strategy: Replace dq-frame NSO with αβ-frame speed estimators:
  1. LPF dθ/dt (1st/2nd order)  
  2. 2nd-order PLL on AF angle
  3. αβ-frame Luenberger (no Park transform)
  
Also sweep speed-loop zeta (1~10) since lower zeta reduces coupling gain.
"""
import matplotlib; matplotlib.use('Agg')
from pylab import np, plt
import copy, time as _time

from tutorials_ep6_svpwm import (The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental)
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator
from demo_sensorless_active_flux import angle_error, angle_diff_scalar
from eval_staged_tuning import get_motor_preset


class PLLSpeedTracker:
    """2nd-order PLL to track AF angle → estimate speed. No Park transform needed."""
    def __init__(self, bw, dt, npp):
        self.dt = dt; self.npp = npp
        wn = 2*np.pi*bw
        self.kp = 2*wn     # proportional (rad/s per rad error)
        self.ki = wn**2     # integral (rad/s^2 per rad error)
        self.theta_est = 0.0
        self.omega_est = 0.0  # electrical rad/s
        self.integrator = 0.0
        
    def step(self, theta_af):
        """Track AF angle, output speed in elec rad/s."""
        err = angle_diff_scalar(theta_af, self.theta_est)
        self.integrator += self.ki * err * self.dt
        self.omega_est = self.kp * err + self.integrator
        self.theta_est += self.omega_est * self.dt
        self.theta_est = np.arctan2(np.sin(self.theta_est), np.cos(self.theta_est))
        return self.omega_est
    
    @property
    def omega_rpm(self):
        return self.omega_est / (2*np.pi*self.npp) * 60


class CascadedLPF:
    """2nd-order LPF via two cascaded 1st-order LPFs."""
    def __init__(self, tau, dt):
        self.alpha = dt / (tau + dt)
        self.y1 = 0.0
        self.y2 = 0.0
    
    def step(self, x):
        self.y1 = self.y1*(1-self.alpha) + x*self.alpha
        self.y2 = self.y2*(1-self.alpha) + self.y1*self.alpha
        return self.y2


def run_ab_speed_s3(d, CLBW_Hz, zeta, af_Kp, af_Ki, R_mis, v_off,
                    cmd_rpm, load_step,
                    speed_method='pll', speed_bw=50.0, lpf_tau=0.005,
                    verbose=True):
    """
    Stage 3 sensorless with αβ-frame speed estimation (no NSO Park coupling).
    
    speed_method:
      'lpf1' - 1st order LPF on dθ/dt
      'lpf2' - 2nd order (cascaded) LPF on dθ/dt  
      'pll'  - 2nd order PLL tracking AF angle
    """
    dd = copy.deepcopy(d)
    R_true = dd['init_R']; R_obs = R_true * R_mis
    L = dd['init_Lq']; J_s = dd['init_Js']; n_pp = dd['init_npp']
    KE = dd['init_KE']; CL_TS = dd['CL_TS']
    VL_TS = CL_TS * dd['VL_EXE_PER_CL_EXE']

    cKp, cKi = get_coeffs_dc_motor_current_regulator(R_true, L, CLBW_Hz)
    sKp, sKi = get_coeffs_dc_motor_SPEED_regulator(J_s, n_pp, KE, zeta, cKp/L)
    dd['CL_SERIES_KP'] = cKp; dd['CL_SERIES_KI'] = cKi
    dd['VL_SERIES_KP'] = sKp; dd['VL_SERIES_KI'] = sKi

    CTRL = The_Motor_Controller(CL_TS=CL_TS, VL_TS=VL_TS, init_npp=n_pp,
        init_IN=dd['init_IN'], init_R=R_true, init_Ld=dd['init_Ld'],
        init_Lq=L, init_KE=KE, init_Rreq=0.0, init_Js=J_s,
        DC_BUS_VOLTAGE=dd['DC_BUS_VOLTAGE'])
    CTRL.bool_apply_decoupling_voltages_to_current_regulation = False
    CTRL.bool_apply_sweeping_frequency_excitation = False
    CTRL.bool_overwrite_speed_commands = True
    CTRL.bool_zero_id_control = True
    CTRL.bool_apply_speed_closed_loop_control = True

    # Full sensorless: AF angle for Park, αβ speed for speed loop
    CTRL.index_separate_speed_estimation = 1
    CTRL.ell1 = CTRL.ell2 = CTRL.ell3 = CTRL.ell4 = 0.0
    CTRL.bool_use_sensorless_theta = 1
    CTRL.use_disturbance_feedforward_rejection = 0

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1)
    ki_f = dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False', 10)
    vl_lim = dd.get('VL_LIMIT_OVERLOAD_FACTOR', 10.0)
    Vm = dd['DC_BUS_VOLTAGE']/1.732; Im = vl_lim*1.414*dd['init_IN']
    reg_id = The_PID_Regulator(cKp, cKp*cKi*ki_f, 0,0, Vm, Vm, CL_TS)
    reg_iq = The_PID_Regulator(cKp, cKp*cKi*ki_f, 0,0, Vm, Vm, CL_TS)
    reg_spd = The_PID_Regulator(sKp, sKp*sKi, 0,0, Im, Im, VL_TS)

    # AF state
    af_psi_s = np.array([KE, 0.0]); af_corr_int = np.zeros(2)
    
    # Speed estimator (no Park transform!)
    if speed_method == 'pll':
        pll = PLLSpeedTracker(bw=speed_bw, dt=CL_TS, npp=n_pp)
    elif speed_method == 'lpf2':
        lpf = CascadedLPF(tau=lpf_tau, dt=CL_TS)
    
    lpf1_alpha = CL_TS / (lpf_tau + CL_TS)
    omega_lpf1 = 0.0
    theta_prev = 0.0

    total_time = dd['TIME_SLICE'] * dd['NUMBER_OF_SLICES']
    n_ctrl = int(total_time / CL_TS) + 100
    t_tr = np.zeros(n_ctrl); th_true = np.zeros(n_ctrl); w_true = np.zeros(n_ctrl)
    th_af = np.zeros(n_ctrl); w_est = np.zeros(n_ctrl); cmd_tr = np.zeros(n_ctrl)

    CTRL.cmd_rpm = 0.0; ci = 0

    for sl in range(int(total_time / CL_TS)):
        t0 = sl * CL_TS; tm = t0 + CL_TS*0.5
        if tm<0.3: CTRL.cmd_rpm = cmd_rpm*min(tm/0.3,1); ACM.TLoad = 0
        elif tm<0.6: CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0
        elif tm<0.8: CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = load_step
        elif tm<1.0: CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0
        elif tm<1.3: CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0
        else: CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0

        if ci > 0:
            CTRL.vartheta_d = th_af[ci-1]
            # Inject αβ-frame speed estimate (no dq Park!)
            CTRL.xS[1] = omega_elec_est

        mt, wd = ACMSimPyIncremental(t0=t0, TIME=CL_TS, ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_spd)

        if ci >= n_ctrl: break
        t_tr[ci] = mt[0]; th_true[ci] = wd[0][0]; w_true[ci] = wd[1][0]; cmd_tr[ci] = wd[12][0]
        ia = wd[6][0]; ib = wd[7][0]; ua = wd[28][0]+v_off[0]; ub = wd[29][0]+v_off[1]

        # AF observer (unchanged)
        af_a = af_psi_s[0]-L*ia; af_b = af_psi_s[1]-L*ib
        af_amp = np.sqrt(af_a**2+af_b**2)
        psi_err = KE-af_amp; af_corr_int += af_Ki*psi_err*CL_TS
        if af_amp>1e-10: un_a,un_b = af_a/af_amp, af_b/af_amp
        else: un_a,un_b = 1.0, 0.0
        ca = (af_Kp*psi_err+af_corr_int[0])*un_a; cb = (af_Kp*psi_err+af_corr_int[1])*un_b
        af_psi_s[0] += CL_TS*(ua-R_obs*ia+ca); af_psi_s[1] += CL_TS*(ub-R_obs*ib+cb)
        af_a2=af_psi_s[0]-L*ia; af_b2=af_psi_s[1]-L*ib
        theta_af_now = np.arctan2(af_b2, af_a2)
        th_af[ci] = theta_af_now

        # αβ-frame speed estimation (NO PARK TRANSFORM)
        if speed_method == 'pll':
            omega_elec_est = pll.step(theta_af_now)
            w_est[ci] = pll.omega_rpm
        elif speed_method == 'lpf1':
            dth = angle_diff_scalar(theta_af_now, theta_prev)
            raw_omega = dth / CL_TS
            omega_lpf1 = omega_lpf1*(1-lpf1_alpha) + raw_omega*lpf1_alpha
            omega_elec_est = omega_lpf1
            w_est[ci] = omega_elec_est / (2*np.pi*n_pp) * 60
        elif speed_method == 'lpf2':
            dth = angle_diff_scalar(theta_af_now, theta_prev)
            raw_omega = dth / CL_TS
            omega_elec_est = lpf.step(raw_omega)
            w_est[ci] = omega_elec_est / (2*np.pi*n_pp) * 60
        
        theta_prev = theta_af_now
        ci += 1

    n = ci; t = t_tr[:n]; ss = t > t[-1]*0.7
    ae = angle_error(th_true[:n], th_af[:n])
    ang = np.degrees(np.sqrt(np.mean(ae[ss]**2)))
    trk = np.sqrt(np.mean((w_true[:n][ss]-cmd_tr[:n][ss])**2))
    return {'t':t,'theta_true':th_true[:n],'omega_true':w_true[:n],
            'theta_af':th_af[:n],'omega_est':w_est[:n],'cmd':cmd_tr[:n],
            'ang_rms':ang,'trk_rms':trk}


# ============ MAIN: Systematic sweep ============
if __name__ == '__main__':
    p = get_motor_preset('servo')
    cmd = p['cmd_rpm']
    
    print('='*75)
    print('  servo: αβ-frame Speed Estimation (break NSO-AF coupling)')
    print('='*75)

    configs = []
    # PLL speed tracker with various bw and zeta
    for bw in [5, 10, 20, 30, 50, 80, 100, 150]:
        for z in [1, 2, 3, 5, 8, 10]:
            for cl in [50, 100, 200]:
                configs.append(('pll', bw, 0, cl, z, f'PLL bw={bw} z={z} CL{cl}'))
    
    # LPF1 with various tau and zeta
    for tau_ms in [1, 2, 5, 10, 20]:
        for z in [1, 2, 3, 5, 10]:
            for cl in [50, 100, 200]:
                configs.append(('lpf1', 0, tau_ms/1000, cl, z, f'LPF1 t={tau_ms}ms z={z} CL{cl}'))
    
    # LPF2 (cascaded) 
    for tau_ms in [2, 5, 10, 20]:
        for z in [1, 2, 3, 5]:
            for cl in [50, 100]:
                configs.append(('lpf2', 0, tau_ms/1000, cl, z, f'LPF2 t={tau_ms}ms z={z} CL{cl}'))

    header = f'{"Label":<28s} | {"AngSS":>7s} {"TrkSS":>7s}'
    print(header); print('-' * 48)

    best_ang = 999; best_trk = 999; best_label = ''; best_res = None
    results_list = []
    
    for method, bw, tau, cl, z, label in configs:
        try:
            r = run_ab_speed_s3(
                p['d'], CLBW_Hz=cl, zeta=z,
                af_Kp=p['af_kp'], af_Ki=p['af_ki'],
                R_mis=1.0, v_off=(0,0),
                cmd_rpm=cmd, load_step=p['load_step'],
                speed_method=method, speed_bw=bw, lpf_tau=tau,
                verbose=False)
            flag = ' ★' if r['ang_rms']<5 else (' ◆' if r['ang_rms']<15 else (' ·' if r['ang_rms']<30 else ''))
            if r['ang_rms'] < 30 or r['trk_rms'] < 100:
                print(f'{label:<28s} | {r["ang_rms"]:>7.2f} {r["trk_rms"]:>7.2f}{flag}')
            results_list.append((label, r['ang_rms'], r['trk_rms']))
            # Best by combined metric: angle < 20 AND tracking < 50
            if r['ang_rms'] < best_ang:
                best_ang = r['ang_rms']; best_label = label; best_res = r; best_trk = r['trk_rms']
        except Exception as e:
            pass
    
    print(f'\nBest overall: {best_label} -> Ang={best_ang:.2f} Trk={best_trk:.2f}')
    
    # Find best with tracking < 50
    good = [(l,a,t) for l,a,t in results_list if t < 50]
    if good:
        good.sort(key=lambda x: x[1])
        print(f'\nBest with Trk<50: {good[0][0]} -> Ang={good[0][1]:.2f} Trk={good[0][2]:.2f}')
    
    good2 = [(l,a,t) for l,a,t in results_list if t < 100]
    if good2:
        good2.sort(key=lambda x: x[1])
        print(f'Best with Trk<100: {good2[0][0]} -> Ang={good2[0][1]:.2f} Trk={good2[0][2]:.2f}')

    # Plot best result
    if best_res and best_ang < 50:
        fig, axes = plt.subplots(3,1,figsize=(14,9),sharex=True)
        t = best_res['t']
        axes[0].plot(t, best_res['cmd'],'r--',alpha=0.5,label='cmd')
        axes[0].plot(t, best_res['omega_true'],'b',label='true')
        axes[0].plot(t, best_res['omega_est'],'g',alpha=0.7,label='est')
        axes[0].set_ylabel('Speed [rpm]'); axes[0].legend()
        axes[0].set_title(f'servo S3 αβ-Speed | {best_label} | Ang={best_ang:.2f}')
        ae = np.degrees(angle_error(best_res['theta_true'], best_res['theta_af']))
        axes[1].plot(t, ae, 'r'); axes[1].set_ylabel('Angle err [deg]')
        axes[2].plot(t, best_res['omega_true']-best_res['cmd'],'b')
        axes[2].set_ylabel('Track err [rpm]'); axes[2].set_xlabel('Time [s]')
        plt.tight_layout(); plt.savefig('fig_ab_speed_servo_best.png',dpi=150); plt.close()
        print('Saved fig_ab_speed_servo_best.png')
