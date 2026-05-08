#!/usr/bin/env python
"""
Test: Stage 3 with LPF speed (from AF angle derivative) instead of NSO.
This eliminates the NSO-AF coupling completely.
"""
import matplotlib; matplotlib.use('Agg')
from pylab import np, plt
import copy, time as _time

from tutorials_ep6_svpwm import (The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental)
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator
from demo_sensorless_active_flux import angle_error, angle_diff_scalar
from eval_staged_tuning import get_motor_preset


def run_lpf_speed_s3(d, CLBW_Hz, zeta, af_Kp, af_Ki, R_mis, v_off,
                     cmd_rpm, load_step, lpf_tau=0.01, verbose=True):
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

    # Full sensorless
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

    af_psi_s = np.array([KE, 0.0]); af_corr_int = np.zeros(2)
    lpf_alpha = CL_TS / (lpf_tau + CL_TS)
    omega_lpf = 0.0; theta_prev = 0.0

    total_time = dd['TIME_SLICE'] * dd['NUMBER_OF_SLICES']
    n_ctrl = int(total_time / CL_TS) + 100
    t_tr = np.zeros(n_ctrl); th_true = np.zeros(n_ctrl); w_true = np.zeros(n_ctrl)
    th_af = np.zeros(n_ctrl); w_lpf = np.zeros(n_ctrl); cmd_tr = np.zeros(n_ctrl)

    CTRL.cmd_rpm = 0.0; ci = 0

    for sl in range(int(total_time/CL_TS)):
        t0 = sl * CL_TS; tm = t0 + CL_TS*0.5
        if tm<0.3: CTRL.cmd_rpm = cmd_rpm*min(tm/0.3,1); ACM.TLoad = 0
        elif tm<0.6: CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0
        elif tm<0.8: CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = load_step
        elif tm<1.0: CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0
        elif tm<1.3: CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0
        else: CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0

        if ci > 0:
            CTRL.vartheta_d = th_af[ci-1]
            CTRL.xS[1] = omega_lpf  # LPF speed in elec rad/s

        mt, wd = ACMSimPyIncremental(t0=t0, TIME=CL_TS, ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_spd)

        if ci >= n_ctrl: break
        t_tr[ci] = mt[0]; th_true[ci] = wd[0][0]; w_true[ci] = wd[1][0]; cmd_tr[ci] = wd[12][0]
        ia = wd[6][0]; ib = wd[7][0]; ua = wd[28][0]+v_off[0]; ub = wd[29][0]+v_off[1]

        # AF observer
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

        # LPF speed from AF angle derivative
        dth = angle_diff_scalar(theta_af_now, theta_prev)
        raw_omega = dth / CL_TS
        omega_lpf = omega_lpf * (1-lpf_alpha) + raw_omega * lpf_alpha
        theta_prev = theta_af_now
        w_lpf[ci] = omega_lpf / (2*np.pi*n_pp) * 60

        ci += 1

    n = ci; t = t_tr[:n]; ss = t > t[-1]*0.75
    ae = angle_error(th_true[:n], th_af[:n])
    ang = np.degrees(np.sqrt(np.mean(ae[ss]**2)))
    trk = np.sqrt(np.mean((w_true[:n][ss]-cmd_tr[:n][ss])**2))
    return {'t':t,'theta_true':th_true[:n],'omega_true':w_true[:n],
            'theta_af':th_af[:n],'omega_lpf':w_lpf[:n],'cmd':cmd_tr[:n],
            'ang_rms':ang,'trk_rms':trk,'n_pp':n_pp}


# ============ MAIN ============
print('='*70)
print('  LPF Speed (from AF angle) Stage 3 on servo & big_L')
print('='*70)

for motor_name in ['servo', 'big_L']:
    p = get_motor_preset(motor_name)
    cmd = p['cmd_rpm']
    print(f'\n--- {motor_name} ---')
    header = f'{"Label":<24s} | {"AngSS":>7s} {"TrkSS":>7s}'
    print(header); print('-' * 44)

    taus = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05]
    clbws = [30, 50, 100, 200] if motor_name == 'servo' else [10, 20, 50]
    
    best_ang = 999; best_label = ''; best_res = None
    for clbw in clbws:
        for tau in taus:
            label = f'CL{clbw} tau={tau*1e3:.0f}ms'
            try:
                r = run_lpf_speed_s3(
                    p['d'], CLBW_Hz=clbw, zeta=10,
                    af_Kp=p['af_kp'], af_Ki=p['af_ki'],
                    R_mis=1.0, v_off=(0,0),
                    cmd_rpm=cmd, load_step=p['load_step'],
                    lpf_tau=tau, verbose=False)
                flag = ' ★' if r['ang_rms'] < 5 else (' ◆' if r['ang_rms'] < 15 else '')
                print(f'{label:<24s} | {r["ang_rms"]:>7.2f} {r["trk_rms"]:>7.2f}{flag}')
                if r['ang_rms'] < best_ang:
                    best_ang = r['ang_rms']; best_label = label; best_res = r
            except Exception as e:
                print(f'{label:<24s} | ERROR: {str(e)[:35]}')
    
    print(f'Best: {best_label} -> {best_ang:.2f} deg')
    
    if best_res and best_ang < 30:
        fig, axes = plt.subplots(3,1,figsize=(14,9),sharex=True)
        t = best_res['t']
        axes[0].plot(t, best_res['cmd'],'r--',alpha=0.5,label='cmd')
        axes[0].plot(t, best_res['omega_true'],'b',label='true')
        axes[0].plot(t, best_res['omega_lpf'],'g',alpha=0.7,label='LPF')
        axes[0].set_ylabel('Speed [rpm]'); axes[0].legend()
        axes[0].set_title(f'{motor_name} S3 LPF-Speed | {best_label} | Ang={best_ang:.2f}')
        ae = np.degrees(angle_error(best_res['theta_true'], best_res['theta_af']))
        axes[1].plot(t, ae, 'r'); axes[1].set_ylabel('Angle err [deg]')
        axes[2].plot(t, best_res['omega_true']-best_res['cmd'],'b'); axes[2].set_ylabel('Track err [rpm]')
        axes[2].set_xlabel('Time [s]')
        plt.tight_layout(); plt.savefig(f'fig_lpf_s3_{motor_name}.png',dpi=150); plt.close()
        print(f'Saved fig_lpf_s3_{motor_name}.png')
