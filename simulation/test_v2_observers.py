#!/usr/bin/env python
"""Test EMF-PLL and Adaptive SMO on servo motor."""
import matplotlib; matplotlib.use('Agg')
from pylab import np, plt
import copy, time as _time

from tutorials_ep6_svpwm import (The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental)
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator
from demo_sensorless_active_flux import angle_error, NaturalSpeedObserver
from observers_v2 import EMFDirectPLL, AdaptiveSMO
from eval_staged_tuning import get_motor_preset


def run_v2_sim(d, obs_obj, CLBW_Hz, zeta, omega_ob, cmd_rpm, load_step,
               R_mismatch_factor=1.0, v_off=(0,0)):
    dd = copy.deepcopy(d)
    R_true = dd['init_R']; L = dd['init_Lq']; J_s = dd['init_Js']
    n_pp = dd['init_npp']; KE = dd['init_KE']; CL_TS = dd['CL_TS']
    VL_TS = CL_TS * dd['VL_EXE_PER_CL_EXE']
    R_obs = R_true * R_mismatch_factor

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
    CTRL.index_separate_speed_estimation = 1
    CTRL.ell1 = CTRL.ell2 = CTRL.ell3 = CTRL.ell4 = 0.0
    CTRL.bool_use_sensorless_theta = 1
    CTRL.use_disturbance_feedforward_rejection = 1

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1)
    ki_f = dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False', 10)
    vl_lim = dd.get('VL_LIMIT_OVERLOAD_FACTOR', 10.0)
    Vm = dd['DC_BUS_VOLTAGE']/1.732; Im = vl_lim*1.414*dd['init_IN']
    reg_id = The_PID_Regulator(cKp, cKp*cKi*ki_f, 0,0, Vm, Vm, CL_TS)
    reg_iq = The_PID_Regulator(cKp, cKp*cKi*ki_f, 0,0, Vm, Vm, CL_TS)
    reg_spd = The_PID_Regulator(sKp, sKp*sKi, 0,0, Im, Im, VL_TS)

    nso = NaturalSpeedObserver(omega_ob=omega_ob, npp=n_pp, Js=J_s,
        R=R_obs, Ld=dd['init_Ld'], Lq=L, KE=KE, dt=CL_TS)

    total_time = dd['TIME_SLICE'] * dd['NUMBER_OF_SLICES']
    n_ctrl = int(total_time / CL_TS) + 100
    t_tr = np.zeros(n_ctrl); th_true = np.zeros(n_ctrl); w_true = np.zeros(n_ctrl)
    th_est = np.zeros(n_ctrl); w_est = np.zeros(n_ctrl); cmd_tr = np.zeros(n_ctrl)

    CTRL.cmd_rpm = 0.0; ci = 0
    for sl in range(int(total_time/CL_TS)):
        t0 = sl * CL_TS; tm = t0 + CL_TS*0.5
        if tm<0.3: CTRL.cmd_rpm=cmd_rpm*min(tm/0.3,1); ACM.TLoad=0
        elif tm<0.6: CTRL.cmd_rpm=cmd_rpm; ACM.TLoad=0
        elif tm<0.8: CTRL.cmd_rpm=cmd_rpm; ACM.TLoad=load_step
        elif tm<1.0: CTRL.cmd_rpm=cmd_rpm; ACM.TLoad=0
        elif tm<1.3: CTRL.cmd_rpm=-cmd_rpm; ACM.TLoad=0
        else: CTRL.cmd_rpm=-cmd_rpm; ACM.TLoad=0

        if ci > 0:
            CTRL.vartheta_d = th_est[ci-1]
            CTRL.xS[1] = nso.omega_elec
            CTRL.xS[2] = nso.disturbance_est

        mt, wd = ACMSimPyIncremental(t0=t0, TIME=CL_TS, ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_spd)

        if ci >= n_ctrl: break
        t_tr[ci] = mt[0]; th_true[ci] = wd[0][0]; w_true[ci] = wd[1][0]
        cmd_tr[ci] = wd[12][0]
        ia = wd[6][0]; ib = wd[7][0]
        ua = wd[28][0]+v_off[0]; ub = wd[29][0]+v_off[1]

        th_now, om_now, _ = obs_obj.step(ia, ib, ua, ub)
        cos_th = np.cos(th_now); sin_th = np.sin(th_now)
        nso.step(-ia*sin_th+ib*cos_th, -ua*sin_th+ub*cos_th, ia*cos_th+ib*sin_th)

        th_est[ci] = th_now; w_est[ci] = nso.omega_rpm
        ci += 1

    n = ci
    return {'t':t_tr[:n],'theta_true':th_true[:n],'omega_true':w_true[:n],
            'theta_est':th_est[:n],'omega_est':w_est[:n],'cmd_rpm':cmd_tr[:n],'n_pp':n_pp}


def quick_metrics(res, cmd):
    t = res['t']; ss = t > t[-1]*0.75
    ae = angle_error(res['theta_true'], res['theta_est'])
    ang = np.degrees(np.sqrt(np.mean(ae[ss]**2)))
    trk = np.sqrt(np.mean((res['omega_true'][ss]-res['cmd_rpm'][ss])**2))
    return ang, trk


# ============ MAIN ============
print('='*70)
print('  EMF-PLL & Adaptive SMO on servo')
print('='*70)

p = get_motor_preset('servo')
dd = p['d']; cmd = p['cmd_rpm']; ls = p['load_step']
R = dd['init_R']; L = dd['init_Lq']; KE = dd['init_KE']; npp = dd['init_npp']
CL_TS = dd['CL_TS']

configs = [
    # (type, CLBW, pll_bw, emf_gain/base_gain, lpf_fc, R_mis, v_off, zeta, w_ob, label)
    # EMF-PLL sweeps (ideal)
    ('emf', 100, 50,  None, 200, 1.0, (0,0), 10, 80, 'EMF pll50 CL100'),
    ('emf', 100, 100, None, 200, 1.0, (0,0), 10, 80, 'EMF pll100 CL100'),
    ('emf', 100, 200, None, 500, 1.0, (0,0), 10, 80, 'EMF pll200 CL100'),
    ('emf', 200, 100, None, 500, 1.0, (0,0), 10, 80, 'EMF pll100 CL200'),
    ('emf', 200, 200, None, 500, 1.0, (0,0), 10, 80, 'EMF pll200 CL200'),
    ('emf', 500, 200, None, 500, 1.0, (0,0), 10, 80, 'EMF pll200 CL500'),
    # Different EMF gains
    ('emf', 200, 100, 500, 500, 1.0, (0,0), 10, 80, 'EMF g500 pll100'),
    ('emf', 200, 100, 2000, 500, 1.0, (0,0), 10, 80, 'EMF g2000 pll100'),
    # With mismatch
    ('emf', 200, 100, None, 500, 1.5, (0.5,-0.3), 10, 80, 'EMF mis pll100'),
    ('emf', 200, 200, None, 500, 1.5, (0.5,-0.3), 10, 80, 'EMF mis pll200'),
    # Adaptive SMO
    ('asmo', 100, 50,  2.0, 200, 1.0, (0,0), 10, 80, 'ASMO pll50 CL100'),
    ('asmo', 100, 100, 2.0, 200, 1.0, (0,0), 10, 80, 'ASMO pll100 CL100'),
    ('asmo', 200, 100, 5.0, 300, 1.0, (0,0), 10, 80, 'ASMO g5 pll100'),
    ('asmo', 200, 200, 5.0, 300, 1.0, (0,0), 10, 80, 'ASMO g5 pll200'),
    ('asmo', 200, 100, 5.0, 300, 1.5, (0.5,-0.3), 10, 80, 'ASMO mis pll100'),
]

header = f'{"Label":<22s} | {"AngSS":>7s} {"TrkSS":>7s}'
print(header); print('-' * 42)

best_ang = 999; best_label = ''; best_res = None
for typ, clbw, pll_bw, gain, lpf, rmis, voff, z, wob, label in configs:
    try:
        if typ == 'emf':
            obs = EMFDirectPLL(R*rmis, L, KE, CL_TS, emf_gain=gain, pll_bw=pll_bw, emf_lpf_fc=lpf)
        else:
            obs = AdaptiveSMO(R*rmis, L, KE, npp, CL_TS, base_gain=gain, lpf_fc=lpf, pll_bw=pll_bw)
        res = run_v2_sim(dd, obs, clbw, z, wob, cmd, ls, R_mismatch_factor=rmis, v_off=voff)
        ang, trk = quick_metrics(res, cmd)
        print(f'{label:<22s} | {ang:>7.2f} {trk:>7.2f}')
        if ang < best_ang: best_ang = ang; best_label = label; best_res = res
    except Exception as e:
        print(f'{label:<22s} | ERROR: {str(e)[:40]}')

print(f'\nBest: {best_label} -> {best_ang:.2f} deg')

if best_res is not None and best_ang < 90:
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    t = best_res['t']
    axes[0].plot(t, best_res['cmd_rpm'], '--', alpha=0.5, label='cmd')
    axes[0].plot(t, best_res['omega_true'], label='true')
    axes[0].plot(t, best_res['omega_est'], alpha=0.7, label='est')
    axes[0].set_ylabel('Speed [rpm]'); axes[0].legend()
    axes[0].set_title(f'servo - {best_label} (Angle RMS={best_ang:.1f} deg)')

    ae = np.degrees(angle_error(best_res['theta_true'], best_res['theta_est']))
    axes[1].plot(t, ae); axes[1].set_ylabel('Angle error [deg]')

    axes[2].plot(t, best_res['omega_true'] - best_res['cmd_rpm'])
    axes[2].set_ylabel('Track err [rpm]'); axes[2].set_xlabel('Time [s]')
    plt.tight_layout(); plt.savefig('fig_v2_servo_best.png', dpi=150); plt.close()
    print('Saved fig_v2_servo_best.png')
