#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test PLL and SMO observers on servo and big_L motors.
Replaces the AF angle extraction with PLL/SMO while keeping NSO for speed.
"""
import matplotlib
matplotlib.use('Agg')
from pylab import np, plt
import copy, time as _time

from tutorials_ep6_svpwm import (
    The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental,
)
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator
from demo_sensorless_active_flux import wrap_angle, angle_error, angle_diff_scalar, NaturalSpeedObserver
from observers_alt import PLLFluxObserver, SlidingModeObserver


def run_alt_observer_sim(d, observer_type='pll', CLBW_Hz=100, zeta=10,
                         R_mismatch_factor=1.5, v_off=(0.5, -0.3),
                         omega_ob=80.0, cmd_rpm=500, load_step=1.27,
                         pll_bw=100.0, smo_gain=None, smo_lpf=None,
                         af_Kp=500.0, af_Ki=5000.0, verbose=True):
    """Run simulation with alternative observer (PLL or SMO) for angle estimation."""
    dd = copy.deepcopy(d)
    R_true = dd['init_R']
    R_obs = R_true * R_mismatch_factor
    L = dd['init_Lq']
    J_s = dd['init_Js']
    n_pp = dd['init_npp']
    KE = dd['init_KE']
    CL_TS = dd['CL_TS']
    VL_TS = CL_TS * dd['VL_EXE_PER_CL_EXE']

    currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R_true, L, CLBW_Hz)
    currentBW = currentKp / L
    speedKp, speedKi = get_coeffs_dc_motor_SPEED_regulator(J_s, n_pp, KE, zeta, currentBW)

    dd['CL_SERIES_KP'] = currentKp
    dd['CL_SERIES_KI'] = currentKi
    dd['VL_SERIES_KP'] = speedKp
    dd['VL_SERIES_KI'] = speedKi

    CTRL = The_Motor_Controller(CL_TS=CL_TS, VL_TS=VL_TS, init_npp=n_pp,
        init_IN=dd['init_IN'], init_R=R_true, init_Ld=dd['init_Ld'],
        init_Lq=L, init_KE=KE, init_Rreq=0.0, init_Js=J_s,
        DC_BUS_VOLTAGE=dd['DC_BUS_VOLTAGE'])
    CTRL.bool_apply_decoupling_voltages_to_current_regulation = False
    CTRL.bool_apply_sweeping_frequency_excitation = False
    CTRL.bool_overwrite_speed_commands = True
    CTRL.bool_zero_id_control = True
    CTRL.bool_apply_speed_closed_loop_control = True

    # Full sensorless: angle + speed + torque FF
    CTRL.index_separate_speed_estimation = 1
    CTRL.ell1 = CTRL.ell2 = CTRL.ell3 = CTRL.ell4 = 0.0
    CTRL.bool_use_sensorless_theta = 1
    CTRL.use_disturbance_feedforward_rejection = 1

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=dd['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])

    ki_factor = dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False', 10)
    vl_lim = dd.get('VL_LIMIT_OVERLOAD_FACTOR', 10.0)
    Vmax = dd['DC_BUS_VOLTAGE'] / 1.732
    Imax = vl_lim * 1.414 * dd['init_IN']
    reg_id = The_PID_Regulator(currentKp, currentKp*currentKi*ki_factor, 0,0, Vmax, Vmax, CL_TS)
    reg_iq = The_PID_Regulator(currentKp, currentKp*currentKi*ki_factor, 0,0, Vmax, Vmax, CL_TS)
    reg_speed = The_PID_Regulator(speedKp, speedKp*speedKi, 0,0, Imax, Imax, VL_TS)

    MACHINE_TS = CL_TS / ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    controller_down_sampling = int(CL_TS / MACHINE_TS)

    # Create angle observer
    if observer_type == 'pll':
        angle_obs = PLLFluxObserver(R_obs, L, KE, CL_TS, pll_bw=pll_bw, af_Kp=af_Kp, af_Ki=af_Ki)
    else:
        angle_obs = SlidingModeObserver(R_obs, L, KE, CL_TS, smo_gain=smo_gain, lpf_fc=smo_lpf, pll_bw=pll_bw)

    # Speed observer (NSO)
    nso = NaturalSpeedObserver(omega_ob=omega_ob, npp=n_pp, Js=J_s,
                               R=R_obs, Ld=dd['init_Ld'], Lq=L, KE=KE, dt=CL_TS)

    total_time = dd['TIME_SLICE'] * dd['NUMBER_OF_SLICES']
    n_ctrl = int(total_time / CL_TS) + 100
    t_trace = np.zeros(n_ctrl)
    theta_true = np.zeros(n_ctrl)
    omega_true = np.zeros(n_ctrl)
    theta_est = np.zeros(n_ctrl)
    omega_est = np.zeros(n_ctrl)
    cmd_trace = np.zeros(n_ctrl)
    load_true_arr = np.zeros(n_ctrl)

    CTRL.cmd_rpm = 0.0
    ctrl_idx = 0
    t_start = _time.time()

    sim_TIME_SLICE = CL_TS
    sim_N_SLICES = int(total_time / CL_TS)

    for sl in range(sim_N_SLICES):
        t0 = sl * sim_TIME_SLICE
        t_mid = t0 + sim_TIME_SLICE * 0.5

        if t_mid < 0.3:
            CTRL.cmd_rpm = cmd_rpm * min(t_mid/0.3, 1.0); ACM.TLoad = 0.0
        elif t_mid < 0.6:
            CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0.0
        elif t_mid < 0.8:
            CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = load_step
        elif t_mid < 1.0:
            CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0.0
        elif t_mid < 1.3:
            CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0.0
        else:
            CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0.0

        # Inject observer outputs
        if ctrl_idx > 0:
            CTRL.vartheta_d = theta_est[ctrl_idx - 1]
            CTRL.xS[1] = nso.omega_elec
            CTRL.xS[2] = nso.disturbance_est

        machine_times, watch_data = ACMSimPyIncremental(
            t0=t0, TIME=sim_TIME_SLICE, ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed)

        for k in range(0, len(machine_times), controller_down_sampling):
            if ctrl_idx >= n_ctrl: break
            t_trace[ctrl_idx] = machine_times[k]
            theta_true[ctrl_idx] = watch_data[0][k]
            omega_true[ctrl_idx] = watch_data[1][k]
            cmd_trace[ctrl_idx] = watch_data[12][k]
            load_true_arr[ctrl_idx] = ACM.TLoad

            ia = watch_data[6][k]
            ib = watch_data[7][k]
            ua = watch_data[28][k] + v_off[0]
            ub = watch_data[29][k] + v_off[1]

            # Run angle observer
            th_now, om_now, _ = angle_obs.step(ia, ib, ua, ub)

            # Run NSO (needs dq quantities in AF frame)
            cos_th = np.cos(th_now); sin_th = np.sin(th_now)
            id_m = ia*cos_th + ib*sin_th
            iq_m = -ia*sin_th + ib*cos_th
            uq_m = -ua*sin_th + ub*cos_th
            nso.step(iq_m, uq_m, id_m)

            theta_est[ctrl_idx] = th_now
            omega_est[ctrl_idx] = nso.omega_rpm
            ctrl_idx += 1

    n = ctrl_idx
    elapsed = _time.time() - t_start
    if verbose:
        print(f'  {observer_type.upper()} sim done in {elapsed:.1f}s, {n} steps')

    return {
        't': t_trace[:n], 'theta_true': theta_true[:n], 'omega_true': omega_true[:n],
        'theta_est': theta_est[:n], 'omega_est': omega_est[:n],
        'cmd_rpm': cmd_trace[:n], 'load_true': load_true_arr[:n], 'n_pp': n_pp,
    }


def eval_results(res, cmd_rpm):
    """Compute metrics from results."""
    t = res['t']; n = len(t)
    ss_mask = t > t[-1] * 0.75
    
    ang_err = angle_error(res['theta_true'], res['theta_est'])
    ang_rms_ss = np.degrees(np.sqrt(np.mean(ang_err[ss_mask]**2)))
    
    spd_err = res['omega_est'] - res['omega_true']
    spd_rms_ss = np.sqrt(np.mean(spd_err[ss_mask]**2))
    
    trk_err = res['omega_true'] - res['cmd_rpm']
    trk_rms_ss = np.sqrt(np.mean(trk_err[ss_mask]**2))
    
    # Load dip
    t_load_on = None
    for i in range(1, n):
        if res['load_true'][i] > 0 and res['load_true'][i-1] == 0:
            t_load_on = t[i]; break
    if t_load_on is not None:
        mask_load = (t >= t_load_on) & (t < t_load_on + 0.5)
        if mask_load.any():
            dip = abs(res['omega_true'][mask_load].min() - cmd_rpm) if cmd_rpm > 0 else 0
        else: dip = 999
    else: dip = 999
    
    # Reversal
    t_rev = None
    for i in range(1, n):
        if res['cmd_rpm'][i] < 0 and res['cmd_rpm'][i-1] >= 0:
            t_rev = t[i]; break
    rev_s = 999
    if t_rev is not None:
        mask_rev = (t >= t_rev)
        w_rev = res['omega_true'][mask_rev]
        t_rev_arr = t[mask_rev]
        settled = np.where(np.abs(w_rev - (-cmd_rpm)) < 0.05 * cmd_rpm)[0]
        if len(settled) > 0:
            rev_s = t_rev_arr[settled[0]] - t_rev
    
    return ang_rms_ss, spd_rms_ss, trk_rms_ss, dip, rev_s


# =============================================================================
#  MAIN: Test on servo and big_L
# =============================================================================
if __name__ == '__main__':
    from eval_staged_tuning import get_motor_preset
    
    print('='*75)
    print('  PLL & SMO Observer Test on servo and big_L (Stage 3)')
    print('='*75)

    motors = ['servo', 'big_L']
    obs_types = ['pll', 'smo']
    
    for motor_name in motors:
        p = get_motor_preset(motor_name)
        cmd = p['cmd_rpm']
        KE = p['d']['init_KE']
        L = p['d']['init_Lq']
        R = p['d']['init_R']
        print(f'\n{"="*75}')
        print(f'  Motor: {motor_name} (KE/L={KE/p["d"]["init_Ld"]:.1f}, L/R={p["d"]["init_Ld"]/R*1e3:.1f}ms)')
        print(f'{"="*75}')

        # Sweep configurations
        if motor_name == 'servo':
            configs = [
                # (obs, CLBW, pll_bw, smo_gain, smo_lpf, omega_ob, zeta, label)
                ('pll', 100, 50,  None, None, 80, 10, 'PLL bw=50'),
                ('pll', 100, 100, None, None, 80, 10, 'PLL bw=100'),
                ('pll', 100, 200, None, None, 80, 10, 'PLL bw=200'),
                ('pll', 200, 100, None, None, 80, 10, 'PLL bw=100 CL200'),
                ('pll', 200, 200, None, None, 80, 10, 'PLL bw=200 CL200'),
                ('pll', 300, 100, None, None, 80, 10, 'PLL bw=100 CL300'),
                ('pll', 500, 200, None, None, 80, 10, 'PLL bw=200 CL500'),
                ('smo', 100, 50,  None, None, 80, 10, 'SMO pll=50'),
                ('smo', 100, 100, None, None, 80, 10, 'SMO pll=100'),
                ('smo', 100, 200, None, None, 80, 10, 'SMO pll=200'),
                ('smo', 200, 100, None, None, 80, 10, 'SMO pll=100 CL200'),
                ('smo', 200, 200, None, None, 80, 10, 'SMO pll=200 CL200'),
                ('smo', 300, 200, None, None, 80, 10, 'SMO pll=200 CL300'),
                ('smo', 500, 200, None, None, 80, 10, 'SMO pll=200 CL500'),
            ]
        else:  # big_L
            configs = [
                ('pll', 50, 30,   None, None, 80, 6.5, 'PLL bw=30'),
                ('pll', 50, 50,   None, None, 80, 6.5, 'PLL bw=50'),
                ('pll', 50, 100,  None, None, 80, 6.5, 'PLL bw=100'),
                ('pll', 80, 50,   None, None, 80, 6.5, 'PLL bw=50 CL80'),
                ('pll', 80, 100,  None, None, 80, 6.5, 'PLL bw=100 CL80'),
                ('smo', 50, 30,   None, 50,  80, 6.5, 'SMO lpf=50'),
                ('smo', 50, 50,   None, 100, 80, 6.5, 'SMO lpf=100'),
                ('smo', 80, 50,   None, 100, 80, 6.5, 'SMO lpf=100 CL80'),
                ('smo', 80, 100,  None, 200, 80, 6.5, 'SMO lpf=200 CL80'),
            ]

        header = f'{"Label":<22s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
        print(header)
        print('-' * 68)

        best_ang = 999; best_cfg = ''
        for obs, clbw, pll_bw, smo_g, smo_l, w_ob, z, label in configs:
            try:
                res = run_alt_observer_sim(
                    p['d'], observer_type=obs, CLBW_Hz=clbw, zeta=z,
                    R_mismatch_factor=p['r_mis'], v_off=p['v_off'],
                    omega_ob=w_ob, cmd_rpm=cmd, load_step=p['load_step'],
                    pll_bw=pll_bw, smo_gain=smo_g, smo_lpf=smo_l,
                    af_Kp=p['af_kp'], af_Ki=p['af_ki'], verbose=False)
                ang, spd, trk, dip, rev = eval_results(res, cmd)
                rev_s = f'{rev:.3f}' if rev < 999 else 'FAIL'
                print(f'{label:<22s} | {ang:>7.2f} {spd:>7.2f} {trk:>7.2f} {dip:>7.2f} {rev_s:>6s}')
                if ang < best_ang:
                    best_ang = ang; best_cfg = label; best_res = res
            except Exception as e:
                print(f'{label:<22s} | ERROR: {str(e)[:45]}')

        print(f'\n  Best: {best_cfg} -> Angle RMS = {best_ang:.2f} deg')
        
        # Plot best result
        if best_ang < 90:
            fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
            t = best_res['t']
            
            axes[0].plot(t, best_res['cmd_rpm'], '--', label='cmd', alpha=0.5)
            axes[0].plot(t, best_res['omega_true'], label='true')
            axes[0].plot(t, best_res['omega_est'], label='est', alpha=0.7)
            axes[0].set_ylabel('Speed [rpm]'); axes[0].legend(); axes[0].set_title(f'{motor_name} - {best_cfg}')
            
            ang_err = np.degrees(angle_error(best_res['theta_true'], best_res['theta_est']))
            axes[1].plot(t, ang_err)
            axes[1].set_ylabel('Angle error [deg]'); axes[1].set_title(f'Angle error (RMS SS={best_ang:.2f} deg)')
            
            trk_err = best_res['omega_true'] - best_res['cmd_rpm']
            axes[2].plot(t, trk_err)
            axes[2].set_ylabel('Tracking error [rpm]'); axes[2].set_xlabel('Time [s]')
            
            plt.tight_layout()
            plt.savefig(f'fig_alt_obs_{motor_name}_best.png', dpi=150)
            plt.close()
            print(f'  Saved fig_alt_obs_{motor_name}_best.png')
