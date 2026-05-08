#!/usr/bin/env python
"""
Sweep NSO observer bandwidth (omega_ob) for the servo motor.
Run in Stage 2b mode (observer speed + encoder angle) to isolate speed estimation quality.
"""
import sys, time
import numpy as np
from demo_sensorless_active_flux import run_sensorless_demo, angle_error
from eval_staged_tuning import get_motor_preset

motor = 'servo'
preset = get_motor_preset(motor)
p = preset

omega_obs = [20, 50, 80, 100, 150, 200, 300, 500]

print(f'Sweeping NSO ω_ob for {motor} motor (Stage 2b: observer speed + encoder angle)...')
print(f'{"ω_ob":>8s} | {"SpdRMS_SS":>10s} {"SpdMax_SS":>10s} {"TrkRMS_SS":>10s} {"TrkMax_SS":>10s} {"AngRMS_SS":>10s} {"LoadDip":>8s} {"RevSet":>7s} | {"Stable":>6s}')
print('-' * 100)

for wob in omega_obs:
    try:
        results = run_sensorless_demo(
            p['d'],
            zeta=p['zeta'],
            CLBW_Hz=p['clbw'],
            af_Kp=p['af_kp'],
            af_Ki=p['af_ki'],
            R_mismatch_factor=p['r_mis'],
            voltage_offset_alpha=p['v_off'][0],
            voltage_offset_beta=p['v_off'][1],
            eso_omega_ob=wob,
            speed_observer='nso',
            use_sensorless_speed=True,
            use_sensorless_angle=False,  # encoder angle (Stage 2b)
            use_sensorless_torque_ff=False,
            pure_p_current=False,
            cmd_rpm_ref=p['cmd_rpm'],
            load_step=p['load_step'],
            verbose=False,
        )

        t = results['t']
        omega_true = results['omega_true']
        omega_obs_arr = results['omega_af_eso']
        theta_true = results['theta_true']
        theta_af = results['theta_af']
        cmd_rpm = results['cmd_rpm']
        npp = results['n_pp']

        mask_ss = ((t >= 0.4) & (t <= 0.55)) | ((t >= 1.5) & (t <= 1.9))
        
        err_af_deg = np.abs(angle_error(theta_af, theta_true)) * 180 / np.pi
        spd_err = omega_obs_arr - omega_true
        trk_err = omega_true - cmd_rpm

        spd_rms_ss = np.sqrt(np.mean(spd_err[mask_ss]**2))
        spd_max_ss = np.max(np.abs(spd_err[mask_ss]))
        trk_rms_ss = np.sqrt(np.mean(trk_err[mask_ss]**2))
        trk_max_ss = np.max(np.abs(trk_err[mask_ss]))
        ang_rms_ss = np.sqrt(np.mean(err_af_deg[mask_ss]**2))

        # Load disturbance rejection
        mask_pre = (t >= 0.55) & (t <= 0.60)
        mask_post = (t >= 0.60) & (t <= 0.80)
        if mask_pre.any() and mask_post.any():
            spd_before = np.mean(omega_true[mask_pre])
            load_dip = abs(spd_before - np.min(omega_true[mask_post]))
        else:
            load_dip = 999

        # Reversal settling
        mask_rev = (t >= 1.0) & (t <= 1.9)
        rev_idx = np.where(mask_rev)[0]
        target_rpm = -p['cmd_rpm']
        settled = np.abs(omega_true[rev_idx] - target_rpm) < 0.05 * p['cmd_rpm']
        rev_settle = (t[rev_idx[np.argmax(settled)]] - 1.0) if settled.any() else 999

        stable = not np.any(np.isnan(omega_true))
        stable_str = 'YES' if stable else 'NO'

        print(f'{wob:>8.0f} | {spd_rms_ss:>10.2f} {spd_max_ss:>10.2f} {trk_rms_ss:>10.2f} {trk_max_ss:>10.2f} {ang_rms_ss:>10.2f} {load_dip:>8.2f} {rev_settle:>7.3f} | {stable_str:>6s}')
    except Exception as e:
        print(f'{wob:>8.0f} | ERROR: {str(e)[:60]}')

print('\nDone.')
