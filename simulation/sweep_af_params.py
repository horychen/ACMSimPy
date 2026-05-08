#!/usr/bin/env python
"""
Sweep AF observer Kp/Ki parameters for the servo motor.
Goal: find the (Kp, Ki) pair that minimizes steady-state angle estimation RMS.
"""
import sys, copy, time
import numpy as np

from demo_sensorless_active_flux import run_sensorless_demo, angle_error
from eval_staged_tuning import get_motor_preset

motor = 'servo'
preset = get_motor_preset(motor)

# Parameter grid
af_kps = [100, 200, 500, 1000, 2000, 5000]
af_kis = [500, 1000, 2000, 5000, 10000, 20000, 50000]

print(f'Sweeping AF Kp/Ki for {motor} motor...')
print(f'{"Kp":>8s} {"Ki":>8s} | {"AngleRMS_SS":>12s} {"AngleRMS_all":>12s} {"AngleMax_SS":>12s} | {"SpdRMS_SS":>10s} {"Stable":>6s}')
print('-' * 90)

best_rms = 999
best_kp, best_ki = None, None

t_total = time.time()

for kp in af_kps:
    for ki in af_kis:
        try:
            p = preset
            results = run_sensorless_demo(
                p['d'],
                zeta=p['zeta'],
                CLBW_Hz=p['clbw'],
                af_Kp=kp,
                af_Ki=ki,
                R_mismatch_factor=p['r_mis'],
                voltage_offset_alpha=p['v_off'][0],
                voltage_offset_beta=p['v_off'][1],
                eso_omega_ob=200.0,
                speed_observer='nso',
                use_sensorless_speed=False,
                use_sensorless_torque_ff=False,
                pure_p_current=False,
                cmd_rpm_ref=p['cmd_rpm'],
                load_step=p['load_step'],
                verbose=False,
            )

            t = results['t']
            theta_true = results['theta_true']
            theta_af = results['theta_af']
            omega_true = results['omega_true']
            omega_obs = results['omega_af_eso']
            err_af_deg = np.abs(angle_error(theta_af, theta_true)) * 180 / np.pi

            mask_ss = ((t >= 0.4) & (t <= 0.55)) | ((t >= 1.5) & (t <= 1.9))
            mask_all = (t >= 0.1) & (t <= 1.9)

            rms_ss = np.sqrt(np.mean(err_af_deg[mask_ss]**2))
            rms_all = np.sqrt(np.mean(err_af_deg[mask_all]**2))
            max_ss = np.max(err_af_deg[mask_ss])
            
            spd_err = omega_obs - omega_true
            spd_rms_ss = np.sqrt(np.mean(spd_err[mask_ss]**2))

            stable = not np.any(np.isnan(omega_true))
            stable_str = 'YES' if stable else 'NO'

            marker = ' ***' if rms_ss < best_rms else ''
            print(f'{kp:>8.0f} {ki:>8.0f} | {rms_ss:>12.2f} {rms_all:>12.2f} {max_ss:>12.2f} | {spd_rms_ss:>10.2f} {stable_str:>6s}{marker}')

            if rms_ss < best_rms and stable:
                best_rms = rms_ss
                best_kp, best_ki = kp, ki
        except Exception as e:
            print(f'{kp:>8.0f} {ki:>8.0f} | {"ERROR":>12s} {str(e)[:40]}')

elapsed = time.time() - t_total
print(f'\nTotal sweep time: {elapsed:.0f} s')
print(f'\nBest: Kp={best_kp}, Ki={best_ki}, Angle RMS (SS) = {best_rms:.2f} deg')
