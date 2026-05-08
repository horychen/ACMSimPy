#!/usr/bin/env python
"""
Sweep AF observer PI gains for Stage 3 full sensorless on servo motor.
Looking for a Kp/Ki combo that keeps AF angle stable in closed loop.
"""
from demo_sensorless_active_flux import run_sensorless_demo, angle_error
from eval_staged_tuning import get_motor_preset, compute_metrics
import numpy as np

p = get_motor_preset('servo')

# Sweep configurations: (af_Kp, af_Ki)
kp_vals = [100, 500, 1000, 2000, 5000]
ki_vals = [500, 1000, 5000, 10000, 50000]

header = f'{"Kp":>6s} {"Ki":>8s} | {"AngRMS_SS":>10s} {"AngMax_SS":>10s} {"SpdRMS_SS":>10s} {"TrkRMS_SS":>10s}'
print(header)
print('-' * 62)

best_ang = 999
best_cfg = None

for kp in kp_vals:
    for ki in ki_vals:
        try:
            results = run_sensorless_demo(
                p['d'], zeta=p['zeta'], CLBW_Hz=p['clbw'],
                af_Kp=kp, af_Ki=ki,
                R_mismatch_factor=1.0,  # ideal params first
                voltage_offset_alpha=0.0,
                voltage_offset_beta=0.0,
                eso_omega_ob=80.0, speed_observer='nso',
                use_sensorless_speed=True,
                use_sensorless_angle=True,
                use_sensorless_torque_ff=False,  # no TL FF to isolate
                cmd_rpm_ref=p['cmd_rpm'], load_step=p['load_step'],
                verbose=False)
            m = compute_metrics(results, p['cmd_rpm'])
            ang = m['angle_rms_ss_deg']
            print(f'{kp:>6d} {ki:>8d} | {ang:>10.2f} {m["angle_max_ss_deg"]:>10.2f} {m["spd_est_rms_ss"]:>10.2f} {m["trk_rms_ss"]:>10.2f}')
            if ang < best_ang:
                best_ang = ang
                best_cfg = (kp, ki)
        except Exception as e:
            print(f'{kp:>6d} {ki:>8d} | ERROR: {str(e)[:40]}')

print(f'\nBest: Kp={best_cfg[0]}, Ki={best_cfg[1]}, angle_rms_ss={best_ang:.2f} deg')
