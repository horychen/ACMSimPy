#!/usr/bin/env python
"""
Progressive closure test for servo:
1. Stage 0:   encoder speed + encoder angle          (baseline)
2. Stage 2b:  NSO speed    + encoder angle           (speed only)
3. Stage 2c:  encoder speed + AF angle               (angle only)
4. Stage 3:   NSO speed    + AF angle                (full sensorless)

This identifies which feedback path causes the instability.
Also test with different CLBW to find the stability boundary.
"""
import sys; sys.path.insert(0, '.')
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics

p = get_motor_preset('servo')
cmd = p['cmd_rpm']

configs = [
    # (use_speed, use_angle, use_ff, CLBW, zeta, label)
    # Stage 0: full encoder
    (False, False, False, 100, 10, 'S0 CL100'),
    (False, False, False, 200, 10, 'S0 CL200'),
    (False, False, False, 500, 10, 'S0 CL500'),
    
    # Stage 2b: NSO speed, encoder angle
    (True,  False, False, 50,  10, 'S2b CL50'),
    (True,  False, False, 100, 10, 'S2b CL100'),
    (True,  False, False, 200, 10, 'S2b CL200'),
    (True,  False, False, 500, 10, 'S2b CL500'),
    (True,  False, True,  100, 10, 'S2b+FF CL100'),
    (True,  False, True,  200, 10, 'S2b+FF CL200'),
    
    # Stage 2c: encoder speed, AF angle (NEW - isolate angle path)
    (False, True,  False, 50,  10, 'S2c CL50'),
    (False, True,  False, 100, 10, 'S2c CL100'),
    (False, True,  False, 200, 10, 'S2c CL200'),
    (False, True,  False, 500, 10, 'S2c CL500'),
    (False, True,  False, 30,  10, 'S2c CL30'),
    (False, True,  False, 20,  10, 'S2c CL20'),
    (False, True,  False, 10,  10, 'S2c CL10'),
    
    # Stage 3: full sensorless
    (True,  True,  False, 50,  10, 'S3 CL50'),
    (True,  True,  False, 100, 10, 'S3 CL100'),
    (True,  True,  False, 200, 10, 'S3 CL200'),
    (True,  True,  True,  100, 10, 'S3+FF CL100'),
    (True,  True,  True,  50,  10, 'S3+FF CL50'),
    (True,  True,  False, 30,  10, 'S3 CL30'),
    (True,  True,  False, 20,  10, 'S3 CL20'),
    (True,  True,  False, 10,  10, 'S3 CL10'),
]

print(f'servo progressive closure (ideal params)')
header = f'{"Label":<16s} | {"AngSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print(header); print('-' * 52)

for use_spd, use_ang, use_ff, clbw, z, label in configs:
    try:
        r = run_sensorless_demo(
            p['d'], zeta=z, CLBW_Hz=clbw,
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mismatch_factor=1.0,
            voltage_offset_alpha=0, voltage_offset_beta=0,
            eso_omega_ob=80, speed_observer='nso',
            use_sensorless_speed=use_spd,
            use_sensorless_angle=use_ang,
            use_sensorless_torque_ff=use_ff,
            cmd_rpm_ref=cmd, load_step=p['load_step'], verbose=False)
        m = compute_metrics(r, cmd)
        rev = f'{m["reversal_settling_s"]:.3f}' if m['reversal_settling_s'] < 999 else 'FAIL'
        flag = ' ★' if m['angle_rms_ss_deg'] < 5 else (' ◆' if m['angle_rms_ss_deg'] < 15 else '')
        print(f'{label:<16s} | {m["angle_rms_ss_deg"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}{flag}')
    except Exception as e:
        print(f'{label:<16s} | ERROR: {str(e)[:40]}')
