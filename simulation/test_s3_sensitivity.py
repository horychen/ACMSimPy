#!/usr/bin/env python
"""Test Stage 3 sensitivity to R mismatch and voltage offset."""
from demo_sensorless_active_flux import run_sensorless_demo, angle_error
from eval_staged_tuning import get_motor_preset, compute_metrics
import numpy as np

p = get_motor_preset('servo')
configs = [
    ('ideal',        1.0, 0.0, 0.0),
    ('R_only_1.5',   1.5, 0.0, 0.0),
    ('R_only_1.2',   1.2, 0.0, 0.0),
    ('R_only_1.1',   1.1, 0.0, 0.0),
    ('Voff_only',    1.0, 0.5, -0.3),
    ('Voff_small',   1.0, 0.1, -0.05),
    ('R1.1+Voff_sm', 1.1, 0.1, -0.05),
    ('R1.5+Voff',    1.5, 0.5, -0.3),
]

header = f'{"Config":<16s} | {"AngRMS_SS":>10s} {"SpdRMS_SS":>10s} {"TrkRMS_SS":>10s} {"RevSettle":>10s}'
print(header)
print('-' * 70)

for name, r_mis, v_a, v_b in configs:
    try:
        results = run_sensorless_demo(
            p['d'], zeta=p['zeta'], CLBW_Hz=p['clbw'],
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mismatch_factor=r_mis,
            voltage_offset_alpha=v_a, voltage_offset_beta=v_b,
            eso_omega_ob=80.0, speed_observer='nso',
            use_sensorless_speed=True,
            use_sensorless_angle=True,
            use_sensorless_torque_ff=True,
            cmd_rpm_ref=p['cmd_rpm'], load_step=p['load_step'],
            verbose=False)
        m = compute_metrics(results, p['cmd_rpm'])
        rev = f'{m["reversal_settling_s"]:.3f}' if m['reversal_settling_s'] < 999 else 'FAIL'
        print(f'{name:<16s} | {m["angle_rms_ss_deg"]:>10.2f} {m["spd_est_rms_ss"]:>10.2f} {m["trk_rms_ss"]:>10.2f} {rev:>10s}')
    except Exception as e:
        print(f'{name:<16s} | ERROR: {str(e)[:50]}')
