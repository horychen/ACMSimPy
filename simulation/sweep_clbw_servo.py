#!/usr/bin/env python
"""
Sweep CLBW for servo motor to find Stage 3 stable region.
Keep DC_BUS_VOLTAGE=200V (already set).
"""
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy

p_base = get_motor_preset('servo')

clbws = [20, 30, 50, 80, 100, 150, 200, 300, 500]

header = f'{"CLBW":>6s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print(f'servo motor: DC_BUS=200V, cmd_rpm=500')
print(header)
print('-' * 58)

for clbw in clbws:
    p = copy.deepcopy(p_base)
    try:
        results = run_sensorless_demo(
            p['d'], zeta=p['zeta'], CLBW_Hz=clbw,
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mismatch_factor=p['r_mis'],
            voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
            eso_omega_ob=p.get('omega_ob', 80.0), speed_observer='nso',
            use_sensorless_speed=True,
            use_sensorless_angle=True,
            use_sensorless_torque_ff=True,
            cmd_rpm_ref=500, load_step=p['load_step'],
            verbose=False)
        m = compute_metrics(results, 500)
        rev_s = m['reversal_settling_s']
        rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
        print(f'{clbw:>6d} | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
    except Exception as e:
        print(f'{clbw:>6d} | ERROR: {str(e)[:50]}')
