#!/usr/bin/env python
"""
Test small_L at 20V with different speed loop bandwidth (FOC_desired_VLBW_HZ) 
and current loop bandwidth (CLBW) to find stable gains.
"""
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy, numpy as np

p_base = get_motor_preset('small_L')

# Servo-like: lower speed loop BW to be less aggressive
configs = [
    # (label,   DC_BUS, CLBW, VLBW_HZ, FOC_delta)
    ('20V_baseline',    20,  1000, 120, 15),
    ('20V_VLBW60',      20,  1000,  60, 15),
    ('20V_VLBW30',      20,  1000,  30, 15),
    ('20V_VLBW15',      20,  1000,  15, 15),
    ('20V_CL500_VL60',  20,   500,  60, 15),
    ('20V_CL500_VL30',  20,   500,  30, 15),
    ('20V_CL200_VL30',  20,   200,  30, 15),
    ('10V_baseline',    10,  1000, 120, 15),
    ('10V_VLBW60',      10,  1000,  60, 15),
    ('10V_VLBW30',      10,  1000,  30, 15),
]

header = f'{"Config":<18s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print(header)
print('-' * 66)

for name, vbus, clbw, vlbw, delta in configs:
    p = copy.deepcopy(p_base)
    p['d']['DC_BUS_VOLTAGE'] = vbus
    p['d']['FOC_desired_VLBW_HZ'] = vlbw
    p['d']['FOC_delta'] = delta
    try:
        results = run_sensorless_demo(
            p['d'], zeta=p['zeta'], CLBW_Hz=clbw,
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mismatch_factor=p['r_mis'],
            voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
            eso_omega_ob=p.get('omega_ob', 200.0), speed_observer='nso',
            use_sensorless_speed=True,
            use_sensorless_angle=True,
            use_sensorless_torque_ff=True,
            cmd_rpm_ref=200, load_step=p['load_step'],
            verbose=False)
        m = compute_metrics(results, 200)
        rev_s = m['reversal_settling_s']
        rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
        print(f'{name:<18s} | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
    except Exception as e:
        print(f'{name:<18s} | ERROR: {str(e)[:50]}')
