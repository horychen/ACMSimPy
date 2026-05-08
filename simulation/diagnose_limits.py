#!/usr/bin/env python
"""
Find the minimal DC_BUS_VOLTAGE where small_L motor doesn't saturate at 200 rpm.
Use moderate VL_LIMIT_OVERLOAD_FACTOR (3, 5, 10) to avoid numerical blow-up.
"""
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy, numpy as np

p_base = get_motor_preset('small_L')

# Restore VL factor to reasonable values
configs = [
    # (label,   DC_BUS,  VL_OVERLOAD, cmd_rpm)
    ('5V_3x',      5,    3.0,  200),   # original baseline
    ('10V_3x',    10,    3.0,  200),   # 2x bus
    ('15V_3x',    15,    3.0,  200),   # 3x bus
    ('20V_3x',    20,    3.0,  200),   # 4x bus
    ('10V_5x',    10,    5.0,  200),   # 2x bus, higher IL
    ('20V_5x',    20,    5.0,  200),   # 4x bus, higher IL
    ('10V_10x',   10,   10.0,  200),   # 2x bus, 10x IL
    ('20V_10x',   20,   10.0,  200),   # 4x bus, 10x IL
]

header = f'{"Config":<12s} {"Vbus":>4s} {"IL_lim":>6s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print(header)
print('-' * 68)

for name, vbus, vl_factor, cmd_rpm in configs:
    p = copy.deepcopy(p_base)
    p['d']['DC_BUS_VOLTAGE'] = vbus
    p['d']['VL_LIMIT_OVERLOAD_FACTOR'] = vl_factor
    il_max = vl_factor * 1.414 * p['d']['init_IN']
    try:
        results = run_sensorless_demo(
            p['d'], zeta=p['zeta'], CLBW_Hz=p['clbw'],
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mismatch_factor=p['r_mis'],
            voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
            eso_omega_ob=p.get('omega_ob', 200.0), speed_observer='nso',
            use_sensorless_speed=True,
            use_sensorless_angle=True,
            use_sensorless_torque_ff=True,
            cmd_rpm_ref=cmd_rpm, load_step=p['load_step'],
            verbose=False)
        m = compute_metrics(results, cmd_rpm)
        rev_s = m['reversal_settling_s']
        rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
        print(f'{name:<12s} {vbus:>4d} {il_max:>6.1f} | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
    except Exception as e:
        print(f'{name:<12s} {vbus:>4d} {il_max:>6.1f} | ERROR: {str(e)[:50]}')
