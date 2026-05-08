#!/usr/bin/env python
"""
Test servo motor Stage 3 with IDEAL parameters (no R mismatch, no voltage offsets).
Sweep CLBW to see if the AF observer can work under perfect conditions.
"""
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy

p_base = get_motor_preset('servo')

clbws = [50, 100, 200, 300, 500]
vbuses = [60, 100, 200]

header = f'{"Vbus":>4s} {"CLBW":>5s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print('servo motor, IDEAL params (R_mis=1.0, v_off=0)')
print(header)
print('-' * 58)

for vbus in vbuses:
    for clbw in clbws:
        p = copy.deepcopy(p_base)
        p['d']['DC_BUS_VOLTAGE'] = vbus
        try:
            results = run_sensorless_demo(
                p['d'], zeta=p['zeta'], CLBW_Hz=clbw,
                af_Kp=p['af_kp'], af_Ki=p['af_ki'],
                R_mismatch_factor=1.0,     # no mismatch
                voltage_offset_alpha=0.0,   # no offset
                voltage_offset_beta=0.0,
                eso_omega_ob=p.get('omega_ob', 80.0), speed_observer='nso',
                use_sensorless_speed=True,
                use_sensorless_angle=True,
                use_sensorless_torque_ff=True,
                cmd_rpm_ref=500, load_step=p['load_step'],
                verbose=False)
            m = compute_metrics(results, 500)
            rev_s = m['reversal_settling_s']
            rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
            print(f'{vbus:>4d} {clbw:>5d} | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
        except Exception as e:
            print(f'{vbus:>4d} {clbw:>5d} | ERROR: {str(e)[:50]}')
