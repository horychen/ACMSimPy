#!/usr/bin/env python
"""Fine-tune around CLBW=300-400 for small_L at 20V."""
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy

p_base = get_motor_preset('small_L')

clbws = [250, 280, 300, 320, 350, 380, 400, 450]

header = f'{"CLBW":>6s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print('DC_BUS_VOLTAGE = 20V, Stage 3 only')
print(header)
print('-' * 54)

for clbw in clbws:
    p = copy.deepcopy(p_base)
    p['d']['DC_BUS_VOLTAGE'] = 20
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
        print(f'{clbw:>6d} | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
    except Exception as e:
        print(f'{clbw:>6d} | ERROR: {str(e)[:50]}')
