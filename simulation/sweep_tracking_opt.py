#!/usr/bin/env python
"""
=======================================================================
  Future Work Item 5: Tracking Optimization for small_L @ CLBW=320
=======================================================================
Current: Tracking RMS = 4.28 rpm (strict threshold = 4 rpm)
Try: 1) different zeta  2) different FOC_delta  3) decoupling ON
"""
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy

p_base = get_motor_preset('small_L')

configs = [
    # (label,   zeta, FOC_delta, decoupling, CLBW)
    ('baseline',        15,   15,  False, 320),
    # -- zeta sweep --
    ('zeta=10',         10,   15,  False, 320),
    ('zeta=12',         12,   15,  False, 320),
    ('zeta=18',         18,   15,  False, 320),
    ('zeta=20',         20,   15,  False, 320),
    # -- FOC_delta sweep --
    ('delta=10',        15,   10,  False, 320),
    ('delta=20',        15,   20,  False, 320),
    ('delta=25',        15,   25,  False, 320),
    # -- decoupling ON --
    ('decouple_ON',     15,   15,  True,  320),
    ('decouple_z10',    10,   15,  True,  320),
    # -- slightly higher CLBW with adjusted zeta --
    ('CL350_z12',       12,   15,  False, 350),
    ('CL340_z13',       13,   15,  False, 340),
    ('CL330_z14',       14,   15,  False, 330),
]

header = f'{"Config":<16s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"TrkMax":>7s} {"Dip":>7s} {"Rev":>6s}'
print('small_L @ 20V, Stage 3 (full sensorless)')
print(header)
print('-' * 75)

for name, zeta, delta, decouple, clbw in configs:
    p = copy.deepcopy(p_base)
    p['d']['FOC_delta'] = delta
    p['d']['CTRL.bool_apply_decoupling_voltages_to_current_regulation'] = decouple
    try:
        results = run_sensorless_demo(
            p['d'], zeta=zeta, CLBW_Hz=clbw,
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
        print(f'{name:<16s} | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["trk_max_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
    except Exception as e:
        print(f'{name:<16s} | ERROR: {str(e)[:50]}')
