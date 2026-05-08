#!/usr/bin/env python
"""
Full staged evaluation for small_L at 20V bus, CLBW=300Hz — the "no saturation" config.
"""
from demo_sensorless_active_flux import run_sensorless_demo, plot_sensorless_results
from eval_staged_tuning import get_motor_preset, compute_metrics, print_scorecard
import copy

p = copy.deepcopy(get_motor_preset('small_L'))
p['d']['DC_BUS_VOLTAGE'] = 20
clbw = 300

stages_cfg = [
    ('Stage0',  False, False, False),
    ('Stage2b', True,  False, False),
    ('Stage3',  True,  True,  True),
]

all_results = []
for name, use_spd, use_ang, use_tqff in stages_cfg:
    results = run_sensorless_demo(
        p['d'], zeta=p['zeta'], CLBW_Hz=clbw,
        af_Kp=p['af_kp'], af_Ki=p['af_ki'],
        R_mismatch_factor=p['r_mis'],
        voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
        eso_omega_ob=p.get('omega_ob', 200.0), speed_observer='nso',
        use_sensorless_speed=use_spd, use_sensorless_angle=use_ang,
        use_sensorless_torque_ff=use_tqff,
        cmd_rpm_ref=200, load_step=p['load_step'], verbose=True)
    plot_sensorless_results(results, save_path=f'fig_smallL_20V_CL300_{name}')
    import matplotlib.pyplot as plt
    plt.close('all')
    m = compute_metrics(results, 200)
    all_results.append((name, m))

print_scorecard(all_results, 200)
