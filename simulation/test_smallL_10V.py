#!/usr/bin/env python
"""Test small_L with doubled DC bus voltage (5V -> 10V)."""
from demo_sensorless_active_flux import run_sensorless_demo, plot_sensorless_results
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy

p = copy.deepcopy(get_motor_preset('small_L'))
p['d']['DC_BUS_VOLTAGE'] = 10  # 5V -> 10V

stages_cfg = [
    ('S0_10V',  False, False, False),
    ('S2b_10V', True,  False, False),
    ('S3_10V',  True,  True,  True),
]

for name, use_spd, use_ang, use_tqff in stages_cfg:
    results = run_sensorless_demo(
        p['d'], zeta=p['zeta'], CLBW_Hz=p['clbw'],
        af_Kp=p['af_kp'], af_Ki=p['af_ki'],
        R_mismatch_factor=p['r_mis'],
        voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
        eso_omega_ob=p.get('omega_ob', 200.0), speed_observer='nso',
        use_sensorless_speed=use_spd, use_sensorless_angle=use_ang,
        use_sensorless_torque_ff=use_tqff,
        cmd_rpm_ref=p['cmd_rpm'], load_step=p['load_step'], verbose=False)
    plot_sensorless_results(results, save_path=f'fig_smallL_10V_{name}')
    import matplotlib.pyplot as plt
    plt.close('all')
    m = compute_metrics(results, p['cmd_rpm'])
    rev_s = m['reversal_settling_s']
    rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
    ang = m['angle_rms_ss_deg']
    spd = m['spd_est_rms_ss']
    trk = m['trk_rms_ss']
    dip = m['load_dip_rpm']
    print(f'{name:<12s} | AngSS={ang:>6.2f} SpdSS={spd:>6.2f} TrkSS={trk:>6.2f} Dip={dip:>6.2f} Rev={rev}')
