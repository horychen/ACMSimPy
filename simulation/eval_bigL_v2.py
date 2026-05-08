#!/usr/bin/env python
"""big_L: Stage 3 with properly tuned CLBW (50-100 Hz range)."""
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics, print_scorecard
import copy

p_base = get_motor_preset('big_L')

# Sweep CLBW for Stage 3
clbws = [10, 15, 20, 30, 50, 60, 70, 80, 100]
header = f'{"CLBW":>6s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print('big_L Stage 3, omega_ob=200')
print(header)
print('-' * 55)

for clbw in clbws:
    p = copy.deepcopy(p_base)
    try:
        results = run_sensorless_demo(
            p['d'], zeta=p['zeta'], CLBW_Hz=clbw,
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mismatch_factor=p['r_mis'],
            voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
            eso_omega_ob=200.0, speed_observer='nso',
            use_sensorless_speed=True, use_sensorless_angle=True,
            use_sensorless_torque_ff=True,
            cmd_rpm_ref=150, load_step=p['load_step'], verbose=False)
        m = compute_metrics(results, 150)
        rev_s = m['reversal_settling_s']
        rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
        print(f'{clbw:>6d} | {m["angle_rms_ss_deg"]:>7.02f} {m["spd_est_rms_ss"]:>7.02f} {m["trk_rms_ss"]:>7.02f} {m["load_dip_rpm"]:>7.02f} {rev:>6s}')
    except Exception as e:
        print(f'{clbw:>6d} | ERROR: {str(e)[:50]}')

# Also try lower omega_ob
print('\nbig_L Stage 3, CLBW=80, sweep omega_ob')
obs_bws = [10, 20, 30, 50, 80, 100, 150, 200]
header = f'{"w_ob":>6s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print(header)
print('-' * 55)

for w_ob in obs_bws:
    p = copy.deepcopy(p_base)
    try:
        results = run_sensorless_demo(
            p['d'], zeta=p['zeta'], CLBW_Hz=80,
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mismatch_factor=p['r_mis'],
            voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
            eso_omega_ob=w_ob, speed_observer='nso',
            use_sensorless_speed=True, use_sensorless_angle=True,
            use_sensorless_torque_ff=True,
            cmd_rpm_ref=150, load_step=p['load_step'], verbose=False)
        m = compute_metrics(results, 150)
        rev_s = m['reversal_settling_s']
        rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
        print(f'{w_ob:>6d} | {m["angle_rms_ss_deg"]:>7.02f} {m["spd_est_rms_ss"]:>7.02f} {m["trk_rms_ss"]:>7.02f} {m["load_dip_rpm"]:>7.02f} {rev:>6s}')
    except Exception as e:
        print(f'{w_ob:>6d} | ERROR: {str(e)[:50]}')

# Full eval with best config
print('\n\nbig_L FULL staged eval: CLBW=80, omega_ob=80')
stages_cfg = [
    ('Stage0',  False, False, False),
    ('Stage2b', True,  False, False),
    ('Stage3',  True,  True,  True),
]
all_results = []
for name, use_spd, use_ang, use_tqff in stages_cfg:
    p = copy.deepcopy(p_base)
    results = run_sensorless_demo(
        p['d'], zeta=p['zeta'], CLBW_Hz=80,
        af_Kp=p['af_kp'], af_Ki=p['af_ki'],
        R_mismatch_factor=p['r_mis'],
        voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
        eso_omega_ob=80.0, speed_observer='nso',
        use_sensorless_speed=use_spd, use_sensorless_angle=use_ang,
        use_sensorless_torque_ff=use_tqff,
        cmd_rpm_ref=150, load_step=p['load_step'], verbose=True)
    m = compute_metrics(results, 150)
    all_results.append((name, m))

print_scorecard(all_results, 150)
