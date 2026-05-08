#!/usr/bin/env python
"""
=======================================================================
  Future Work Item 1: big_L Motor Full Pipeline
=======================================================================
Step 1: Run baseline Stage 0/2b/3 with current preset
Step 2: If S3 fails, sweep CLBW to find stable region  
Step 3: Sweep NSO omega_ob to optimize
"""
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics, print_scorecard
import copy, numpy as np

# ---- big_L motor parameters ----
# npp=24, R=1.97Ω, L=103.5mH, KE=0.0745Wb, J=0.0765
# L/R = 52.5 ms (very slow!), KE/L = 0.72 (very low!)
# DC_BUS = 800V, cmd_rpm = 150

p = get_motor_preset('big_L')
cmd_rpm = p['cmd_rpm']
KE = p['d']['init_KE']
npp = p['d']['init_npp']
Vbus = p['d']['DC_BUS_VOLTAGE']

omega_mech = cmd_rpm * 2 * np.pi / 60
emf_pk = KE * omega_mech * npp
print(f'big_L: EMF_pk={emf_pk:.1f}V, EMF/Bus={emf_pk/Vbus*100:.1f}%, L/R={p["d"]["init_Ld"]/p["d"]["init_R"]*1e3:.1f}ms, KE/L={KE/p["d"]["init_Ld"]:.1f}')

# =====================================================
print('\n' + '='*72)
print('  STEP 1: Baseline eval with preset params')
print('='*72)

stages_cfg = [
    ('Stage0',  False, False, False),
    ('Stage2b', True,  False, False),
    ('Stage3',  True,  True,  True),
]

all_results = []
for name, use_spd, use_ang, use_tqff in stages_cfg:
    print(f'\n--- {name} ---')
    try:
        results = run_sensorless_demo(
            p['d'], zeta=p['zeta'], CLBW_Hz=p['clbw'],
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mismatch_factor=p['r_mis'],
            voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
            eso_omega_ob=p.get('omega_ob', 200.0), speed_observer='nso',
            use_sensorless_speed=use_spd, use_sensorless_angle=use_ang,
            use_sensorless_torque_ff=use_tqff,
            cmd_rpm_ref=cmd_rpm, load_step=p['load_step'], verbose=True)
        m = compute_metrics(results, cmd_rpm)
        all_results.append((name, m))
    except Exception as e:
        print(f'  ERROR: {e}')
        all_results.append((name, None))

print_scorecard(all_results, cmd_rpm)

# =====================================================
print('\n' + '='*72)
print('  STEP 2: CLBW sweep for Stage 3')
print('='*72)

clbws = [10, 20, 30, 50, 80, 100, 150, 200, 300]
header = f'{"CLBW":>6s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print(header)
print('-' * 55)

for clbw in clbws:
    p2 = copy.deepcopy(p)
    try:
        results = run_sensorless_demo(
            p2['d'], zeta=p2['zeta'], CLBW_Hz=clbw,
            af_Kp=p2['af_kp'], af_Ki=p2['af_ki'],
            R_mismatch_factor=p2['r_mis'],
            voltage_offset_alpha=p2['v_off'][0], voltage_offset_beta=p2['v_off'][1],
            eso_omega_ob=p2.get('omega_ob', 200.0), speed_observer='nso',
            use_sensorless_speed=True,
            use_sensorless_angle=True,
            use_sensorless_torque_ff=True,
            cmd_rpm_ref=cmd_rpm, load_step=p2['load_step'],
            verbose=False)
        m = compute_metrics(results, cmd_rpm)
        rev_s = m['reversal_settling_s']
        rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
        print(f'{clbw:>6d} | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
    except Exception as e:
        print(f'{clbw:>6d} | ERROR: {str(e)[:50]}')

# =====================================================
print('\n' + '='*72)
print('  STEP 3: NSO omega_ob sweep for Stage 3 (using best CLBW from preset)')
print('='*72)

obs_bws = [10, 20, 50, 80, 100, 150, 200, 300]
header = f'{"w_ob":>6s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print(header)
print('-' * 55)

for w_ob in obs_bws:
    p2 = copy.deepcopy(p)
    try:
        results = run_sensorless_demo(
            p2['d'], zeta=p2['zeta'], CLBW_Hz=p2['clbw'],
            af_Kp=p2['af_kp'], af_Ki=p2['af_ki'],
            R_mismatch_factor=p2['r_mis'],
            voltage_offset_alpha=p2['v_off'][0], voltage_offset_beta=p2['v_off'][1],
            eso_omega_ob=w_ob, speed_observer='nso',
            use_sensorless_speed=True,
            use_sensorless_angle=True,
            use_sensorless_torque_ff=True,
            cmd_rpm_ref=cmd_rpm, load_step=p2['load_step'],
            verbose=False)
        m = compute_metrics(results, cmd_rpm)
        rev_s = m['reversal_settling_s']
        rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
        print(f'{w_ob:>6d} | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.02f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
    except Exception as e:
        print(f'{w_ob:>6d} | ERROR: {str(e)[:50]}')
