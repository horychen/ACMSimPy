#!/usr/bin/env python
"""
Test small_L with high bus voltage (no saturation) at different cmd_rpm.
Goal: find cmd_rpm where EMF/Bus ratio is high enough for stable Stage 3.
"""
from demo_sensorless_active_flux import run_sensorless_demo, angle_error
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy, numpy as np

p = copy.deepcopy(get_motor_preset('small_L'))
# DC bus already set to 150V in preset (no saturation)

KE = p['d']['init_KE']
npp = p['d']['init_npp']
Vbus = p['d']['DC_BUS_VOLTAGE']

rpms = [200, 500, 1000, 1500, 2000, 3000, 5000]

header = f'{"RPM":>6s} {"EMF_pk":>7s} {"EMF/Bus":>7s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print(header)
print('-' * 72)

for rpm in rpms:
    omega_mech = rpm * 2 * np.pi / 60
    emf_pk = KE * omega_mech * npp
    emf_ratio = emf_pk / Vbus * 100

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
            cmd_rpm_ref=rpm, load_step=p['load_step'],
            verbose=False)
        m = compute_metrics(results, rpm)
        rev_s = m['reversal_settling_s']
        rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
        print(f'{rpm:>6d} {emf_pk:>7.2f} {emf_ratio:>6.1f}% | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
    except Exception as e:
        print(f'{rpm:>6d} {emf_pk:>7.2f} {emf_ratio:>6.1f}% | ERROR: {str(e)[:40]}')

print(f'\nVbus={Vbus}V, KE={KE}, npp={npp}')
