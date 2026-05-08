#!/usr/bin/env python
"""
Sweep servo motor: simultaneous DC_BUS and CLBW optimization.
Find the sweet spot where bus voltage is high enough (no saturation) but low enough (good EMF/Bus).
"""
from demo_sensorless_active_flux import run_sensorless_demo
from eval_staged_tuning import get_motor_preset, compute_metrics
import copy, numpy as np

p_base = get_motor_preset('servo')
KE = p_base['d']['init_KE']
npp = p_base['d']['init_npp']
cmd_rpm = 500

omega_mech = cmd_rpm * 2 * np.pi / 60
emf_pk = KE * omega_mech * npp

# At 500 rpm: EMF_pk ≈ 20.9 V
# Need at least ~25V bus for some headroom (back-EMF + IR drop)
configs = [
    # (DC_BUS,  CLBW)
    (30,    500),
    (30,    300),
    (30,    200),
    (30,    100),
    (30,     50),
    (40,    500),
    (40,    300),
    (40,    200),
    (40,    100),
    (50,    500),
    (50,    300),
    (50,    200),
    (50,    100),
    (60,    500),
    (60,    300),
    (60,    200),
    (80,    500),
    (80,    300),
    (80,    200),
    (80,    100),
    (100,   500),
    (100,   300),
    (100,   200),
]

header = f'{"Vbus":>4s} {"CLBW":>5s} {"EMF%":>5s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print(f'servo motor, cmd_rpm={cmd_rpm}, EMF_pk={emf_pk:.1f}V')
print(header)
print('-' * 65)

for vbus, clbw in configs:
    p = copy.deepcopy(p_base)
    p['d']['DC_BUS_VOLTAGE'] = vbus
    emf_ratio = emf_pk / vbus * 100
    try:
        results = run_sensorless_demo(
            p['d'], zeta=p['zeta'], CLBW_Hz=clbw,
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mismatch_factor=p['r_mis'],
            voltage_offset_alpha=p['v_off'][0], voltage_offset_beta=p['v_off'][1],
            eso_omega_ob=p.get('omega_ob', 80.0), speed_observer='nso',
            use_sensorless_speed=True,
            use_sensorless_angle=True,
            use_sensorless_torque_ff=True,
            cmd_rpm_ref=cmd_rpm, load_step=p['load_step'],
            verbose=False)
        m = compute_metrics(results, cmd_rpm)
        rev_s = m['reversal_settling_s']
        rev = f'{rev_s:.3f}' if rev_s < 999 else 'FAIL'
        print(f'{vbus:>4d} {clbw:>5d} {emf_ratio:>4.0f}% | {m["angle_rms_ss_deg"]:>7.2f} {m["spd_est_rms_ss"]:>7.2f} {m["trk_rms_ss"]:>7.2f} {m["load_dip_rpm"]:>7.2f} {rev:>6s}')
    except Exception as e:
        print(f'{vbus:>4d} {clbw:>5d} {emf_ratio:>4.0f}% | ERROR: {str(e)[:50]}')
