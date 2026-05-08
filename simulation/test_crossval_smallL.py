#!/usr/bin/env python
"""Cross-validate: PLL and SMO on small_L (should work since AF works)."""
import sys; sys.path.insert(0, '.')
from test_alt_observers import run_alt_observer_sim, eval_results
from eval_staged_tuning import get_motor_preset
import copy

p = get_motor_preset('small_L')
cmd = p['cmd_rpm']

configs = [
    # (obs, CLBW, pll_bw, smo_gain, smo_lpf, omega_ob, zeta, label)
    ('pll', 320, 100, None, None, 200, 15, 'PLL bw100 CL320'),
    ('pll', 320, 200, None, None, 200, 15, 'PLL bw200 CL320'),
    ('pll', 320, 500, None, None, 200, 15, 'PLL bw500 CL320'),
    ('pll', 320, 1000, None, None, 200, 15, 'PLL bw1000 CL320'),
    ('smo', 320, 100, None, 500, 200, 15, 'SMO lpf500 CL320'),
    ('smo', 320, 200, None, 1000, 200, 15, 'SMO lpf1k CL320'),
    ('smo', 320, 500, None, 2000, 200, 15, 'SMO lpf2k CL320'),
]

header = f'{"Label":<22s} | {"AngSS":>7s} {"SpdSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print('small_L cross-validation: PLL & SMO')
print(header); print('-' * 68)

for obs, clbw, pll_bw, sg, sl, w_ob, z, label in configs:
    try:
        res = run_alt_observer_sim(
            p['d'], observer_type=obs, CLBW_Hz=clbw, zeta=z,
            R_mismatch_factor=p['r_mis'], v_off=p['v_off'],
            omega_ob=w_ob, cmd_rpm=cmd, load_step=p['load_step'],
            pll_bw=pll_bw, smo_gain=sg, smo_lpf=sl,
            af_Kp=p['af_kp'], af_Ki=p['af_ki'], verbose=False)
        ang, spd, trk, dip, rev = eval_results(res, cmd)
        rev_s = f'{rev:.3f}' if rev < 999 else 'FAIL'
        print(f'{label:<22s} | {ang:>7.2f} {spd:>7.2f} {trk:>7.2f} {dip:>7.2f} {rev_s:>6s}')
    except Exception as e:
        print(f'{label:<22s} | ERROR: {str(e)[:45]}')
