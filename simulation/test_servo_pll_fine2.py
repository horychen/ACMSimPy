#!/usr/bin/env python
"""Fine-tune servo PLL around CL=30, pll=50 (ideal) with different zeta & omega_ob."""
import sys; sys.path.insert(0, '.')
from test_alt_observers import run_alt_observer_sim, eval_results
from eval_staged_tuning import get_motor_preset

p = get_motor_preset('servo')
cmd = p['cmd_rpm']

configs = [
    # (CLBW, pll_bw, zeta, omega_ob, label)
    # Best region from previous: CL=20-40, pll=30-80
    (20, 40, 5, 30, 'z5 w30'),
    (20, 40, 10, 30, 'z10 w30'),
    (25, 40, 5, 30, 'z5 CL25'),
    (25, 50, 5, 30, 'z5 p50'),
    (25, 50, 8, 30, 'z8 p50'),
    (30, 40, 5, 30, 'CL30 p40 z5'),
    (30, 50, 5, 30, 'CL30 p50 z5'),
    (30, 50, 5, 50, 'CL30 p50 z5 w50'),
    (30, 50, 8, 30, 'CL30 p50 z8'),
    (30, 50, 10, 30, 'CL30 p50 z10'),
    (30, 60, 5, 30, 'CL30 p60 z5'),
    (30, 70, 5, 30, 'CL30 p70 z5'),
    (35, 50, 5, 30, 'CL35 p50 z5'),
    (35, 60, 5, 30, 'CL35 p60 z5'),
    (40, 50, 5, 30, 'CL40 p50 z5'),
    (40, 60, 5, 30, 'CL40 p60 z5'),
    (40, 80, 5, 30, 'CL40 p80 z5'),
    # Very low zeta
    (30, 50, 3, 30, 'CL30 p50 z3'),
    (30, 50, 2, 30, 'CL30 p50 z2'),
    (25, 40, 3, 30, 'CL25 p40 z3'),
]

header = f'{"Label":<18s} | {"AngSS":>7s} {"TrkSS":>7s} {"Dip":>7s} {"Rev":>6s}'
print('servo PLL fine-tune (IDEAL params, R_mis=1.0, v_off=0)')
print(header); print('-' * 55)

for clbw, pll_bw, zeta, w_ob, label in configs:
    try:
        r = run_alt_observer_sim(
            p['d'], observer_type='pll', CLBW_Hz=clbw, zeta=zeta,
            R_mismatch_factor=1.0, v_off=(0, 0),
            omega_ob=w_ob, cmd_rpm=cmd, load_step=p['load_step'],
            pll_bw=pll_bw, af_Kp=p['af_kp'], af_Ki=p['af_ki'], verbose=False)
        ang, spd, trk, dip, rev = eval_results(r, cmd)
        rev_s = f'{rev:.3f}' if rev < 999 else 'FAIL'
        flag = ' ★' if ang < 10 else (' ◆' if ang < 30 else (' ·' if ang < 60 else ''))
        print(f'{label:<18s} | {ang:>7.2f} {trk:>7.2f} {dip:>7.2f} {rev_s:>6s}{flag}')
    except Exception as e:
        print(f'{label:<18s} | ERROR: {str(e)[:40]}')
