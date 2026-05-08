#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fine-tuned PLL speed S3 for servo:
1. PLL bw fine sweep: 20~200 Hz
2. zeta fine sweep: 0.5~5
3. CLBW: 30~200 Hz
4. Also test: PLL speed with additional 1st-order LPF smoothing
"""
import matplotlib; matplotlib.use('Agg')
from pylab import np, plt
import copy

from test_ab_speed import run_ab_speed_s3, PLLSpeedTracker
from demo_sensorless_active_flux import angle_error
from eval_staged_tuning import get_motor_preset

p = get_motor_preset('servo')
cmd = p['cmd_rpm']

configs = []
# Fine sweep around best regions
for bw in [30, 50, 80, 100, 120, 150, 200]:
    for z in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]:
        for cl in [30, 50, 80, 100, 150, 200]:
            configs.append(('pll', bw, 0, cl, z, f'bw{bw} z{z} CL{cl}'))

print(f'servo PLL speed S3 fine-tune ({len(configs)} configs)')
print(f'{"Label":<22s} | {"Ang":>6s} {"Trk":>6s} {"Qscore":>7s}')
print('-' * 48)

results = []
for method, bw, tau, cl, z, label in configs:
    try:
        r = run_ab_speed_s3(
            p['d'], CLBW_Hz=cl, zeta=z,
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            R_mis=1.0, v_off=(0,0),
            cmd_rpm=cmd, load_step=p['load_step'],
            speed_method='pll', speed_bw=bw, verbose=False)
        # Quality score: lower is better, penalize both angle and tracking
        q = r['ang_rms'] + 0.02 * r['trk_rms']
        results.append((label, r['ang_rms'], r['trk_rms'], q, r))
        if r['ang_rms'] < 15 or r['trk_rms'] < 100:
            flag = ' ★' if r['ang_rms']<5 and r['trk_rms']<50 else \
                   (' ◆' if r['ang_rms']<10 else (' ·' if r['ang_rms']<20 else ''))
            print(f'{label:<22s} | {r["ang_rms"]:>6.2f} {r["trk_rms"]:>6.1f} {q:>7.2f}{flag}')
    except:
        pass

# Sort by quality score
results.sort(key=lambda x: x[3])
print(f'\n--- Top 10 by quality score ---')
for i, (l, a, t, q, _) in enumerate(results[:10]):
    print(f'  {i+1:2d}. {l:<22s} | Ang={a:>6.2f} Trk={t:>6.1f} Q={q:>7.2f}')

# Sort by angle
results.sort(key=lambda x: x[1])
print(f'\n--- Top 10 by angle ---')
for i, (l, a, t, q, _) in enumerate(results[:10]):
    print(f'  {i+1:2d}. {l:<22s} | Ang={a:>6.2f} Trk={t:>6.1f}')

# Sort by tracking
results.sort(key=lambda x: x[2])
print(f'\n--- Top 10 by tracking ---')
for i, (l, a, t, q, _) in enumerate(results[:10]):
    print(f'  {i+1:2d}. {l:<22s} | Ang={a:>6.2f} Trk={t:>6.1f}')

# Plot best by quality score
results.sort(key=lambda x: x[3])
if results:
    best = results[0]
    r = best[4]
    fig, axes = plt.subplots(3,1,figsize=(14,9),sharex=True)
    t = r['t']
    axes[0].plot(t, r['cmd'],'r--',alpha=0.5,label='cmd')
    axes[0].plot(t, r['omega_true'],'b',label='true')
    axes[0].plot(t, r['omega_est'],'g',alpha=0.5,label='est')
    axes[0].set_ylabel('Speed [rpm]'); axes[0].legend()
    axes[0].set_title(f'servo S3 PLL | {best[0]} | Ang={best[1]:.2f} Trk={best[2]:.1f}')
    ae = np.degrees(angle_error(r['theta_true'], r['theta_af']))
    axes[1].plot(t, ae, 'r'); axes[1].set_ylabel('Angle err [deg]')
    axes[2].plot(t, r['omega_true']-r['cmd'],'b'); axes[2].set_ylabel('Track err [rpm]')
    axes[2].set_xlabel('Time [s]')
    plt.tight_layout(); plt.savefig('fig_pll_fine_best.png',dpi=150); plt.close()
    print('Saved fig_pll_fine_best.png')
    
    # Also plot the best tracking one
    results.sort(key=lambda x: x[2])
    best_trk = results[0]
    r = best_trk[4]
    fig, axes = plt.subplots(3,1,figsize=(14,9),sharex=True)
    t = r['t']
    axes[0].plot(t, r['cmd'],'r--',alpha=0.5,label='cmd')
    axes[0].plot(t, r['omega_true'],'b',label='true')
    axes[0].plot(t, r['omega_est'],'g',alpha=0.5,label='est')
    axes[0].set_ylabel('Speed [rpm]'); axes[0].legend()
    axes[0].set_title(f'servo S3 PLL | {best_trk[0]} | Ang={best_trk[1]:.2f} Trk={best_trk[2]:.1f}')
    ae = np.degrees(angle_error(r['theta_true'], r['theta_af']))
    axes[1].plot(t, ae, 'r'); axes[1].set_ylabel('Angle err [deg]')
    axes[2].plot(t, r['omega_true']-r['cmd'],'b'); axes[2].set_ylabel('Track err [rpm]')
    axes[2].set_xlabel('Time [s]')
    plt.tight_layout(); plt.savefig('fig_pll_fine_besttrk.png',dpi=150); plt.close()
    print('Saved fig_pll_fine_besttrk.png')
