#!/usr/bin/env python
"""
Plot: AF Gain Sweep for all three motors - OE vs Angle Error vs Gain
"""
import matplotlib; matplotlib.use('Agg')
from pylab import np, plt
from diag_openloop import run_open_loop_diag
from demo_sensorless_active_flux import angle_error
from eval_staged_tuning import get_motor_preset

gains = [
    (0, 0), (10, 100), (20, 200), (50, 500), (100, 1000),
    (200, 2000), (500, 5000), (1000, 10000), (2000, 20000),
    (5000, 50000), (10000, 100000),
]

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Open-Loop AF Observer: OE vs Angle Error vs Gain (Ideal Params)', fontsize=14, fontweight='bold')

for col, motor_name in enumerate(['servo', 'small_L', 'big_L']):
    p = get_motor_preset(motor_name)
    cmd = p['cmd_rpm']
    KE = p['d']['init_KE']
    
    kps = []; ang_list = []; oe_list = []; corr_list = []
    
    for kp, ki in gains:
        r = run_open_loop_diag(
            p['d'], CLBW_Hz=p['clbw'], zeta=p['zeta'],
            R_mismatch_factor=1.0, v_off=(0,0),
            af_Kp=kp, af_Ki=ki, pll_bw=100,
            cmd_rpm=cmd, load_step=p['load_step'], verbose=False)
        ss = r['t'] > 0.4
        ae = np.degrees(angle_error(r['theta_true'], r['af_theta']))
        ang_rms = np.sqrt(np.mean(ae[ss]**2))
        oe_rms = np.sqrt(np.mean(r['af_oe'][ss]**2)) * 1e3  # mWb
        corr_rms = np.sqrt(np.mean(r['af_corr_mag'][ss]**2))
        kps.append(kp if kp > 0 else 0.1)
        ang_list.append(ang_rms)
        oe_list.append(oe_rms)
        corr_list.append(corr_rms)

    kps = np.array(kps)
    # Top row: Angle error & OE vs Kp
    ax1 = axes[0, col]
    ax1.semilogx(kps, ang_list, 'ro-', label='Angle RMS (deg)', markersize=5)
    ax1.set_ylabel('Angle RMS [deg]', color='r')
    ax1.tick_params(axis='y', labelcolor='r')
    ax1r = ax1.twinx()
    ax1r.semilogx(kps, oe_list, 'bs-', label='OE RMS (mWb)', markersize=5)
    ax1r.set_ylabel('OE RMS [mWb]', color='b')
    ax1r.tick_params(axis='y', labelcolor='b')
    ax1.set_title(f'{motor_name} (KE={KE*1e3:.1f} mWb)')
    ax1.set_xlabel('Kp')

    # Bottom row: Correction magnitude vs Kp
    ax2 = axes[1, col]
    ax2.semilogx(kps, corr_list, 'g^-', label='|Correction| RMS (V)', markersize=5)
    ax2.set_ylabel('Correction RMS [V]', color='g')
    ax2.set_xlabel('Kp')
    # Also show OE
    ax2r = ax2.twinx()
    ax2r.semilogx(kps, oe_list, 'bs-', alpha=0.5, markersize=4)
    ax2r.set_ylabel('OE RMS [mWb]', color='b')

plt.tight_layout()
plt.savefig('fig_gain_sweep_3motors.png', dpi=150)
plt.close()
print('Saved fig_gain_sweep_3motors.png')
