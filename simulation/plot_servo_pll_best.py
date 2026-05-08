#!/usr/bin/env python
"""Plot the best servo PLL result to understand what 30deg means."""
import matplotlib; matplotlib.use('Agg')
from pylab import np, plt
from test_alt_observers import run_alt_observer_sim, eval_results
from demo_sensorless_active_flux import angle_error
from eval_staged_tuning import get_motor_preset

p = get_motor_preset('servo')
res = run_alt_observer_sim(
    p['d'], observer_type='pll', CLBW_Hz=30, zeta=2,
    R_mismatch_factor=1.0, v_off=(0,0),
    omega_ob=30, cmd_rpm=500, load_step=p['load_step'],
    pll_bw=50, af_Kp=p['af_kp'], af_Ki=p['af_ki'], verbose=True)

ang, spd, trk, dip, rev = eval_results(res, 500)
print(f'Angle={ang:.1f}, Tracking={trk:.1f}, Dip={dip:.1f}')

t = res['t']
fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

axes[0].plot(t, res['cmd_rpm'], 'r--', label='cmd', alpha=0.5)
axes[0].plot(t, res['omega_true'], 'b', label='true')
axes[0].plot(t, res['omega_est'], 'g', alpha=0.6, label='est (NSO)')
axes[0].set_ylabel('Speed [rpm]'); axes[0].legend(loc='upper right')
axes[0].set_title(f'servo PLL CL=30 pll=50 z=2 (ideal) | Ang={ang:.1f} deg')

ae = np.degrees(angle_error(res['theta_true'], res['theta_est']))
axes[1].plot(t, ae, 'r')
axes[1].axhline(0, color='k', lw=0.5)
axes[1].set_ylabel('Angle error [deg]')

axes[2].plot(t, res['omega_true'] - res['cmd_rpm'], 'b')
axes[2].axhline(0, color='k', lw=0.5)
axes[2].set_ylabel('Track err [rpm]'); axes[2].set_xlabel('Time [s]')

plt.tight_layout()
plt.savefig('fig_servo_pll_best_v2.png', dpi=150)
plt.close()
print('Saved fig_servo_pll_best_v2.png')
