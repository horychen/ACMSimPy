"""
验证实验：改变 omega_ob (100, 200, 400)，观察凹坑位置是否跟着 ESO 带宽移动。
如果是真实的离散干涉效应，凹坑应在 ~omega_ob/(2*pi) 附近。
"""
import numpy as np, copy, time, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sim_bode_sweep_v2 import MOTOR_PARAMS, run_single_frequency

d = copy.deepcopy(MOTOR_PARAMS)
zeta, CLBW_Hz = 15, 1000

# 三组 omega_ob
omega_obs = [100, 200, 400]
colors = ['#e74c3c', '#3498db', '#2ecc71']
labels = [r'$\omega_{ob}$=100 (15.9 Hz)', r'$\omega_{ob}$=200 (31.8 Hz)', r'$\omega_{ob}$=400 (63.7 Hz)']

freqs = np.linspace(5, 100, 50)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

for idx, (omega_ob, color, label) in enumerate(zip(omega_obs, colors, labels)):
    print(f'\n=== omega_ob = {omega_ob} (BW = {omega_ob/(2*np.pi):.1f} Hz) ===')
    sim_mag = []
    t0 = time.time()
    for f in freqs:
        m, p = run_single_frequency(d, f, zeta, CLBW_Hz, True, omega_ob,
                                     mode='disturbance', rpm_0=500.0)
        sim_mag.append(m)
    sim_mag = np.array(sim_mag)
    print(f'  Done in {time.time()-t0:.1f}s')

    ax1.plot(freqs, sim_mag, 'o-', color=color, ms=3, lw=1.2, label=label)
    ax1.axvline(omega_ob/(2*np.pi), color=color, ls=':', lw=1, alpha=0.6)

# Also plot no-ESO as reference
print('\n=== No ESO (reference) ===')
sim_mag_ref = []
for f in freqs:
    m, p = run_single_frequency(d, f, zeta, CLBW_Hz, False, 0,
                                 mode='disturbance', rpm_0=500.0)
    sim_mag_ref.append(m)
sim_mag_ref = np.array(sim_mag_ref)
ax1.plot(freqs, sim_mag_ref, 'k--', lw=1, alpha=0.5, label='No ESO (reference)')

ax1.set_ylabel('Magnitude [dB]')
ax1.set_title('Disturbance Channel: ESO Magnitude Dip vs Observer Bandwidth')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# Difference from no-ESO
for idx, (omega_ob, color, label) in enumerate(zip(omega_obs, colors, labels)):
    sim_mag = []
    for f in freqs:
        m, p = run_single_frequency(d, f, zeta, CLBW_Hz, True, omega_ob,
                                     mode='disturbance', rpm_0=500.0)
        sim_mag.append(m)
    sim_mag = np.array(sim_mag)
    ax2.plot(freqs, sim_mag - sim_mag_ref, 'o-', color=color, ms=3, lw=1.2, label=label)
    ax2.axvline(omega_ob/(2*np.pi), color=color, ls=':', lw=1, alpha=0.6)

ax2.axhline(0, color='k', ls='-', lw=0.5)
ax2.set_ylabel('ESO - NoESO [dB]')
ax2.set_xlabel('Frequency [Hz]')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('fig_verify_dip_vs_omega_ob.png', dpi=150, bbox_inches='tight')
print('\nSaved: fig_verify_dip_vs_omega_ob.png')
