"""
在 15~50 Hz 密集扫频，诊断蓝色 ESO 曲线在 disturbance 通道的幅值跌落。
"""
import numpy as np, copy, time, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import control

from sim_bode_sweep_v2 import MOTOR_PARAMS, run_single_frequency
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator

d = copy.deepcopy(MOTOR_PARAMS)
zeta, CLBW_Hz, omega_ob = 15, 1000, 200

# 密集扫频
freqs = np.linspace(15, 50, 40)
print(f"Sweeping {len(freqs)} points in [{freqs[0]:.1f}, {freqs[-1]:.1f}] Hz")
print(f"ESO bandwidth: {omega_ob/(2*np.pi):.1f} Hz")
print()

t0 = time.time()
sim_mag, sim_phase = [], []
for f in freqs:
    m, p = run_single_frequency(d, f, zeta, CLBW_Hz, True, omega_ob,
                                 mode='disturbance', rpm_0=500.0)
    sim_mag.append(m)
    sim_phase.append(p)
    print(f"  f={f:6.2f} Hz  mag={m:7.2f} dB  phase={p:7.1f} deg")

sim_mag = np.array(sim_mag)
sim_phase = np.array(sim_phase)
print(f"\nTotal simulation time: {time.time()-t0:.1f}s")

# ---- 解析对比 ----
R, L, J_s, n_pp, KE = d['init_R'], d['init_Lq'], d['init_Js'], d['init_npp'], d['init_KE']
CL_TS = d['CL_TS']
VL_TS = CL_TS * d['VL_EXE_PER_CL_EXE']
currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R, L, CLBW_Hz)
currentBandwidth = currentKp / L
speedKp, speedKi = get_coeffs_dc_motor_SPEED_regulator(J_s, n_pp, KE, zeta, currentBandwidth)

s = control.tf('s')
Kp_s = speedKp
Ki_s = speedKp * speedKi
G_pi_speed = Kp_s + Ki_s / s
G_current_cl = currentBandwidth / (s + currentBandwidth)
G_plant_speed = (n_pp * KE) / (J_s * s)
L_speed = G_pi_speed * G_current_cl * G_plant_speed
S_d = 1 / (1 + L_speed)

w = omega_ob
G_eso_residual = s**4 / (s + w)**4
G_dist_ESO = S_d * G_eso_residual * (-1 / (J_s * s)) * (30 / np.pi)

# 也算无 ESO 的
G_dist_noESO = S_d * (-1 / (J_s * s)) * (30 / np.pi)

freqs_dense = np.linspace(15, 50, 200)
ana_mag_eso, ana_mag_noeso = [], []
for f in freqs_dense:
    jw = 1j * 2 * np.pi * f
    ana_mag_eso.append(20 * np.log10(float(np.abs(G_dist_ESO(jw)))))
    ana_mag_noeso.append(20 * np.log10(float(np.abs(G_dist_noESO(jw)))))
ana_mag_eso = np.array(ana_mag_eso)
ana_mag_noeso = np.array(ana_mag_noeso)

# 也跑无 ESO 的仿真作为参照
sim_mag_noeso = []
for f in freqs:
    m, p = run_single_frequency(d, f, zeta, CLBW_Hz, False, 0,
                                 mode='disturbance', rpm_0=500.0)
    sim_mag_noeso.append(m)
    print(f"  [NoESO] f={f:6.2f} Hz  mag={m:7.2f} dB")
sim_mag_noeso = np.array(sim_mag_noeso)

# ---- 画图 ----
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Magnitude
ax1.plot(freqs_dense, ana_mag_noeso, 'g--', lw=1.5, alpha=0.5, label='Analytical (No ESO)')
ax1.plot(freqs_dense, ana_mag_eso, 'b--', lw=1.5, alpha=0.5, label='Analytical (ESO)')
ax1.plot(freqs, sim_mag_noeso, 'g^-', ms=4, lw=1, label='Sim (No ESO)')
ax1.plot(freqs, sim_mag, 'bo-', ms=4, lw=1, label='Sim (ESO)')
ax1.axvline(omega_ob/(2*np.pi), color='red', ls=':', lw=1, alpha=0.7, label=f'ESO BW = {omega_ob/(2*np.pi):.1f} Hz')
ax1.set_ylabel('Magnitude [dB]')
ax1.set_title(f'Disturbance Channel: Dense Sweep 15~50 Hz (zeta={zeta}, CLBW={CLBW_Hz}Hz, omega_ob={omega_ob})')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# Difference
ax2.plot(freqs, sim_mag - np.interp(freqs, freqs_dense, ana_mag_eso), 'bo-', ms=4, lw=1,
         label='Sim - Analytical (ESO)')
ax2.plot(freqs, sim_mag_noeso - np.interp(freqs, freqs_dense, ana_mag_noeso), 'g^-', ms=4, lw=1,
         label='Sim - Analytical (No ESO)')
ax2.axhline(0, color='k', ls='-', lw=0.5)
ax2.axvline(omega_ob/(2*np.pi), color='red', ls=':', lw=1, alpha=0.7)
ax2.set_ylabel('Sim - Analytical [dB]')
ax2.set_xlabel('Frequency [Hz]')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('fig_debug_eso_dip.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: fig_debug_eso_dip.png")
