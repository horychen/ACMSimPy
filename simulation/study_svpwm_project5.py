"""
Coding Project 5 - Zero-Sequence Voltage & Voltage Utilization Study

Sign convention:
  u_aG, u_bG, u_cG : terminal-to-ground potentials (0 or Vdc from inverter switching)
  u_nG : motor neutral-to-ground potential = (u_aG + u_bG + u_cG) / 3
       This is NOT zero-mean; it is biased around Vdc/2.

Two ways to compute u_nG:
  Method 1 (theoretical, from smooth commands):
    u_a* = u_alpha*,  u_b* = -u_alpha*/2 + sqrt3/2 * u_beta*,  u_c* = ...
    u_0  = -(max(u_a*,u_b*,u_c*) + min(u_a*,u_b*,u_c*)) / 2
    u_nG_th = Vdc/2 + u_0

  Method 2 (numerical, from switching terminal potentials):
    u_nG_num = (u_aG + u_bG + u_cG) / 3

Voltage utilization:
  eta_Vdc = (max(u_a*,u_b*,u_c*) - min(u_a*,u_b*,u_c*)) / 2
"""
from tutorials_ep6_svpwm import *
from collections import OrderedDict as OD

# ======================== Simulation Parameters ========================
d = {
    'CL_TS': 1e-4,
    'VL_EXE_PER_CL_EXE': 5,
    'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 500,
    'TIME_SLICE': 0.2,
    'NUMBER_OF_SLICES': 6,
    'init_npp': 4,
    'init_IN': 1.3*6/1.414,
    'init_R': 1.35,
    'init_Ld': 3.6e-3,
    'init_Lq': 3.6e-3,
    'init_KE': 0.125,
    'init_Rreq': 0.0,
    'init_Js': 0.44e-5,
    'DC_BUS_VOLTAGE': 48,
    'user_system_input_code': (
        "if ii < 1: CTRL.cmd_idq[0] = 0.0; CTRL.cmd_rpm = 50\n"
        "elif ii < 5: ACM.TLoad = 0.2\n"
        "elif ii < 100: CTRL.cmd_rpm = -50"
    ),
    'CTRL.bool_apply_speed_closed_loop_control': True,
    'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
    'CTRL.bool_apply_sweeping_frequency_excitation': False,
    'CTRL.bool_overwrite_speed_commands': True,
    'CTRL.bool_zero_id_control': True,
    'FOC_delta': 15,
    'FOC_desired_VLBW_HZ': 120,
    'FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False': 10,
    'CL_SERIES_KP': None, 'CL_SERIES_KI': None,
    'VL_SERIES_KP': None, 'VL_SERIES_KI': None,
    'VL_LIMIT_OVERLOAD_FACTOR': 3.0,
    'disp.Kp': 0.0, 'disp.Ki': 0.0, 'disp.Kd': 0.0,
    'disp.tau': 0.0, 'disp.OutLimit': 0.0, 'disp.IntLimit': 0.0,
}

Vdc = d['DC_BUS_VOLTAGE']
CL_TS = d['CL_TS']

# ======================== Run simulation (once) ========================
numba__scope_dict = OD([
    (r'Speed [rpm]',   ('CTRL.cmd_rpm', 'CTRL.omega_r_mech',)),
    (r'CTRL.uab [V]',  ('CTRL.cmd_uab[0]', 'CTRL.cmd_uab[1]',)),
    # NOTE: svgen1.S1/S2/S3 labels are misleading -- actual data stored is
    #   [30] = 30 + voltage_potential_at_terminal[0]
    #   [31] =      voltage_potential_at_terminal[1]
    #   [32] = -30 + voltage_potential_at_terminal[2]
    (r'S [1]',         ('svgen1.S1', 'svgen1.S2', 'svgen1.S3',
                         'svgen1.S4', 'svgen1.S5', 'svgen1.S6',)),
    (r'iab [A]',       ('CTRL.iab[0]', 'CTRL.iab[1]',)),
])

sim = Simulation_Benchmark(d, bool_start_simulation=False)
sim.start_simulation_slices(d, numba__scope_dict)
gdd = sim.gdd
t = sim.global_machine_times

# ======================== Extract signals ========================
# Undo display offsets to get true terminal-to-ground potentials
u_aG = gdd['svgen1.S1'] - 30.0
u_bG = gdd['svgen1.S2']
u_cG = gdd['svgen1.S3'] + 30.0

# ---- Method 2 (Numerical): u_nG from switching terminal potentials ----
u_nG_num = (u_aG + u_bG + u_cG) / 3.0

# ---- Method 1 (Theoretical): u_nG from smooth command voltages ----
# Inverse Clarke: reconstruct three-phase commands from alpha-beta
cmd_ua = gdd['CTRL.cmd_uab[0]']
cmd_ub = gdd['CTRL.cmd_uab[1]']
u_a_star = cmd_ua
u_b_star = -0.5 * cmd_ua + 0.8660254 * cmd_ub
u_c_star = -0.5 * cmd_ua - 0.8660254 * cmd_ub

u_max = np.maximum(np.maximum(u_a_star, u_b_star), u_c_star)
u_min = np.minimum(np.minimum(u_a_star, u_b_star), u_c_star)

# Injected zero-sequence (AC, zero-mean)
u_0_injected = -(u_max + u_min) / 2.0

# u_nG with Vdc/2 bias (this is the actual neutral-to-ground potential)
u_nG_th = Vdc / 2.0 + u_0_injected

# ---- Voltage utilization ----
eta_Vdc = (u_max - u_min) / 2.0

# ---- Phase currents from Clarke ----
iAlfa = gdd['CTRL.iab[0]']
iBeta = gdd['CTRL.iab[1]']
ia = iAlfa
ib = -0.5 * iAlfa + 0.8660254 * iBeta
ic = -0.5 * iAlfa - 0.8660254 * iBeta

# ======================== Zoom window ========================
t_center = 0.80  # well past initial transient
zoom_half = 2.5 * CL_TS
t_zs = t_center - zoom_half
t_ze = t_center + zoom_half
mz = (t >= t_zs) & (t <= t_ze)
tz = t[mz]

# ======================== Plotting ========================
plt.style.use('bmh')
mpl.rc('font', family='Times New Roman', size=10.0)
mpl.rc('legend', fontsize=7)
mpl.rcParams['lines.linewidth'] = 0.65
mpl.rcParams['mathtext.fontset'] = 'stix'

def add_sampling_lines(ax):
    for k in range(6):
        ax.axvline(t_zs + k * CL_TS, color='gray', ls='--', lw=0.35, alpha=0.5)

fig, axes = plt.subplots(5, 2, figsize=(14, 20), dpi=150, facecolor='w')

# ─────── F1: Terminal voltages ───────
ax = axes[0, 0]
ax.set_title('F1 - Terminal voltages (full view)', fontsize=9)
ax.plot(t, u_aG, alpha=0.5, label=r'$u_{aG}$')
ax.plot(t, u_bG, alpha=0.5, label=r'$u_{bG}$')
ax.plot(t, u_cG, alpha=0.5, label=r'$u_{cG}$')
ax.plot(t, cmd_ua, lw=1, label=r'$u_\alpha^*$ (cmd)', color='k', alpha=0.7)
ax.set_ylabel('Voltage [V]')
ax.legend(loc='best', fontsize=5)
ax.axvspan(t_zs, t_ze, alpha=0.25, color='red')

ax = axes[0, 1]
ax.set_title(f'F1 - Voltage zoom (5xCL_TS @ {t_center}s)', fontsize=9)
ax.plot(tz, u_aG[mz], label=r'$u_{aG}$', drawstyle='steps-post')
ax.plot(tz, u_bG[mz], label=r'$u_{bG}$', drawstyle='steps-post')
ax.plot(tz, u_cG[mz], label=r'$u_{cG}$', drawstyle='steps-post')
ax.plot(tz, cmd_ua[mz], lw=1.2, label=r'$u_\alpha^*$', color='k')
ax.plot(tz, cmd_ub[mz], lw=1.2, label=r'$u_\beta^*$', color='brown')
ax.set_ylabel('Voltage [V]')
ax.legend(loc='best', fontsize=5)
add_sampling_lines(ax)

# ─────── F2: Phase currents ───────
ax = axes[1, 0]
ax.set_title(r'F2 - Phase currents $i_a, i_b, i_c$ (full view)', fontsize=9)
ax.plot(t, ia, alpha=0.6, label=r'$i_a$')
ax.plot(t, ib, alpha=0.6, label=r'$i_b$')
ax.plot(t, ic, alpha=0.6, label=r'$i_c$')
ax.set_ylabel('Current [A]')
ax.legend(loc='best', fontsize=5)
ax.axvspan(t_zs, t_ze, alpha=0.25, color='red')

ax = axes[1, 1]
ax.set_title(f'F2 - Current zoom (5xCL_TS @ {t_center}s)', fontsize=9)
ax.plot(tz, ia[mz], label=r'$i_a$')
ax.plot(tz, ib[mz], label=r'$i_b$')
ax.plot(tz, ic[mz], label=r'$i_c$')
ax.set_ylabel('Current [A]')
ax.legend(loc='best', fontsize=5)
add_sampling_lines(ax)

# ─────── F3: u_nG comparison (both methods) ───────
ax = axes[2, 0]
ax.set_title(r'F3 - $u_{nG}$: Method 1 (theoretical) vs Method 2 (numerical)  [full]', fontsize=9)
ax.plot(t, u_nG_num, alpha=0.35, lw=0.4, color='tab:blue',
        label=r'Method 2: $u_{nG,\mathrm{num}} = (u_{aG}+u_{bG}+u_{cG})/3$')
ax.plot(t, u_nG_th, alpha=0.8, lw=0.8, color='tab:red',
        label=r'Method 1: $u_{nG,\mathrm{th}} = V_{dc}/2 - (\max+\min)/2$')
ax.axhline(Vdc/2, color='gray', ls=':', lw=0.5, alpha=0.5, label=r'$V_{dc}/2$')
ax.set_ylabel('Voltage [V]')
ax.legend(loc='best', fontsize=5)
ax.axvspan(t_zs, t_ze, alpha=0.25, color='red')

ax = axes[2, 1]
ax.set_title(f'F3 - $u_{{nG}}$ zoom (5xCL_TS @ {t_center}s)', fontsize=9)
ax.plot(tz, u_nG_num[mz], color='tab:blue', drawstyle='steps-post', lw=0.8,
        label=r'Method 2: $u_{nG,\mathrm{num}}$')
ax.plot(tz, u_nG_th[mz], color='tab:red', ls='--', drawstyle='steps-post', lw=1.2,
        label=r'Method 1: $u_{nG,\mathrm{th}}$')
ax.axhline(Vdc/2, color='gray', ls=':', lw=0.5, alpha=0.5, label=r'$V_{dc}/2$')
ax.set_ylabel('Voltage [V]')
ax.legend(loc='best', fontsize=5)
add_sampling_lines(ax)

# ─────── F4: Injected zero-sequence u_0 (the AC part, zero-mean) ───────
ax = axes[3, 0]
ax.set_title(r'F4 - Injected zero-sequence $u_0 = -(\max+\min)/2$ (full)', fontsize=9)
ax.plot(t, u_0_injected, color='purple', alpha=0.6, label=r'$u_0$ (injected)')
ax.plot(t, u_nG_num - Vdc/2, color='tab:blue', alpha=0.3, lw=0.4,
        label=r'$u_{nG,\mathrm{num}} - V_{dc}/2$ (AC part)')
ax.axhline(0, color='gray', ls=':', lw=0.5, alpha=0.5)
ax.set_ylabel('Voltage [V]')
ax.legend(loc='best', fontsize=5)
ax.axvspan(t_zs, t_ze, alpha=0.25, color='red')

ax = axes[3, 1]
ax.set_title(f'F4 - Injected $u_0$ zoom (5xCL_TS @ {t_center}s)', fontsize=9)
ax.plot(tz, u_0_injected[mz], color='purple', drawstyle='steps-post',
        label=r'$u_0 = -(\max+\min)/2$', lw=1.2)
ax.plot(tz, (u_nG_num - Vdc/2)[mz], color='tab:blue', drawstyle='steps-post',
        alpha=0.6, lw=0.7, label=r'$u_{nG,\mathrm{num}} - V_{dc}/2$')
ax.axhline(0, color='gray', ls=':', lw=0.5, alpha=0.5)
ax.set_ylabel('Voltage [V]')
ax.legend(loc='best', fontsize=5)
add_sampling_lines(ax)

# ─────── F5: Voltage utilization eta_Vdc ───────
ax = axes[4, 0]
ax.set_title(r'F5 - Voltage utilization $\eta_{V_{dc}} = (\max - \min)/2$ (full)', fontsize=9)
ax.plot(t, eta_Vdc, color='darkgreen', alpha=0.7, lw=0.7,
        label=r'$\eta_{V_{dc}} = (\max(u_a^*,u_b^*,u_c^*) - \min(...))/2$')
ax.axhline(Vdc/2, color='red', ls=':', lw=0.6, alpha=0.6,
           label=r'$V_{dc}/2 = %.0f$ V (SPWM limit)' % (Vdc/2))
ax.axhline(Vdc/np.sqrt(3), color='blue', ls=':', lw=0.6, alpha=0.6,
           label=r'$V_{dc}/\sqrt{3} = %.1f$ V (SVPWM limit)' % (Vdc/np.sqrt(3)))
ax.set_ylabel('Voltage [V]')
ax.set_xlabel('Time [s]')
ax.legend(loc='best', fontsize=5)
ax.axvspan(t_zs, t_ze, alpha=0.25, color='red')

ax = axes[4, 1]
ax.set_title(f'F5 - Voltage utilization zoom (5xCL_TS @ {t_center}s)', fontsize=9)
ax.plot(tz, eta_Vdc[mz], color='darkgreen', drawstyle='steps-post', lw=1.2,
        label=r'$\eta_{V_{dc}}$')
ax.axhline(Vdc/2, color='red', ls=':', lw=0.6, alpha=0.6,
           label=r'$V_{dc}/2$')
ax.axhline(Vdc/np.sqrt(3), color='blue', ls=':', lw=0.6, alpha=0.6,
           label=r'$V_{dc}/\sqrt{3}$')
ax.set_ylabel('Voltage [V]')
ax.set_xlabel('Time [s]')
ax.legend(loc='best', fontsize=5)
add_sampling_lines(ax)

for row in axes:
    for a in row:
        a.grid(True)

fig.tight_layout()
fig.savefig('codingProject5_zero_sequence.png', dpi=200, bbox_inches='tight')
try:
    fig.savefig('codingProject5_zero_sequence.pdf', dpi=400, bbox_inches='tight')
except PermissionError:
    fig.savefig('codingProject5_zero_sequence_v2.pdf', dpi=400, bbox_inches='tight')
    print("  (PDF saved as _v2.pdf due to file lock)")
print("\nDone! Saved: codingProject5_zero_sequence.png / .pdf")
plt.show()
