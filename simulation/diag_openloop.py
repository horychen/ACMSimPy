#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
=======================================================================
  Open-Loop Observer Diagnostic
=======================================================================
Stage 0 control (encoder angle + encoder speed), all observers run in
parallel as PURE OBSERVERS (no feedback to controller).

For each observer, record:
  - angle estimation error vs true encoder angle
  - flux amplitude vs expected KE
  - AF PI correction magnitude (output error OE)
  - effect of observer gains (af_Kp, af_Ki) on OE and angle error

This isolates observer performance from control loop feedback.
"""
import matplotlib
matplotlib.use('Agg')
from pylab import np, plt
import copy, time as _time

from tutorials_ep6_svpwm import (
    The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental,
)
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator
from demo_sensorless_active_flux import wrap_angle, angle_error, angle_diff_scalar
from eval_staged_tuning import get_motor_preset
from observers_alt import PLLFluxObserver, SlidingModeObserver


def run_open_loop_diag(d, CLBW_Hz=100, zeta=10,
                       R_mismatch_factor=1.0, v_off=(0.0, 0.0),
                       af_Kp=500.0, af_Ki=5000.0,
                       pll_bw=100.0, smo_gain=None, smo_lpf=200.0,
                       cmd_rpm=500, load_step=1.27, verbose=True):
    """
    Run Stage 0 (encoder control) with AF, PLL, SMO observers in parallel.
    Returns detailed per-step diagnostics.
    """
    dd = copy.deepcopy(d)
    R_true = dd['init_R']
    R_obs = R_true * R_mismatch_factor
    L = dd['init_Lq']
    J_s = dd['init_Js']
    n_pp = dd['init_npp']
    KE = dd['init_KE']
    CL_TS = dd['CL_TS']
    VL_TS = CL_TS * dd['VL_EXE_PER_CL_EXE']

    cKp, cKi = get_coeffs_dc_motor_current_regulator(R_true, L, CLBW_Hz)
    sKp, sKi = get_coeffs_dc_motor_SPEED_regulator(J_s, n_pp, KE, zeta, cKp / L)
    dd['CL_SERIES_KP'] = cKp; dd['CL_SERIES_KI'] = cKi
    dd['VL_SERIES_KP'] = sKp; dd['VL_SERIES_KI'] = sKi

    # Stage 0: encoder control, no sensorless injection
    CTRL = The_Motor_Controller(CL_TS=CL_TS, VL_TS=VL_TS, init_npp=n_pp,
        init_IN=dd['init_IN'], init_R=R_true, init_Ld=dd['init_Ld'],
        init_Lq=L, init_KE=KE, init_Rreq=0.0, init_Js=J_s,
        DC_BUS_VOLTAGE=dd['DC_BUS_VOLTAGE'])
    CTRL.bool_apply_decoupling_voltages_to_current_regulation = False
    CTRL.bool_apply_sweeping_frequency_excitation = False
    CTRL.bool_overwrite_speed_commands = True
    CTRL.bool_zero_id_control = True
    CTRL.bool_apply_speed_closed_loop_control = True

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1)
    ki_f = dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False', 10)
    vl_lim = dd.get('VL_LIMIT_OVERLOAD_FACTOR', 10.0)
    Vm = dd['DC_BUS_VOLTAGE'] / 1.732
    Im = vl_lim * 1.414 * dd['init_IN']
    reg_id = The_PID_Regulator(cKp, cKp*cKi*ki_f, 0,0, Vm, Vm, CL_TS)
    reg_iq = The_PID_Regulator(cKp, cKp*cKi*ki_f, 0,0, Vm, Vm, CL_TS)
    reg_spd = The_PID_Regulator(sKp, sKp*sKi, 0,0, Im, Im, VL_TS)

    # ---- Observer instances (all run open-loop, no feedback) ----
    # AF observer state
    af_psi_s = np.array([KE, 0.0])
    af_corr_int = np.zeros(2)

    # PLL observer
    pll_obs = PLLFluxObserver(R_obs, L, KE, CL_TS, pll_bw=pll_bw, af_Kp=af_Kp, af_Ki=af_Ki)

    # SMO observer
    smo_obs = SlidingModeObserver(R_obs, L, KE, CL_TS, smo_gain=smo_gain, lpf_fc=smo_lpf, pll_bw=pll_bw)

    # ---- Storage ----
    total_time = dd['TIME_SLICE'] * dd['NUMBER_OF_SLICES']
    n_ctrl = int(total_time / CL_TS) + 100

    t_tr = np.zeros(n_ctrl)
    theta_true_arr = np.zeros(n_ctrl)
    omega_true_arr = np.zeros(n_ctrl)

    # AF diagnostics
    af_theta = np.zeros(n_ctrl)
    af_psi_amp = np.zeros(n_ctrl)        # |ψ_AF|
    af_oe = np.zeros(n_ctrl)             # output error = KE - |ψ_AF|
    af_corr_mag = np.zeros(n_ctrl)       # |correction| magnitude
    af_corr_int_mag = np.zeros(n_ctrl)   # |integrator| magnitude

    # PLL diagnostics
    pll_theta = np.zeros(n_ctrl)
    pll_omega = np.zeros(n_ctrl)

    # SMO diagnostics
    smo_theta = np.zeros(n_ctrl)
    smo_omega = np.zeros(n_ctrl)

    # Voltage/current info
    u_mag_arr = np.zeros(n_ctrl)
    i_mag_arr = np.zeros(n_ctrl)
    emf_true_arr = np.zeros(n_ctrl)

    CTRL.cmd_rpm = 0.0
    ci = 0
    t0_sim = _time.time()

    sim_TS = dd['TIME_SLICE']
    sim_N = dd['NUMBER_OF_SLICES']

    for sl in range(sim_N):
        t0 = sl * sim_TS
        tm = t0 + sim_TS * 0.5

        if tm < 0.3:
            CTRL.cmd_rpm = cmd_rpm * min(tm / 0.3, 1.0); ACM.TLoad = 0.0
        elif tm < 0.6:
            CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0.0
        elif tm < 0.8:
            CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = load_step
        elif tm < 1.0:
            CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0.0
        elif tm < 1.3:
            CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0.0
        else:
            CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0.0

        mt, wd = ACMSimPyIncremental(t0=t0, TIME=sim_TS, ACM=ACM, CTRL=CTRL,
                                      reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_spd)

        ds = int(CL_TS / (CL_TS / 1))  # downsampling = 1 since MACHINE_SIM=1
        for k in range(0, len(mt), ds):
            if ci >= n_ctrl:
                break

            t_tr[ci] = mt[k]
            theta_true_arr[ci] = wd[0][k]  # true elec angle
            omega_true_arr[ci] = wd[1][k]  # true speed rpm

            ia = wd[6][k]; ib = wd[7][k]
            ua = wd[28][k] + v_off[0]
            ub = wd[29][k] + v_off[1]

            u_mag_arr[ci] = np.sqrt(ua**2 + ub**2)
            i_mag_arr[ci] = np.sqrt(ia**2 + ib**2)
            # True EMF magnitude: KE * omega_elec
            omega_elec = omega_true_arr[ci] / 60 * 2 * np.pi * n_pp
            emf_true_arr[ci] = KE * abs(omega_elec)

            # ---- AF Observer (same logic as demo_sensorless_active_flux.py) ----
            af_a = af_psi_s[0] - L * ia
            af_b = af_psi_s[1] - L * ib
            af_amp = np.sqrt(af_a**2 + af_b**2)

            psi_err = KE - af_amp  # THIS is the Output Error (OE)
            af_corr_int[0] += af_Ki * psi_err * CL_TS
            af_corr_int[1] += af_Ki * psi_err * CL_TS

            if af_amp > 1e-10:
                u_a_n = af_a / af_amp; u_b_n = af_b / af_amp
            else:
                u_a_n = 1.0; u_b_n = 0.0

            corr_a = (af_Kp * psi_err + af_corr_int[0]) * u_a_n
            corr_b = (af_Kp * psi_err + af_corr_int[1]) * u_b_n

            af_psi_s[0] += CL_TS * (ua - R_obs * ia + corr_a)
            af_psi_s[1] += CL_TS * (ub - R_obs * ib + corr_b)

            af_a2 = af_psi_s[0] - L * ia
            af_b2 = af_psi_s[1] - L * ib
            af_theta[ci] = np.arctan2(af_b2, af_a2)
            af_psi_amp[ci] = np.sqrt(af_a2**2 + af_b2**2)
            af_oe[ci] = psi_err
            af_corr_mag[ci] = np.sqrt(corr_a**2 + corr_b**2)
            af_corr_int_mag[ci] = np.sqrt(af_corr_int[0]**2 + af_corr_int[1]**2)

            # ---- PLL Observer ----
            pll_th, pll_om, _ = pll_obs.step(ia, ib, ua, ub)
            pll_theta[ci] = pll_th
            pll_omega[ci] = pll_om

            # ---- SMO Observer ----
            smo_th, smo_om, _ = smo_obs.step(ia, ib, ua, ub)
            smo_theta[ci] = smo_th
            smo_omega[ci] = smo_om

            ci += 1

    n = ci
    elapsed = _time.time() - t0_sim
    if verbose:
        print(f'  Open-loop diag done in {elapsed:.1f}s, {n} steps')

    return {
        't': t_tr[:n], 'theta_true': theta_true_arr[:n], 'omega_true': omega_true_arr[:n],
        'af_theta': af_theta[:n], 'af_psi_amp': af_psi_amp[:n], 'af_oe': af_oe[:n],
        'af_corr_mag': af_corr_mag[:n], 'af_corr_int_mag': af_corr_int_mag[:n],
        'pll_theta': pll_theta[:n], 'pll_omega': pll_omega[:n],
        'smo_theta': smo_theta[:n], 'smo_omega': smo_omega[:n],
        'u_mag': u_mag_arr[:n], 'i_mag': i_mag_arr[:n], 'emf_true': emf_true_arr[:n],
        'n_pp': n_pp, 'KE': KE, 'L': L, 'R': R_true,
    }


def plot_diag(res, motor_name, suffix=''):
    """Plot comprehensive open-loop diagnostics."""
    t = res['t']; n_pp = res['n_pp']; KE = res['KE']
    ss = t > 0.4

    # Angle errors
    ae_af = np.degrees(angle_error(res['theta_true'], res['af_theta']))
    ae_pll = np.degrees(angle_error(res['theta_true'], res['pll_theta']))
    ae_smo = np.degrees(angle_error(res['theta_true'], res['smo_theta']))

    fig, axes = plt.subplots(5, 1, figsize=(14, 16), sharex=True)
    fig.suptitle(f'{motor_name} Open-Loop Observer Diagnostics{suffix}', fontsize=14, fontweight='bold')

    # 1. Speed
    axes[0].plot(t, res['omega_true'], 'b', lw=1)
    axes[0].set_ylabel('True Speed [rpm]')
    axes[0].set_title('Speed Profile (encoder-controlled)')

    # 2. Angle estimation errors
    axes[1].plot(t, ae_af, 'r', lw=0.8, label=f'AF (RMS_SS={np.sqrt(np.mean(ae_af[ss]**2)):.2f} deg)')
    axes[1].plot(t, ae_pll, 'g', lw=0.8, label=f'PLL (RMS_SS={np.sqrt(np.mean(ae_pll[ss]**2)):.2f} deg)')
    axes[1].plot(t, ae_smo, 'purple', lw=0.8, alpha=0.7, label=f'SMO (RMS_SS={np.sqrt(np.mean(ae_smo[ss]**2)):.2f} deg)')
    axes[1].set_ylabel('Angle Error [deg]')
    axes[1].legend(fontsize=9)
    axes[1].set_title('Angle Estimation Error (open-loop, no feedback)')

    # 3. AF flux amplitude and OE
    axes[2].plot(t, res['af_psi_amp'] * 1e3, 'b', lw=0.8, label='|psi_AF| [mWb]')
    axes[2].axhline(KE * 1e3, color='r', ls='--', lw=1, label=f'KE = {KE*1e3:.2f} mWb')
    axes[2].set_ylabel('Flux [mWb]')
    ax2r = axes[2].twinx()
    ax2r.plot(t, res['af_oe'] * 1e3, 'orange', lw=0.6, alpha=0.7, label='OE = KE - |psi_AF|')
    ax2r.set_ylabel('OE [mWb]', color='orange')
    axes[2].legend(loc='upper left', fontsize=9)
    ax2r.legend(loc='upper right', fontsize=9)
    axes[2].set_title('AF Flux Amplitude & Output Error (OE)')

    # 4. AF Correction magnitude (proportional + integral)
    axes[3].plot(t, res['af_corr_mag'], 'r', lw=0.8, label='|correction| (Kp*OE + int)')
    axes[3].plot(t, res['af_corr_int_mag'], 'b', lw=0.8, alpha=0.7, label='|integrator|')
    axes[3].set_ylabel('Correction [V]')
    axes[3].legend(fontsize=9)
    axes[3].set_title('AF PI Correction Magnitude vs Observer Gain')

    # 5. Signal ratios
    axes[4].plot(t, res['emf_true'], 'g', lw=0.8, label='|EMF| = KE*omega_e')
    axes[4].plot(t, res['u_mag'], 'b', lw=0.8, alpha=0.6, label='|u_ab|')
    axes[4].plot(t, res['i_mag'] * res['R'], 'r', lw=0.8, alpha=0.6, label='|i|*R (resistive drop)')
    axes[4].plot(t, res['af_corr_mag'], 'orange', lw=0.8, alpha=0.7, label='|AF correction|')
    axes[4].set_ylabel('[V]')
    axes[4].legend(fontsize=9)
    axes[4].set_xlabel('Time [s]')
    axes[4].set_title('Signal Magnitudes: EMF vs Voltage vs R*I vs Correction')

    plt.tight_layout()
    fname = f'fig_diag_{motor_name}{suffix}.png'
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f'  Saved {fname}')
    return fname


# ======================= MAIN =======================
if __name__ == '__main__':
    for motor_name in ['servo', 'small_L', 'big_L']:
        p = get_motor_preset(motor_name)
        cmd = p['cmd_rpm']
        KE = p['d']['init_KE']; L = p['d']['init_Lq']; R = p['d']['init_R']
        print(f'\n{"="*70}')
        print(f'  {motor_name}: KE={KE}, L={L*1e3:.2f}mH, R={R}, KE/L={KE/L:.1f}')
        print(f'{"="*70}')

        # 1. Ideal parameters (no mismatch)
        print('\n  --- Ideal (R_mis=1.0, v_off=0) ---')
        res_ideal = run_open_loop_diag(
            p['d'], CLBW_Hz=p['clbw'], zeta=p['zeta'],
            R_mismatch_factor=1.0, v_off=(0, 0),
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            pll_bw=100.0, cmd_rpm=cmd, load_step=p['load_step'])
        plot_diag(res_ideal, motor_name, '_ideal')

        # 2. With R mismatch + voltage offset
        print(f'\n  --- Mismatch (R_mis={p["r_mis"]}, v_off={p["v_off"]}) ---')
        res_mis = run_open_loop_diag(
            p['d'], CLBW_Hz=p['clbw'], zeta=p['zeta'],
            R_mismatch_factor=p['r_mis'], v_off=p['v_off'],
            af_Kp=p['af_kp'], af_Ki=p['af_ki'],
            pll_bw=100.0, cmd_rpm=cmd, load_step=p['load_step'])
        plot_diag(res_mis, motor_name, '_mismatch')

        # 3. Sweep AF gains and report OE vs angle error
        print(f'\n  --- AF Gain Sweep (ideal params) ---')
        gains = [
            (0,    0,     'no_corr'),
            (50,   500,   'low'),
            (100,  1000,  'med_low'),
            (500,  5000,  'default'),
            (1000, 10000, 'med_hi'),
            (2000, 50000, 'high'),
            (5000, 100000,'very_hi'),
        ]
        header = f'{"Kp":>6s} {"Ki":>7s} | {"AngRMS":>7s} {"OE_RMS":>8s} {"OE_Max":>8s} {"Corr_RMS":>9s} {"IntMag":>8s}'
        print(f'    {header}')
        print(f'    {"-"*65}')
        for kp, ki, label in gains:
            r = run_open_loop_diag(
                p['d'], CLBW_Hz=p['clbw'], zeta=p['zeta'],
                R_mismatch_factor=1.0, v_off=(0, 0),
                af_Kp=kp, af_Ki=ki, pll_bw=100.0,
                cmd_rpm=cmd, load_step=p['load_step'], verbose=False)
            ss = r['t'] > 0.4
            ae = np.degrees(angle_error(r['theta_true'], r['af_theta']))
            ang_rms = np.sqrt(np.mean(ae[ss]**2))
            oe_rms = np.sqrt(np.mean(r['af_oe'][ss]**2)) * 1e3
            oe_max = np.max(np.abs(r['af_oe'][ss])) * 1e3
            corr_rms = np.sqrt(np.mean(r['af_corr_mag'][ss]**2))
            int_mag = np.mean(r['af_corr_int_mag'][ss])
            print(f'    {kp:>6d} {ki:>7d} | {ang_rms:>7.2f} {oe_rms:>8.4f} {oe_max:>8.4f} {corr_rms:>9.4f} {int_mag:>8.3f}')
