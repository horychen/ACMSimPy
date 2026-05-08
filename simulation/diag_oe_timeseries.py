#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
=======================================================================
  OE (Output Error) 时域波形诊断
=======================================================================
核心问题：AF/PLL/SMO 各观测器在开环（编码器控制）条件下，
          OE = KE - |ψ_AF| 能否在零附近稳定波动？

实验设计：
  - Stage 0 编码器控制（速度+角度均用编码器真值）
  - 三台电机：servo / small_L / big_L
  - 两种参数：理想（无失配）/ 失配（R×1.5, 电压偏移）
  - 三种观测器：AF / PLL / SMO 并行运行
  - 输出：OE 时域波形 + 角度误差 + 统计表

对于 PLL 和 SMO：
  - PLL 的 "OE" 定义为：KE - |ψ_PLL_AF| （PLL 内部也用 AF 积分器）
  - SMO 的 "OE" 定义为：角度误差的正弦分量 * KE（反映磁链方向误差）
"""
import matplotlib
matplotlib.use('Agg')
from pylab import np, plt
import copy, os

from tutorials_ep6_svpwm import (The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental)
from tuner import get_coeffs_dc_motor_current_regulator, get_coeffs_dc_motor_SPEED_regulator
from demo_sensorless_active_flux import wrap_angle, angle_error, angle_diff_scalar
from eval_staged_tuning import get_motor_preset
from observers_alt import PLLFluxObserver, SlidingModeObserver


def run_oe_diag(d, CLBW_Hz, zeta, R_mis, v_off, af_Kp, af_Ki,
                pll_bw, cmd_rpm, load_step):
    """Run open-loop diag returning detailed OE time series."""
    dd = copy.deepcopy(d)
    R_true = dd['init_R']; R_obs = R_true * R_mis
    L = dd['init_Lq']; n_pp = dd['init_npp']; KE = dd['init_KE']
    CL_TS = dd['CL_TS']; J_s = dd['init_Js']
    VL_TS = CL_TS * dd['VL_EXE_PER_CL_EXE']

    cKp, cKi = get_coeffs_dc_motor_current_regulator(R_true, L, CLBW_Hz)
    sKp, sKi = get_coeffs_dc_motor_SPEED_regulator(J_s, n_pp, KE, zeta, cKp/L)
    dd['CL_SERIES_KP'] = cKp; dd['CL_SERIES_KI'] = cKi
    dd['VL_SERIES_KP'] = sKp; dd['VL_SERIES_KI'] = sKi

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
    Vm = dd['DC_BUS_VOLTAGE']/1.732; Im = vl_lim*1.414*dd['init_IN']
    reg_id = The_PID_Regulator(cKp, cKp*cKi*ki_f, 0,0, Vm, Vm, CL_TS)
    reg_iq = The_PID_Regulator(cKp, cKp*cKi*ki_f, 0,0, Vm, Vm, CL_TS)
    reg_spd = The_PID_Regulator(sKp, sKp*sKi, 0,0, Im, Im, VL_TS)

    # Observers
    af_psi_s = np.array([KE, 0.0]); af_corr_int = np.zeros(2)
    pll_obs = PLLFluxObserver(R_obs, L, KE, CL_TS, pll_bw=pll_bw, af_Kp=af_Kp, af_Ki=af_Ki)
    smo_obs = SlidingModeObserver(R_obs, L, KE, CL_TS, smo_gain=None, lpf_fc=200, pll_bw=pll_bw)

    total_time = dd['TIME_SLICE'] * dd['NUMBER_OF_SLICES']
    N = int(total_time / CL_TS) + 100
    t = np.zeros(N); omega = np.zeros(N)
    # AF
    af_oe_arr = np.zeros(N); af_ae_arr = np.zeros(N); af_amp_arr = np.zeros(N)
    af_corr_arr = np.zeros(N); af_int_arr = np.zeros(N)
    # PLL
    pll_oe_arr = np.zeros(N); pll_ae_arr = np.zeros(N)
    # SMO
    smo_oe_arr = np.zeros(N); smo_ae_arr = np.zeros(N)

    ci = 0
    for sl in range(dd['NUMBER_OF_SLICES']):
        t0 = sl * dd['TIME_SLICE']; tm = t0 + dd['TIME_SLICE']*0.5
        if   tm<0.3: CTRL.cmd_rpm = cmd_rpm*min(tm/0.3,1); ACM.TLoad = 0
        elif tm<0.6: CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0
        elif tm<0.8: CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = load_step
        elif tm<1.0: CTRL.cmd_rpm = cmd_rpm; ACM.TLoad = 0
        elif tm<1.3: CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0
        else:        CTRL.cmd_rpm = -cmd_rpm; ACM.TLoad = 0

        mt, wd = ACMSimPyIncremental(t0=t0, TIME=dd['TIME_SLICE'], ACM=ACM, CTRL=CTRL,
                                      reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_spd)
        for k in range(len(mt)):
            if ci >= N: break
            t[ci] = mt[k]; omega[ci] = wd[1][k]
            th_true = wd[0][k]
            ia = wd[6][k]; ib = wd[7][k]
            ua = wd[28][k]+v_off[0]; ub = wd[29][k]+v_off[1]

            # === AF ===
            af_a = af_psi_s[0]-L*ia; af_b = af_psi_s[1]-L*ib
            af_amp = np.sqrt(af_a**2+af_b**2)
            psi_err = KE-af_amp
            af_corr_int[0] += af_Ki*psi_err*CL_TS
            af_corr_int[1] += af_Ki*psi_err*CL_TS
            if af_amp>1e-10: un_a,un_b = af_a/af_amp, af_b/af_amp
            else: un_a,un_b = 1.0, 0.0
            ca = (af_Kp*psi_err+af_corr_int[0])*un_a
            cb = (af_Kp*psi_err+af_corr_int[1])*un_b
            af_psi_s[0] += CL_TS*(ua-R_obs*ia+ca)
            af_psi_s[1] += CL_TS*(ub-R_obs*ib+cb)
            af_a2 = af_psi_s[0]-L*ia; af_b2 = af_psi_s[1]-L*ib
            th_af = np.arctan2(af_b2, af_a2)
            af_oe_arr[ci] = psi_err  # OE = KE - |ψ_AF|
            af_ae_arr[ci] = angle_diff_scalar(th_true, th_af)
            af_amp_arr[ci] = af_amp
            af_corr_arr[ci] = np.sqrt(ca**2+cb**2)
            af_int_arr[ci] = np.sqrt(af_corr_int[0]**2+af_corr_int[1]**2)

            # === PLL ===
            pll_th, pll_om, pll_psi = pll_obs.step(ia, ib, ua, ub)
            pll_ae_arr[ci] = angle_diff_scalar(th_true, pll_th)
            # PLL 的 OE：PLL 内部也维护了 AF 积分器
            pll_amp = np.sqrt(pll_psi[0]**2 + pll_psi[1]**2) if pll_psi is not None else KE
            pll_oe_arr[ci] = KE - pll_amp

            # === SMO ===
            smo_th, smo_om, smo_psi = smo_obs.step(ia, ib, ua, ub)
            smo_ae_arr[ci] = angle_diff_scalar(th_true, smo_th)
            smo_amp = np.sqrt(smo_psi[0]**2 + smo_psi[1]**2) if smo_psi is not None else KE
            smo_oe_arr[ci] = KE - smo_amp

            ci += 1

    n = ci
    return {
        't': t[:n], 'omega': omega[:n], 'KE': KE, 'L': L, 'n_pp': n_pp,
        'af_oe': af_oe_arr[:n], 'af_ae': af_ae_arr[:n], 'af_amp': af_amp_arr[:n],
        'af_corr': af_corr_arr[:n], 'af_int': af_int_arr[:n],
        'pll_oe': pll_oe_arr[:n], 'pll_ae': pll_ae_arr[:n],
        'smo_oe': smo_oe_arr[:n], 'smo_ae': smo_ae_arr[:n],
    }


def plot_oe_comparison(res, motor_name, suffix, out_dir):
    """时域 OE 对比图 + 角度误差对比图"""
    t = res['t']; KE = res['KE']; ss = t > 0.4

    fig, axes = plt.subplots(4, 1, figsize=(16, 14), sharex=True)
    fig.suptitle(f'{motor_name} — 开环 OE 时域波形 ({suffix})', fontsize=14, fontweight='bold')

    # 1. Speed profile
    axes[0].plot(t, res['omega'], 'b', lw=0.8)
    axes[0].set_ylabel('Speed [rpm]')
    axes[0].set_title('编码器控制速度剖面')
    axes[0].grid(True, alpha=0.3)

    # 2. OE time series (mWb)
    axes[1].plot(t, res['af_oe']*1e3, 'r', lw=0.6, label=f'AF OE (RMS={np.sqrt(np.mean(res["af_oe"][ss]**2))*1e3:.3f} mWb)')
    axes[1].plot(t, res['pll_oe']*1e3, 'g', lw=0.6, label=f'PLL OE (RMS={np.sqrt(np.mean(res["pll_oe"][ss]**2))*1e3:.3f} mWb)')
    axes[1].plot(t, res['smo_oe']*1e3, 'purple', lw=0.6, alpha=0.7, label=f'SMO OE (RMS={np.sqrt(np.mean(res["smo_oe"][ss]**2))*1e3:.3f} mWb)')
    axes[1].axhline(0, color='k', ls='-', lw=0.5)
    axes[1].set_ylabel('OE [mWb]')
    axes[1].legend(fontsize=9)
    axes[1].set_title(f'Output Error = KE − |ψ_AF| (KE={KE*1e3:.1f} mWb)')
    axes[1].grid(True, alpha=0.3)
    # 标注 ±5% KE 范围
    pct5 = KE * 0.05 * 1e3
    axes[1].axhline(pct5, color='gray', ls=':', lw=0.5, alpha=0.5)
    axes[1].axhline(-pct5, color='gray', ls=':', lw=0.5, alpha=0.5)
    axes[1].text(t[-1]*0.98, pct5, f'±5% KE', fontsize=7, ha='right', va='bottom', color='gray')

    # 3. Angle error time series (deg)
    axes[2].plot(t, np.degrees(res['af_ae']), 'r', lw=0.6, label=f'AF (RMS={np.degrees(np.sqrt(np.mean(res["af_ae"][ss]**2))):.2f}°)')
    axes[2].plot(t, np.degrees(res['pll_ae']), 'g', lw=0.6, label=f'PLL (RMS={np.degrees(np.sqrt(np.mean(res["pll_ae"][ss]**2))):.2f}°)')
    axes[2].plot(t, np.degrees(res['smo_ae']), 'purple', lw=0.6, alpha=0.7, label=f'SMO (RMS={np.degrees(np.sqrt(np.mean(res["smo_ae"][ss]**2))):.2f}°)')
    axes[2].axhline(0, color='k', ls='-', lw=0.5)
    axes[2].set_ylabel('Angle Error [deg]')
    axes[2].legend(fontsize=9)
    axes[2].set_title('角度估计误差')
    axes[2].grid(True, alpha=0.3)

    # 4. AF 修正量 vs EMF 信号比
    omega_e = res['omega']/60*2*np.pi*res['n_pp']
    emf = KE * np.abs(omega_e)
    axes[3].plot(t, res['af_corr'], 'r', lw=0.6, label='|AF correction| [V]')
    axes[3].plot(t, emf, 'g', lw=0.8, label='|EMF|=KE·ω_e [V]')
    ratio = np.zeros_like(emf)
    mask = emf > 0.1
    ratio[mask] = res['af_corr'][mask] / emf[mask] * 100
    ax3r = axes[3].twinx()
    ax3r.plot(t[mask], ratio[mask], 'b', lw=0.4, alpha=0.5, label='correction/EMF [%]')
    ax3r.set_ylabel('correction/EMF [%]', color='b')
    axes[3].set_ylabel('[V]')
    axes[3].legend(loc='upper left', fontsize=9)
    ax3r.legend(loc='upper right', fontsize=9)
    axes[3].set_title('AF修正量 vs EMF强度')
    axes[3].set_xlabel('Time [s]')
    axes[3].grid(True, alpha=0.3)

    plt.tight_layout()
    fname = os.path.join(out_dir, f'oe_timeseries_{motor_name}_{suffix}.png')
    plt.savefig(fname, dpi=150)
    plt.close()
    return fname


def compute_stats(res):
    """计算稳态OE统计量"""
    t = res['t']; ss = t > 0.4; KE = res['KE']
    stats = {}
    for name, oe, ae in [('AF', res['af_oe'], res['af_ae']),
                          ('PLL', res['pll_oe'], res['pll_ae']),
                          ('SMO', res['smo_oe'], res['smo_ae'])]:
        oe_ss = oe[ss]; ae_ss = ae[ss]
        stats[name] = {
            'oe_rms_mWb': np.sqrt(np.mean(oe_ss**2)) * 1e3,
            'oe_max_mWb': np.max(np.abs(oe_ss)) * 1e3,
            'oe_mean_mWb': np.mean(oe_ss) * 1e3,
            'oe_pct_KE': np.sqrt(np.mean(oe_ss**2)) / KE * 100,
            'ae_rms_deg': np.degrees(np.sqrt(np.mean(ae_ss**2))),
            'ae_max_deg': np.degrees(np.max(np.abs(ae_ss))),
        }
    return stats


# ===================== MAIN =====================
if __name__ == '__main__':
    out_dir = os.path.join(os.path.dirname(__file__), 'docs', 'sensorless_report_assets')
    os.makedirs(out_dir, exist_ok=True)
    
    # Also copy to repo root
    repo_dir = os.path.join(os.path.dirname(__file__), 'docs')
    os.makedirs(repo_dir, exist_ok=True)

    all_stats = {}
    
    for motor_name in ['servo', 'small_L', 'big_L']:
        p = get_motor_preset(motor_name)
        cmd = p['cmd_rpm']; KE = p['d']['init_KE']; L = p['d']['init_Lq']
        R = p['d']['init_R']
        print(f'\n{"="*70}')
        print(f'  {motor_name}: KE={KE*1e3:.1f}mWb, L={L*1e3:.2f}mH, R={R}Ω, KE/L={KE/L:.1f}')
        print(f'{"="*70}')

        for suffix, r_mis, voff in [('ideal', 1.0, (0,0)),
                                     ('mismatch', p['r_mis'], p['v_off'])]:
            print(f'\n  --- {suffix} (R_mis={r_mis}, v_off={voff}) ---')
            res = run_oe_diag(
                p['d'], CLBW_Hz=p['clbw'], zeta=p['zeta'],
                R_mis=r_mis, v_off=voff,
                af_Kp=p['af_kp'], af_Ki=p['af_ki'],
                pll_bw=100, cmd_rpm=cmd, load_step=p['load_step'])

            fname = plot_oe_comparison(res, motor_name, suffix, out_dir)
            print(f'  Saved {fname}')

            stats = compute_stats(res)
            all_stats[f'{motor_name}_{suffix}'] = stats

            # Print stats table
            print(f'  {"Observer":>6s} | {"OE_RMS":>8s} {"OE_Max":>8s} {"OE/%KE":>7s} | {"AngleRMS":>9s} {"AngleMax":>9s} | {"OE≈0?":>6s}')
            print(f'  {"-"*70}')
            for obs in ['AF', 'PLL', 'SMO']:
                s = stats[obs]
                ok = '✓' if s['oe_pct_KE'] < 5 else ('~' if s['oe_pct_KE'] < 20 else '✗')
                print(f'  {obs:>6s} | {s["oe_rms_mWb"]:>7.3f}m {s["oe_max_mWb"]:>7.3f}m {s["oe_pct_KE"]:>6.1f}% | {s["ae_rms_deg"]:>8.2f}° {s["ae_max_deg"]:>8.2f}° | {ok:>6s}')

    # Print summary table
    print(f'\n\n{"="*90}')
    print(f'  OE 稳态统计汇总 (OE = KE - |ψ_obs|)')
    print(f'{"="*90}')
    print(f'  {"Motor":>8s} {"Cond":>10s} | {"AF OE%":>7s} {"AF°":>6s} | {"PLL OE%":>8s} {"PLL°":>6s} | {"SMO OE%":>8s} {"SMO°":>6s} |')
    print(f'  {"-"*80}')
    for key in sorted(all_stats.keys()):
        parts = key.rsplit('_', 1)
        motor = parts[0]; cond = parts[1]
        s = all_stats[key]
        af = s['AF']; pll = s['PLL']; smo = s['SMO']
        af_ok = '✓' if af['oe_pct_KE']<5 else '✗'
        pll_ok = '✓' if pll['oe_pct_KE']<5 else '✗'
        smo_ok = '✓' if smo['oe_pct_KE']<5 else '✗'
        print(f'  {motor:>8s} {cond:>10s} | {af["oe_pct_KE"]:>5.1f}%{af_ok} {af["ae_rms_deg"]:>5.1f}° | {pll["oe_pct_KE"]:>6.1f}%{pll_ok} {pll["ae_rms_deg"]:>5.1f}° | {smo["oe_pct_KE"]:>6.1f}%{smo_ok} {smo["ae_rms_deg"]:>5.1f}° |')
