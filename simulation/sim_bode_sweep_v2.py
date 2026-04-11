# -*- coding: utf-8 -*-
"""
数值仿真波特图 v2 — 跟踪与抗扰频域分析
=========================================
相比 v1 的改进：
  1. 对比"不同控制器参数 (zeta, CLBW)" 与"不同观测器带宽 (omega_ob)" 的配置组合，
     而不仅仅是 ESO on/off，更贴合 sim_ramp_load_comparison.py 的多参数对比思路。
  2. 叠加解析传递函数 Bode 曲线作为参照基准（来自 python-control 库），
     让数值仿真结果和理论推导可以对比验证。
  3. DFT 提取增强：整周期截断 + 去直流 + 去重复时间戳，避免坏点。
  4. 所有参数均通过 CLI 传入，方便快速出图。

使用示例：
  python sim_bode_sweep_v2.py
  python sim_bode_sweep_v2.py --num-points 30 --freq-end 500
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import time
import control

from tuner import (
    get_coeffs_dc_motor_current_regulator,
    get_coeffs_dc_motor_SPEED_regulator,
)
from tutorials_ep6_svpwm import (
    The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental,
)


# ======================================================================
# 基础电机参数（与 sim_ramp_load_comparison.py 一致）
# ======================================================================
MOTOR_PARAMS = {
    'CL_TS': 1e-4,
    'VL_EXE_PER_CL_EXE': 5,
    'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
    'init_npp': 22,
    'init_IN': 1.3 * 6 / 1.414,
    'init_R': 0.035,
    'init_Ld': 1 * 0.036e-3,
    'init_Lq': 1 * 0.036e-3,
    'init_KE': 0.0125,
    'init_Rreq': 0.0,
    'init_Js': 0.44e-4,
    'DC_BUS_VOLTAGE': 48,
    'CTRL.bool_apply_speed_closed_loop_control': True,
    'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
    'CTRL.bool_apply_sweeping_frequency_excitation': False,
    'CTRL.bool_overwrite_speed_commands': True,
    'CTRL.bool_zero_id_control': True,
    'FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False': 10,
    'VL_LIMIT_OVERLOAD_FACTOR': 3.0,
}


# ======================================================================
# DFT 提取
# ======================================================================
def extract_mag_phase_dft(t_arr, y_arr, u_arr, freq_Hz):
    """
    精确单频 DFT 提取幅值和相位。
    做了整周期截断 + 去直流偏置 + 容错保护。
    """
    period = 1.0 / freq_Hz
    total_time = t_arr[-1] - t_arr[0]
    num_periods = max(1, int(np.floor(total_time / period)))

    valid_duration = num_periods * period
    valid_start_time = t_arr[-1] - valid_duration

    idx = t_arr >= valid_start_time
    t_v = t_arr[idx]
    y_v = y_arr[idx]
    u_v = u_arr[idx]

    # 去直流
    y_v = y_v - np.mean(y_v)
    u_v = u_v - np.mean(u_v)

    if len(t_v) < 2:
        return 0.0, 0.0

    dt = t_v[1] - t_v[0]
    if dt <= 0:
        return 0.0, 0.0

    basis = np.exp(-1j * 2 * np.pi * freq_Hz * t_v)
    Y_dft = np.sum(y_v * basis) * dt
    U_dft = np.sum(u_v * basis) * dt

    if np.abs(U_dft) < 1e-15:
        return 0.0, 0.0

    G = Y_dft / U_dft
    mag_db = 20 * np.log10(max(np.abs(G), 1e-30))
    phase_deg = np.degrees(np.angle(G))
    return mag_db, phase_deg


# ======================================================================
# 解析传递函数
# ======================================================================
def get_analytical_bode(d, zeta, CLBW_Hz, omega_ob, freqs_Hz, enable_ESO=False):
    """
    基于 python-control 求解析法波特图。
    返回 tracking 和 disturbance 通道的 (mag_dB, phase_deg) 数组。

    跟踪通道 T(s) = Gw_closed = Gw_open / (1 + Gw_open)
    抗扰通道 S_d(s) = (1/(J_s*s)) / (1 + Gw_open) * (RPM转换系数)
    若启用 ESO 前馈，则抗扰通道需额外乘以 (1 - L_ESO(s))

    注意：这里的解析模型是连续域理想化模型，不含离散化、饱和等非线性。
    """
    R    = d['init_R']
    L    = d['init_Lq']
    J_s  = d['init_Js']
    n_pp = d['init_npp']
    KE   = d['init_KE']
    KA   = KE

    currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R, L, CLBW_Hz)
    currentBandwidth_radPerSec = currentKp / L
    speedKp, speedKi = get_coeffs_dc_motor_SPEED_regulator(
        J_s, n_pp, KA, zeta, currentBandwidth_radPerSec
    )

    # 电流闭环（零极点对消后）
    Gi_closed = control.tf([1], [L / currentKp, 1])

    # 运动方程 (电流 → 电角速度)
    KT = 1.5 * n_pp * KA
    # 运动方程输出为 elec.rad/s，但我们的 speed 用 RPM = omega_mech * 30/pi
    # omega_mech = omega_elec / n_pp
    # 所以 Plant: iq → RPM = KT*n_pp/J_s / s * (30/pi/n_pp)
    # 合并后: KT * 30 / (pi * J_s) / s
    rpm_per_radps = 30.0 / np.pi   # 机械角速度 → RPM
    plant_iq_to_rpm = control.tf([KT / (J_s / n_pp) * rpm_per_radps / n_pp], [1, 0])

    # speedPI: 输入=误差(RPM), 输出=iq*
    speedPI = control.tf([speedKp, speedKp * speedKi], [1, 0])

    # 开环: 误差(RPM) → speedPI → iq* → Gi_closed → iq → plant → RPM
    Gw_open = plant_iq_to_rpm * Gi_closed * speedPI

    # 跟踪闭环
    T_s = Gw_open / (1 + Gw_open)

    # 抗扰通道: TLoad → RPM
    # d_omega_mech/dt = (KT*iq - TLoad) * n_pp / J_s  (elec domain)
    # 但对于闭环抗扰灵敏度，扰动从 TLoad 注入到 RPM:
    # S_d = -plant_TL_to_RPM / (1 + Gw_open)
    # plant_TL_to_RPM = n_pp/J_s / s * rpm_per_radps / n_pp = rpm_per_radps / J_s / s
    # (负号是因为 TLoad 对输出 RPM 是负贡献，但取绝对值/增益即可)
    plant_TL_to_rpm = control.tf([rpm_per_radps / (J_s / n_pp) / n_pp], [1, 0])
    S_d = plant_TL_to_rpm / (1 + Gw_open)

    # 如果使用 ESO 前馈补偿
    if enable_ESO:
        # 4阶 ESO 只在这里做简化：ESO 观测到的扰动的传递函数
        # 对于 (s+omega_ob)^4 的 ESO，其扰动估计的传递函数为：
        # L_d(s) = (4*w*s^3 + 6*w^2*s^2 + 4*w^3*s + w^4) / (s + w)^4
        #        = 1 - s^4 / (s + w)^4
        # ESO 补偿后的抗扰 = S_d * (1 - L_d(s)) = S_d * s^4 / (s+w)^4
        w = omega_ob
        eso_residual = control.tf([1, 0, 0, 0, 0],                 # s^4
                                   np.polymul(np.polymul([1, w], [1, w]),
                                              np.polymul([1, w], [1, w])))  # (s+w)^4
        S_d = S_d * eso_residual

    omega_array = 2 * np.pi * freqs_Hz

    # Tracking
    mag_t, phase_t_rad, _ = control.bode_plot(T_s, omega_array, dB=False, Hz=False, deg=True, plot=False)
    mag_t_db = 20 * np.log10(np.maximum(mag_t, 1e-30))
    phase_t_deg = np.degrees(phase_t_rad)

    # Disturbance
    mag_d, phase_d_rad, _ = control.bode_plot(S_d, omega_array, dB=False, Hz=False, deg=True, plot=False)
    mag_d_db = 20 * np.log10(np.maximum(mag_d, 1e-30))
    phase_d_deg = np.degrees(phase_d_rad)

    return {
        'tracking':    {'mag': mag_t_db, 'phase': phase_t_deg},
        'disturbance': {'mag': mag_d_db, 'phase': phase_d_deg},
    }


# ======================================================================
# 单频仿真
# ======================================================================
def run_single_frequency(d, freq_Hz, zeta, CLBW_Hz, enable_ESO, omega_ob, mode='tracking', rpm_0=500.0):
    """
    单个频率点的数值仿真注入。
    mode='tracking' → 正弦注入速度给定, 观测速度响应
    mode='disturbance' → 正弦注入负载转矩, 观测速度响应
    """
    dd = copy.deepcopy(d)

    R    = dd['init_R']
    L    = dd['init_Lq']
    J_s  = dd['init_Js']
    n_pp = dd['init_npp']
    KE   = dd['init_KE']
    KA   = KE
    CL_TS = dd['CL_TS']
    VL_TS = dd['CL_TS'] * dd['VL_EXE_PER_CL_EXE']

    currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R, L, CLBW_Hz)
    currentBandwidth_radPerSec = currentKp / L
    speedKp, speedKi = get_coeffs_dc_motor_SPEED_regulator(
        J_s, n_pp, KA, zeta, currentBandwidth_radPerSec
    )

    CTRL = The_Motor_Controller(
        CL_TS=CL_TS, VL_TS=VL_TS,
        init_npp=n_pp, init_IN=dd['init_IN'],
        init_R=R, init_Ld=dd['init_Ld'], init_Lq=L,
        init_KE=KE, init_Rreq=dd['init_Rreq'], init_Js=J_s,
        DC_BUS_VOLTAGE=dd['DC_BUS_VOLTAGE'],
    )
    CTRL.bool_apply_decoupling_voltages_to_current_regulation = dd['CTRL.bool_apply_decoupling_voltages_to_current_regulation']
    CTRL.bool_apply_sweeping_frequency_excitation = False
    CTRL.bool_overwrite_speed_commands = True
    CTRL.bool_zero_id_control = dd['CTRL.bool_zero_id_control']
    CTRL.bool_apply_speed_closed_loop_control = True

    if enable_ESO:
        CTRL.index_separate_speed_estimation = 1
        CTRL.use_disturbance_feedforward_rejection = 1
        CTRL.ell1 = 4 * omega_ob
        CTRL.ell2 = 6 * omega_ob**2
        CTRL.ell3 = 4 * omega_ob**3 * J_s / n_pp
        CTRL.ell4 =     omega_ob**4 * J_s / n_pp

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=dd['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])

    local_Ki_factor = 1.0 if dd.get('CTRL.bool_apply_decoupling_voltages_to_current_regulation', False) \
        else dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False', 10)

    reg_id = The_PID_Regulator(currentKp, currentKp * currentKi * local_Ki_factor, 0.0, 0.0,
                                dd['DC_BUS_VOLTAGE'] / 1.732, dd['DC_BUS_VOLTAGE'] / 1.732, CL_TS)
    reg_iq = The_PID_Regulator(currentKp, currentKp * currentKi * local_Ki_factor, 0.0, 0.0,
                                dd['DC_BUS_VOLTAGE'] / 1.732, dd['DC_BUS_VOLTAGE'] / 1.732, CL_TS)
    limit_factor = 5.0
    reg_speed = The_PID_Regulator(speedKp, speedKp * speedKi, 0.0, 0.0,
                                   limit_factor * dd['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * dd['init_IN'],
                                   limit_factor * dd['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * dd['init_IN'],
                                   VL_TS)

    # 仿真时间
    slice_dt = VL_TS
    periods_transient = 5.0
    periods_record    = 5.0
    # 下限必须足够让电机从零加速到 rpm_0 并让 ESO 完全收敛
    # ESO 的 4 阶滤波器需要 ~5/(omega_ob) 秒收敛 (omega_ob=200 → 0.025s)
    # 但 PI 速度环的稳态建立需要更久，设 1.0s 保底
    transient_time = max(1.0, periods_transient / freq_Hz)
    record_time    = max(0.1, periods_record / freq_Hz)
    total_time = transient_time + record_time
    n_slices = int(total_time / slice_dt)
    record_start_slice = int(transient_time / slice_dt)

    t_arr, y_arr, u_arr = [], [], []
    T_0   = 0.0
    rpm_amp = 10.0
    T_amp   = 1.0

    for si in range(n_slices):
        t_now = si * slice_dt
        is_recording = (si >= record_start_slice)

        if not is_recording:
            # ===== 暖机阶段：恒定指令，零扰动，让电机 + ESO 稳定到工作点 =====
            CTRL.cmd_rpm = rpm_0
            ACM.TLoad = T_0
        else:
            # ===== 录数据阶段：注入正弦激励 =====
            if mode == 'tracking':
                u_val = rpm_amp * np.sin(2 * np.pi * freq_Hz * t_now)
                CTRL.cmd_rpm = rpm_0 + u_val
                ACM.TLoad = T_0
            else:
                u_val = T_amp * np.sin(2 * np.pi * freq_Hz * t_now)
                CTRL.cmd_rpm = rpm_0
                ACM.TLoad = T_0 + u_val

        machine_times, watch_data = ACMSimPyIncremental(
            t0=t_now, TIME=slice_dt,
            ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed
        )

        if is_recording:
            t_arr.extend(machine_times)
            y_arr.extend(watch_data[1])  # speed_rpm
            if mode == 'tracking':
                u_arr.extend(rpm_amp * np.sin(2 * np.pi * freq_Hz * np.array(machine_times)))
            else:
                u_arr.extend(T_amp * np.sin(2 * np.pi * freq_Hz * np.array(machine_times)))

    t_arr = np.array(t_arr)
    y_arr = np.array(y_arr)
    u_arr = np.array(u_arr)

    # 去重复时间戳
    t_arr, unique_idx = np.unique(t_arr, return_index=True)
    y_arr = y_arr[unique_idx]
    u_arr = u_arr[unique_idx]

    y_ac = y_arr - rpm_0
    return extract_mag_phase_dft(t_arr, y_ac, u_arr, freq_Hz)


# ======================================================================
# 完整扫频
# ======================================================================
def sweep_bode(d, freqs, zeta, CLBW_Hz, enable_ESO, omega_ob, label="", rpm_0=500.0):
    """
    对一组完整频率做 tracking + disturbance 扫频。
    """
    result = {
        'tracking':    {'mag': [], 'phase': []},
        'disturbance': {'mag': [], 'phase': []},
    }

    for mode in ['tracking', 'disturbance']:
        for i, f in enumerate(freqs):
            m, p = run_single_frequency(d, f, zeta, CLBW_Hz, enable_ESO, omega_ob, mode=mode, rpm_0=rpm_0)
            result[mode]['mag'].append(m)
            # 相位连续性
            if len(result[mode]['phase']) > 0:
                prev_p = result[mode]['phase'][-1]
                while p - prev_p > 180: p -= 360
                while p - prev_p < -180: p += 360
            result[mode]['phase'].append(p)

    for mode in ['tracking', 'disturbance']:
        result[mode]['mag']   = np.array(result[mode]['mag'])
        result[mode]['phase'] = np.array(result[mode]['phase'])

    return result


# ======================================================================
# MAIN
# ======================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Numerical Bode Sweep v2 — 跟踪/抗扰频域分析，"
                    "对比不同控制器与观测器参数配置",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--num-points', type=int, default=20,
                        help="频率扫描点数 (default: 20)")
    parser.add_argument('--freq-start', type=float, default=1.0,
                        help="起始频率 Hz (default: 1.0)")
    parser.add_argument('--freq-end', type=float, default=200.0,
                        help="终止频率 Hz (default: 200.0)")
    parser.add_argument('--rpm', type=float, default=500.0,
                        help="工作点转速 RPM (default: 500.0)")
    args = parser.parse_args()

    d = copy.deepcopy(MOTOR_PARAMS)
    freqs = np.logspace(np.log10(args.freq_start), np.log10(args.freq_end), args.num_points)
    # 解析用的更密集的频率数组
    freqs_dense = np.logspace(np.log10(args.freq_start), np.log10(args.freq_end), 500)

    # ==================================================================
    # 定义对比配置
    # ==================================================================
    # 我们对比三种典型配置（和 sim_ramp_load_comparison.py 里的三组一致）：
    #   配置 A: zeta=5,  CLBW=500Hz,  无 ESO    — 低带宽/保守
    #   配置 B: zeta=15, CLBW=1000Hz, 无 ESO    — 中带宽/工程典型
    #   配置 C: zeta=15, CLBW=1000Hz, ESO ω=200 — 中带宽 + ESO 前馈补偿
    #
    # 这样既能对比「不同控制器参数」(A vs B)，
    # 又能对比「有无观测器」(B vs C)，图上三条曲线一目了然。

    configs = [
        {
            'label':  r'$\zeta$=5, CLBW=500Hz, No ESO',
            'zeta': 5, 'CLBW_Hz': 500, 'enable_ESO': False, 'omega_ob': 0,
            'color': '#e74c3c', 'ls': '-', 'marker': 's',
        },
        {
            'label':  r'$\zeta$=15, CLBW=1000Hz, No ESO',
            'zeta': 15, 'CLBW_Hz': 1000, 'enable_ESO': False, 'omega_ob': 0,
            'color': '#2ecc71', 'ls': '-', 'marker': '^',
        },
        {
            'label':  r'$\zeta$=15, CLBW=1000Hz, ESO $\omega_{ob}$=200',
            'zeta': 15, 'CLBW_Hz': 1000, 'enable_ESO': True, 'omega_ob': 200,
            'color': '#3498db', 'ls': '-', 'marker': 'o',
        },
    ]

    # ==================================================================
    # 1. 数值仿真扫频
    # ==================================================================
    print('=' * 65)
    print(' 数值仿真波特图 v2 — 跟踪/抗扰频域分析')
    print('=' * 65)
    print(f'  频率范围: {args.freq_start} ~ {args.freq_end} Hz, {args.num_points} 点')
    print(f'  工作点转速: {args.rpm} RPM')
    print(f'  DC_BUS_VOLTAGE: {d["DC_BUS_VOLTAGE"]} V, VL_LIMIT_OVERLOAD_FACTOR: {d["VL_LIMIT_OVERLOAD_FACTOR"]}')
    print(f'  对比配置: {len(configs)} 组')
    print()

    sim_results = {}
    t_start = time.time()

    for ci, cfg in enumerate(configs):
        tag = f"[{ci+1}/{len(configs)}] {cfg['label']}"
        print(f'  {tag}')
        sim_results[ci] = sweep_bode(
            d, freqs,
            cfg['zeta'], cfg['CLBW_Hz'], cfg['enable_ESO'], cfg['omega_ob'],
            label=cfg['label'], rpm_0=args.rpm,
        )
        elapsed_so_far = time.time() - t_start
        print(f'    完成 ({elapsed_so_far:.1f}s elapsed)')

    total_elapsed = time.time() - t_start
    print(f'\n  总耗时: {total_elapsed:.1f}s')

    # ==================================================================
    # 2. 解析传递函数 Bode（密集频率，画虚线作为理论参照）
    # ==================================================================
    print('\n  计算解析传递函数 Bode...')
    ana_results = {}
    for ci, cfg in enumerate(configs):
        ana_results[ci] = get_analytical_bode(
            d, cfg['zeta'], cfg['CLBW_Hz'], cfg['omega_ob'], freqs_dense,
            enable_ESO=cfg['enable_ESO'],
        )

    # ==================================================================
    # 3. 绘图：2×2 布局 (Tracking Mag/Phase, Disturbance Mag/Phase)
    # ==================================================================
    print('  绘图...')
    plt.style.use('bmh')
    mpl.rc('font', family='Times New Roman', size=10.0)
    mpl.rc('legend', fontsize=7.5)
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.unicode_minus'] = False

    fig, axes = plt.subplots(2, 2, figsize=(15, 10), dpi=150)

    titles = {
        (0,0): 'Tracking: Magnitude  (speed / cmd)',
        (1,0): 'Tracking: Phase',
        (0,1): 'Disturbance: Magnitude  (speed / TLoad)',
        (1,1): 'Disturbance: Phase',
    }
    ylabels = {
        (0,0): 'Magnitude [dB]', (1,0): 'Phase [deg]',
        (0,1): 'Magnitude [dB]', (1,1): 'Phase [deg]',
    }
    xlabels = {(1,0): 'Frequency [Hz]', (1,1): 'Frequency [Hz]'}
    chan_map = {0: 'tracking', 1: 'disturbance'}
    row_map  = {0: 'mag',     1: 'phase'}

    for ci, cfg in enumerate(configs):
        color   = cfg['color']
        marker  = cfg['marker']
        label   = cfg['label']

        for col, chan in chan_map.items():
            for row, field in row_map.items():
                ax = axes[row, col]

                ana_data = ana_results[ci][chan][field]
                sim_data = sim_results[ci][chan][field]

                # 对解析相位做匹配仿真的 unwrap
                if field == 'phase':
                    # 用仿真在第一个频率点的相位来对齐解析相位的起始偏移
                    sim_start = sim_data[0] if len(sim_data) > 0 else 0
                    # 找到 ana_data 中第一个点与 sim_start 最近的 360 度倍数偏移
                    ana_start = ana_data[0]
                    offset = round((sim_start - ana_start) / 360.0) * 360.0
                    ana_data = ana_data + offset

                # 解析（虚线），只在 Mag 行的第一个 config 标一次 "Analytical" 图例
                ana_label = None
                if row == 0 and ci == 0:
                    ana_label = 'Analytical (cont. TF)'
                ax.semilogx(
                    freqs_dense, ana_data,
                    color=color, ls='--', lw=1.0, alpha=0.5,
                    label=ana_label,
                )

                # 数值仿真（实线 + markers）
                sim_label = None
                if row == 0:
                    sim_label = label
                ax.semilogx(
                    freqs, sim_data,
                    color=color, ls=cfg['ls'], marker=marker, markersize=4,
                    label=sim_label,
                )

    # 格式化
    for pos, title in titles.items():
        axes[pos].set_title(title, fontweight='bold', fontsize=11)
    for pos, ylabel in ylabels.items():
        axes[pos].set_ylabel(ylabel)
    for pos, xlabel in xlabels.items():
        axes[pos].set_xlabel(xlabel)
    for ax in axes.flat:
        ax.grid(True, which='both', ls='--', alpha=0.4)

    # 图例只在幅值行放
    axes[0,0].legend(loc='lower left', fontsize=7, framealpha=0.9)
    axes[0,1].legend(loc='upper right', fontsize=7, framealpha=0.9)

    # -3dB 参考线（跟踪通道）
    axes[0,0].axhline(-3, color='gray', ls=':', lw=0.8, alpha=0.6, label='-3 dB')

    # 标注说明：虚线 = 解析，实线 = 仿真
    fig.text(0.5, 0.01, 'Dashed = Analytical (continuous-time TF)     Solid + Markers = Numerical Simulation',
             ha='center', fontsize=9, style='italic', color='gray')

    # 参数标注框
    param_text = (
        f"$V_{{dc}}$ = {d['DC_BUS_VOLTAGE']} V,  "
        f"VL_LIMIT_OVERLOAD = {d['VL_LIMIT_OVERLOAD_FACTOR']},  "
        f"$\\omega_{{0}}$ = {args.rpm} RPM"
    )
    fig.text(0.5, 0.96, param_text,
             ha='center', fontsize=9, color='#555',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0', edgecolor='#ccc', alpha=0.8))

    fig.suptitle(
        'Speed Loop Frequency Response: Tracking & Disturbance Rejection\n'
        '(Numerical Sine Sweep vs Analytical Transfer Function)',
        fontsize=13, fontweight='bold', y=1.01
    )

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    save_path = f'fig_sim_bode_sweep_v2_rpm{int(args.rpm)}.png'
    fig.savefig(save_path, dpi=200, bbox_inches='tight')
    print(f'\n  图已保存: {save_path}')
    plt.close('all')
    print('  完成!')


if __name__ == "__main__":
    main()
