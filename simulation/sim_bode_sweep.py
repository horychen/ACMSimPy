# -*- coding: utf-8 -*-
"""
数值仿真扫出 Bode 图：对比传统 PI 与 4阶 ESO
============================================
通过对速度闭环和观测器系统直接注入不同频率的正弦波指令（转矩或速度），
利用仿真跑出稳态后的输出幅值和相位，以数值方式求取 Bode 响应。
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import copy
import time

from tuner import (
    get_coeffs_dc_motor_current_regulator,
    get_coeffs_dc_motor_SPEED_regulator,
)
from tutorials_ep6_svpwm import (
    The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental,
)


def extract_mag_phase_dft(t_arr, y_arr, u_arr, freq_Hz):
    """
    用精确频率的 DFT 提取信号在目标频率下的幅值和相位。
    这里增加了严格的“整周期截断”与“去直流偏置”，防止非整周期导致的频谱泄漏（出现坏点）。
    """
    # 截取正好包含最多整数个周期的最新数据
    period = 1.0 / freq_Hz
    total_time = t_arr[-1] - t_arr[0]
    num_periods = max(1, int(np.floor(total_time / period)))
    
    valid_duration = num_periods * period
    valid_start_time = t_arr[-1] - valid_duration
    
    idx = t_arr >= valid_start_time
    t_valid = t_arr[idx]
    y_valid = y_arr[idx]
    u_valid = u_arr[idx]
    
    # 严防直流偏置导致频谱泄露
    y_valid = y_valid - np.mean(y_valid)
    u_valid = u_valid - np.mean(u_valid)
    
    dt = t_valid[1] - t_valid[0]
    
    # 构建复数正弦波基
    basis = np.exp(-1j * 2 * np.pi * freq_Hz * t_valid)
    
    # 积分形式的 DFT (更稳健)
    Y_dft = np.sum(y_valid * basis) * dt
    U_dft = np.sum(u_valid * basis) * dt
    
    if np.abs(U_dft) < 1e-12:
        return 0.0, 0.0
    
    # 计算复数增益 Y / U
    G = Y_dft / U_dft
    
    mag_db = 20 * np.log10(np.abs(G))
    phase_deg = np.degrees(np.angle(G))
    
    return mag_db, phase_deg


def run_single_frequency(d, freq_Hz, zeta, CLBW_Hz, enable_4th_ESO, omega_ob, mode='tracking'):
    """
    跑一次特定频率的正弦注入。
    mode: 'tracking' 注入到 cmd_rpm, 测 speed_rpm
          'disturbance' 注入到 TLoad, 测 speed_rpm
    """
    dd = copy.deepcopy(d)

    # === 手动 tuning ===
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

    # === 构建仿真对象 ===
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

    # === 配置 4阶 ESO ===
    if enable_4th_ESO:
        CTRL.index_separate_speed_estimation = 1
        CTRL.use_disturbance_feedforward_rejection = 1
        CTRL.ell1 = 4 * omega_ob
        CTRL.ell2 = 6 * omega_ob**2
        CTRL.ell3 = 4 * omega_ob**3 * J_s / n_pp
        CTRL.ell4 =     omega_ob**4 * J_s / n_pp

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=dd['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])

    # PI 调节器
    local_Ki_factor = 1.0 if dd.get('CTRL.bool_apply_decoupling_voltages_to_current_regulation', False) else dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False', 10)

    reg_id = The_PID_Regulator(currentKp, currentKp * currentKi * local_Ki_factor, 0.0, 0.0,
                                dd['DC_BUS_VOLTAGE'] / 1.732, dd['DC_BUS_VOLTAGE'] / 1.732, CL_TS)
    reg_iq = The_PID_Regulator(currentKp, currentKp * currentKi * local_Ki_factor, 0.0, 0.0,
                                dd['DC_BUS_VOLTAGE'] / 1.732, dd['DC_BUS_VOLTAGE'] / 1.732, CL_TS)
    
    # 限制放开一点防止超调饱和导致幅相关系非线性
    limit_factor = 5.0 
    reg_speed = The_PID_Regulator(speedKp, speedKp * speedKi, 0.0, 0.0,
                                   limit_factor * dd['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * dd['init_IN'],
                                   limit_factor * dd['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * dd['init_IN'],
                                   VL_TS)

    # === 仿真时间设置 ===
    # 为了扫高频不被阶梯波影响，这里的时间步需要很小。由于ACMSimPyIncremental里面也是离散的控制周期，
    # 我们就以 VL_TS (速度环周期) 作为一个切片。
    slice_dt = VL_TS  
    
    # 假设需要 5 个周期进入稳态，最后取 3 个周期算 DFT
    periods_transient = 5.0
    periods_record    = 3.0
    
    transient_time = max(0.1, periods_transient / freq_Hz)  # 至少 0.1s 让电机启动进入工作点
    record_time = periods_record / freq_Hz
    
    total_time = transient_time + record_time
    n_slices = int(total_time / slice_dt)
    
    record_start_slice = int(transient_time / slice_dt)

    t_arr = []
    y_arr = []
    u_arr = []

    # 工作点
    rpm_0 = 500.0
    T_0 = 0.0
    
    rpm_amp = 10.0   # 注入的转速扰动幅值 （小信号）
    T_amp   = 1.0    # 注入的负载扰动幅值 （小信号）

    for si in range(n_slices):
        t_now = si * slice_dt
        
        if mode == 'tracking':
            u_val = rpm_amp * np.sin(2 * np.pi * freq_Hz * t_now)
            CTRL.cmd_rpm = rpm_0 + u_val
            ACM.TLoad = T_0
        elif mode == 'disturbance':
            u_val = T_amp * np.sin(2 * np.pi * freq_Hz * t_now)
            CTRL.cmd_rpm = rpm_0
            ACM.TLoad = T_0 + u_val
            
        machine_times, watch_data = ACMSimPyIncremental(
            t0=t_now, TIME=slice_dt,
            ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed
        )
        
        # 只记录稳态段的最后数据，因为 ACMSimPyIncremental 会返回包含多个电流环周期的数组
        # 为了精确，我们把这一个 slice 的结果拼进去，如果超出了 transient
        if si >= record_start_slice:
            t_arr.extend(machine_times)
            y_arr.extend(watch_data[1])  # speed_rpm
            
            # 由于在整个切片里输入是不变的或者我们可以用理想正弦输入作为u_arr
            # 为保证相位极其精准，直接按照 machine_times 重新计算当时的理想正弦输入：
            if mode == 'tracking':
                u_arr.extend(rpm_amp * np.sin(2 * np.pi * freq_Hz * np.array(machine_times)))
            elif mode == 'disturbance':
                u_arr.extend(T_amp * np.sin(2 * np.pi * freq_Hz * np.array(machine_times)))

    t_arr = np.array(t_arr)
    y_arr = np.array(y_arr)
    u_arr = np.array(u_arr)
    
    # 剔除由于分片拼接导致的重复时间点 (t0) 
    t_arr, unique_idx = np.unique(t_arr, return_index=True)
    y_arr = y_arr[unique_idx]
    u_arr = u_arr[unique_idx]
    
    # 提取稳态时的交流分量 (减去均值)
    if mode == 'tracking':
        y_ac = y_arr - rpm_0
    else:
        y_ac = y_arr - rpm_0
        
    mag_db, phase_deg = extract_mag_phase_dft(t_arr, y_ac, u_arr, freq_Hz)
    return mag_db, phase_deg


def main():
    parser = argparse.ArgumentParser(description="Numerical Bode Sweep Simulation")
    parser.add_argument('--num-points', type=int, default=20, help="Number of frequency points to sweep")
    parser.add_argument('--freq-start', type=float, default=1.0, help="Start frequency in Hz")
    parser.add_argument('--freq-end', type=float, default=200.0, help="End frequency in Hz")
    parser.add_argument('--zeta', type=float, default=15.0, help="Controller damping parameter zeta")
    parser.add_argument('--clbw', type=float, default=1000.0, help="Current loop bandwidth in Hz")
    parser.add_argument('--omega-ob', type=float, default=200.0, help="4th-order ESO observer bandwidth in rad/s")
    args = parser.parse_args()

    # 导入默认电机参数
    d = {
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

    freqs = np.logspace(np.log10(args.freq_start), np.log10(args.freq_end), args.num_points)
    
    modes = ['tracking', 'disturbance']
    configs = [
        {'name': 'PI Only',        'eso': False, 'color': 'b'},
        {'name': 'PI + 4th ESO',   'eso': True,  'color': 'r'},
    ]
    
    results = {
        'tracking':    { c['name']: {'mag': [], 'phase': []} for c in configs },
        'disturbance': { c['name']: {'mag': [], 'phase': []} for c in configs }
    }

    print(f"Starting sweep: {args.num_points} points from {args.freq_start} to {args.freq_end} Hz.")
    print(f"Params: zeta={args.zeta}, CLBW={args.clbw}Hz, omega_ob={args.omega_ob}rad/s")
    start_time = time.time()

    for i, f in enumerate(freqs):
        print(f"  [{i+1}/{len(freqs)}] Sweeping {f:.1f} Hz...")
        for mode in modes:
            for c in configs:
                m, p = run_single_frequency(d, f, args.zeta, args.clbw, c['eso'], args.omega_ob, mode=mode)
                results[mode][c['name']]['mag'].append(m)
                # 处理相位的连续性
                if len(results[mode][c['name']]['phase']) > 0:
                    prev_p = results[mode][c['name']]['phase'][-1]
                    while p - prev_p > 180: p -= 360
                    while p - prev_p < -180: p += 360
                results[mode][c['name']]['phase'].append(p)
                
    elapsed = time.time() - start_time
    print(f"Simulation completed in {elapsed:.1f} seconds.")

    # 绘图
    plt.style.use('bmh')
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=150)
    
    # 0,0 Tracking Mag
    # 1,0 Tracking Phase
    # 0,1 Disturbance Mag
    # 1,1 Disturbance Phase
    
    for c in configs:
        name = c['name']
        color = c['color']
        
        # Tracking Magnitude
        axes[0,0].semilogx(freqs, results['tracking'][name]['mag'], f'o-{color}', label=name, markersize=4)
        # Tracking Phase
        axes[1,0].semilogx(freqs, results['tracking'][name]['phase'], f'o-{color}', label=name, markersize=4)
        
        # Disturbance Magnitude
        axes[0,1].semilogx(freqs, results['disturbance'][name]['mag'], f'o-{color}', label=name, markersize=4)
        # Disturbance Phase
        axes[1,1].semilogx(freqs, results['disturbance'][name]['phase'], f'o-{color}', label=name, markersize=4)

    # 格式化
    axes[0,0].set_title("Tracking Channel: Mag(speed / cmd)", fontweight='bold')
    axes[0,0].set_ylabel("Magnitude [dB]")
    axes[0,0].grid(True, which='both', ls='--', alpha=0.5)
    axes[0,0].legend()
    
    axes[1,0].set_title("Tracking Channel: Phase", fontweight='bold')
    axes[1,0].set_xlabel("Frequency [Hz]")
    axes[1,0].set_ylabel("Phase [deg]")
    axes[1,0].grid(True, which='both', ls='--', alpha=0.5)
    
    axes[0,1].set_title("Disturbance Channel: Mag(speed / TLoad)", fontweight='bold')
    axes[0,1].set_ylabel("Magnitude [dB]")
    axes[0,1].grid(True, which='both', ls='--', alpha=0.5)
    axes[0,1].legend()

    axes[1,1].set_title("Disturbance Channel: Phase", fontweight='bold')
    axes[1,1].set_xlabel("Frequency [Hz]")
    axes[1,1].set_ylabel("Phase [deg]")
    axes[1,1].grid(True, which='both', ls='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig("fig_sim_bode_sweep.png", dpi=200, bbox_inches='tight')
    print("Saved -> fig_sim_bode_sweep.png")

if __name__ == "__main__":
    main()
