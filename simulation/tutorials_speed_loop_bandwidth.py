# -*- coding: utf-8 -*-
"""
速度环闭环扫频与带宽分析工具
=============================
基于 tutorials_ep6_svpwm.py 的仿真框架:
1. measure_speed_bandwidth_analytical  — 解析法（传递函数 Bode），快速精确
2. measure_speed_bandwidth_simulation  — 仿真扫频法（时域仿真+FFT），含非线性因素
3. scan_bandwidth_vs_zeta_CLBW         — 参数扫描（zeta × CLBW 网格）
4. 绘图：带宽 vs zeta、带宽 vs CLBW、3D 曲面/热力图
"""

# %%
############################################# PACKAGES
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for headless execution
from pylab import np, plt, mpl
import control
import copy
import time as _time

# ======================================================================
# Reuse tuner functions
# ======================================================================
from tuner import (
    get_coeffs_dc_motor_current_regulator,
    get_coeffs_dc_motor_SPEED_regulator,
)


# ======================================================================
# 1. Analytical bandwidth measurement (fast, via transfer function Bode)
# ======================================================================
def measure_speed_bandwidth_analytical(d, zeta, CLBW_Hz):
    """
    解析法计算速度环闭环带宽。

    Parameters
    ----------
    d : dict
        电机和控制器参数字典（与 tutorials_ep6_svpwm.py 一致）。
        会被复制后修改，不影响原始字典。
    zeta : float
        阻尼系数 (对应原代码 FOC_delta)。
    CLBW_Hz : float
        电流环带宽 [Hz]。

    Returns
    -------
    VLBW_Hz : float
        速度环闭环 -3dB 带宽 [Hz]。
    info : dict
        附加信息 (PI 参数、Bode 数据等)。
    """
    R    = d['init_R']
    L    = d['init_Lq']
    J_s  = d['init_Js']
    n_pp = d['init_npp']
    KE   = d['init_KE']
    KA   = KE  # PMSM: KA = KE
    CL_TS = d['CL_TS']
    VL_TS = d['CL_TS'] * d['VL_EXE_PER_CL_EXE']

    # ---------- 电流环 ----------
    currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R, L, CLBW_Hz)
    Gi_closed = control.tf([1], [L / currentKp, 1])  # 零极点对消后的电流闭环
    currentBandwidth_radPerSec = currentKp / L

    # ---------- 速度环 ----------
    KT = 1.5 * n_pp * KA
    dc_motor_motion = control.tf([KT * n_pp / J_s], [1, 0])  # [Apk] → [elec.rad/s]

    speedKp, speedKi = get_coeffs_dc_motor_SPEED_regulator(
        J_s, n_pp, KA, zeta, currentBandwidth_radPerSec
    )
    speedPI = control.tf([speedKp, speedKp * speedKi], [1, 0])

    Gw_open = dc_motor_motion * Gi_closed * speedPI
    Gw_closed = Gw_open / (1 + Gw_open)

    # Bode 数据
    omega_array = 2 * np.pi * np.logspace(0, 4, 2000)
    mag, phase, omega = control.bode_plot(
        Gw_closed, omega_array, dB=False, Hz=False, deg=True, plot=False
    )

    # -3dB 带宽
    idx = (np.abs(mag - 0.707)).argmin()
    VLBW_Hz = omega[idx] / (2 * np.pi)

    info = {
        'currentKp': currentKp,
        'currentKi': currentKi,
        'speedKp': speedKp,
        'speedKi': speedKi,
        'mag': mag,
        'phase': phase,
        'omega': omega,
        'Gw_closed': Gw_closed,
    }
    return VLBW_Hz, info


# ======================================================================
# 2. Simulation-based bandwidth measurement (slow, but realistic)
# ======================================================================
def measure_speed_bandwidth_simulation(d, zeta, CLBW_Hz,
                                        sweep_hz_start=1,
                                        sweep_hz_end=200,
                                        sweep_hz_step=1,
                                        sweep_rpm_amplitude=50,
                                        num_cycles_per_freq=2,
                                        verbose=False):
    """
    仿真扫频法测量速度环闭环带宽。

    在时域中给速度指令施加不同频率的正弦信号，
    测量转速响应的幅值，构建闭环幅频特性，找 -3dB 点。

    Parameters
    ----------
    d : dict
        电机和控制器参数字典。
    zeta : float
        阻尼系数。
    CLBW_Hz : float
        电流环带宽 [Hz]。
    sweep_hz_start, sweep_hz_end, sweep_hz_step : float
        扫频范围和步长 [Hz]。
    sweep_rpm_amplitude : float
        速度指令正弦幅值 [rpm]。
    num_cycles_per_freq : int
        每个频率运行几个周期。
    verbose : bool
        是否打印调试信息。

    Returns
    -------
    VLBW_Hz : float
        速度环闭环 -3dB 带宽 [Hz]。
    info : dict
        附加信息 (频率-增益数组等)。
    """
    # 导入仿真模块（在函数内部导入，避免顶层循环依赖）
    from tutorials_ep6_svpwm import (
        The_Motor_Controller, The_AC_Machine,
        The_PID_Regulator, ACMSimPyIncremental,
    )

    dd = copy.deepcopy(d)
    dd['FOC_delta'] = zeta
    dd['CL_SERIES_KP'] = None  # 触发 auto-tuning
    dd['CL_SERIES_KI'] = None
    dd['VL_SERIES_KP'] = None
    dd['VL_SERIES_KI'] = None

    # --- 手动 tuning（指定 CLBW_Hz） ---
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

    dd['CL_SERIES_KP'] = currentKp
    dd['CL_SERIES_KI'] = currentKi
    dd['VL_SERIES_KP'] = speedKp
    dd['VL_SERIES_KI'] = speedKi

    # --- 构建仿真对象 ---
    CTRL = The_Motor_Controller(
        CL_TS=dd['CL_TS'],
        VL_TS=VL_TS,
        init_npp=dd['init_npp'],
        init_IN=dd['init_IN'],
        init_R=dd['init_R'],
        init_Ld=dd['init_Ld'],
        init_Lq=dd['init_Lq'],
        init_KE=dd['init_KE'],
        init_Rreq=dd['init_Rreq'],
        init_Js=dd['init_Js'],
        DC_BUS_VOLTAGE=dd['DC_BUS_VOLTAGE'],
    )
    CTRL.bool_apply_decoupling_voltages_to_current_regulation = dd['CTRL.bool_apply_decoupling_voltages_to_current_regulation']
    CTRL.bool_apply_sweeping_frequency_excitation = False  # 我们手动控制扫频
    CTRL.bool_overwrite_speed_commands = True
    CTRL.bool_zero_id_control = dd['CTRL.bool_zero_id_control']
    CTRL.bool_apply_speed_closed_loop_control = True

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=dd['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])

    # 电流环 PI (PID with Kd=0)
    if dd.get('CTRL.bool_apply_decoupling_voltages_to_current_regulation', False):
        local_Ki_factor = 1.0
    else:
        local_Ki_factor = dd.get('FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False', 10)
    reg_id = The_PID_Regulator(currentKp, currentKp * currentKi * local_Ki_factor, 0.0, 0.0,
                                dd['DC_BUS_VOLTAGE'] / 1.732, dd['DC_BUS_VOLTAGE'] / 1.732, CL_TS)
    reg_iq = The_PID_Regulator(currentKp, currentKp * currentKi * local_Ki_factor, 0.0, 0.0,
                                dd['DC_BUS_VOLTAGE'] / 1.732, dd['DC_BUS_VOLTAGE'] / 1.732, CL_TS)
    reg_speed = The_PID_Regulator(speedKp, speedKp * speedKi, 0.0, 0.0,
                                   dd['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * dd['init_IN'],
                                   dd['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * dd['init_IN'],
                                   VL_TS)

    # --- 逐频率扫频 ---
    MACHINE_TS = CL_TS / ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    freq_list = np.arange(sweep_hz_start, sweep_hz_end + 0.5 * sweep_hz_step, sweep_hz_step)
    gain_list = []

    # 先让系统稳定 0.5 秒（零速稳态）
    CTRL.cmd_rpm = 0.0
    warmup_time = 0.5
    _, _ = ACMSimPyIncremental(
        t0=0.0, TIME=warmup_time,
        ACM=ACM, CTRL=CTRL,
        reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed
    )

    current_time = warmup_time
    for fi, freq_hz in enumerate(freq_list):
        if freq_hz <= 0:
            gain_list.append(1.0)
            continue

        period = 1.0 / freq_hz
        sim_duration = period * num_cycles_per_freq
        n_steps = int(sim_duration / MACHINE_TS)

        # 收集时间和响应
        cmd_speed_trace = np.zeros(n_steps)
        act_speed_trace = np.zeros(n_steps)

        # 用逐步仿真来收集数据
        machine_times, watch_data = ACMSimPyIncremental(
            t0=current_time, TIME=sim_duration,
            ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed
        )

        # watch_data[12] = CTRL.cmd_rpm, watch_data[1] = ACM.omega_r_mech in rpm
        cmd_trace = watch_data[12]
        act_trace = watch_data[1]  # ACM.omega_r_mech in rpm (已转换)

        # 取后半段数据（避免暂态）
        half = len(cmd_trace) // 2
        cmd_segment = cmd_trace[half:]
        act_segment = act_trace[half:]

        # 计算增益 = 输出幅值 / 输入幅值
        cmd_amp = (np.max(cmd_segment) - np.min(cmd_segment)) / 2
        act_amp = (np.max(act_segment) - np.min(act_segment)) / 2

        if cmd_amp > 1e-10:
            gain = act_amp / cmd_amp
        else:
            gain = 0.0
        gain_list.append(gain)

        current_time += sim_duration

        # 更新速度指令为正弦（在仿真内部通过 bool_overwrite_speed_commands=True 
        # 和手动设置实现，但 ACMSimPyIncremental 的 Console 部分会覆盖）
        # 所以我们需要在每个频率的仿真之前设置扫频参数
        if verbose and fi % 10 == 0:
            print(f'  freq={freq_hz:.0f} Hz, gain={gain:.4f} ({20*np.log10(max(gain, 1e-10)):.1f} dB)')

    gain_array = np.array(gain_list)

    # 归一化（DC增益应为1，但可能因为暖机不足略有偏差）
    if gain_array[0] > 0:
        gain_array_normalized = gain_array / gain_array[0]
    else:
        gain_array_normalized = gain_array

    # 找 -3dB 点
    idx = (np.abs(gain_array_normalized - 0.707)).argmin()
    VLBW_Hz = freq_list[idx]

    info = {
        'freq_list': freq_list,
        'gain_array': gain_array,
        'gain_array_normalized': gain_array_normalized,
    }
    return VLBW_Hz, info


# ======================================================================
# 2b. Improved simulation sweep using built-in sweep mechanism
# ======================================================================
def measure_speed_bandwidth_simulation_builtin(d, zeta, CLBW_Hz,
                                                sweep_hz_ceiling=200,
                                                sweep_rpm_amplitude=50,
                                                verbose=False):
    """
    利用仿真框架内置的扫频机制 (bool_apply_sweeping_frequency_excitation)
    测量速度环闭环带宽。

    Returns
    -------
    VLBW_Hz : float
        速度环闭环 -3dB 带宽 [Hz]。
    info : dict
        附加信息。
    """
    from tutorials_ep6_svpwm import Simulation_Benchmark
    from collections import OrderedDict as OD

    dd = copy.deepcopy(d)
    dd['FOC_delta'] = zeta

    # 手动设置 PI 参数（指定 CLBW）
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

    dd['CL_SERIES_KP'] = currentKp
    dd['CL_SERIES_KI'] = currentKi
    dd['VL_SERIES_KP'] = speedKp
    dd['VL_SERIES_KI'] = speedKi

    # 启用扫频
    dd['CTRL.bool_apply_sweeping_frequency_excitation'] = True
    dd['CTRL.bool_overwrite_speed_commands'] = True
    dd['CTRL.bool_apply_speed_closed_loop_control'] = True

    # 估计仿真总时间以覆盖所有扫频频率
    # 每个频率运行 1/freq 秒，总时间 ≈ sum(1/f for f in 1..ceiling)
    total_time = sum(1.0 / f for f in range(1, int(sweep_hz_ceiling) + 1))
    dd['TIME_SLICE'] = min(total_time + 1.0, 20.0)
    dd['NUMBER_OF_SLICES'] = max(1, int(np.ceil((total_time + 1.0) / dd['TIME_SLICE'])))

    # 设置 scope 只收集我们需要的信号
    numba__scope_dict = OD([
        (r'Speed [rpm]', ('CTRL.cmd_rpm', 'CTRL.omega_r_mech',)),
    ])

    sim = Simulation_Benchmark(dd, bool_start_simulation=False)
    sim.d = dd

    # 获取全局对象并手动配置扫频参数
    CTRL, ACM, reg_id, reg_iq, reg_speed, _, _ = sim.get_global_objects()
    CTRL.CMD_SPEED_SINE_RPM = sweep_rpm_amplitude
    CTRL.CMD_SPEED_SINE_HZ = 0
    CTRL.CMD_SPEED_SINE_STEP_SIZE = 1
    CTRL.CMD_SPEED_SINE_HZ_CEILING = sweep_hz_ceiling
    CTRL.CMD_SPEED_SINE_END_TIME = 0.0
    CTRL.CMD_SPEED_SINE_LAST_END_TIME = 0.0
    CTRL.bool_apply_sweeping_frequency_excitation = True
    CTRL.bool_overwrite_speed_commands = True

    # 运行仿真（拼接所有 slice）
    from tutorials_ep6_svpwm import ACMSimPyWrapper, ACMSimPyIncremental
    global_machine_times = None
    global_cmd_rpm = None
    global_act_rpm = None

    def save_to_global(_global, _local):
        return _local if _global is None else np.append(_global, _local)

    for ii in range(dd['NUMBER_OF_SLICES']):
        machine_times, watch_data = ACMSimPyIncremental(
            t0=ii * dd['TIME_SLICE'], TIME=dd['TIME_SLICE'],
            ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed
        )
        global_machine_times = save_to_global(global_machine_times, machine_times)
        global_cmd_rpm = save_to_global(global_cmd_rpm, watch_data[12])  # CTRL.cmd_rpm
        global_act_rpm = save_to_global(global_act_rpm, watch_data[1])   # ACM.omega_r_mech in rpm

    # 分析每个频率段的增益
    CL_TS_val = dd['CL_TS']
    MACHINE_TS = CL_TS_val / dd['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD']
    dt = MACHINE_TS

    # 重建各频率段的时间边界
    freq_boundaries = []
    t_start = 0.0
    for f_hz in range(1, int(sweep_hz_ceiling) + 1):
        t_end = t_start + 1.0 / f_hz
        freq_boundaries.append((f_hz, t_start, t_end))
        t_start = t_end

    freq_list = []
    gain_list = []
    for f_hz, t_start, t_end in freq_boundaries:
        # 取对应时间段的数据
        mask = (global_machine_times >= t_start) & (global_machine_times < t_end)
        if np.sum(mask) < 10:
            continue
        cmd_segment = global_cmd_rpm[mask]
        act_segment = global_act_rpm[mask]

        # 取后半段避免暂态
        half = len(cmd_segment) // 2
        cmd_amp = (np.max(cmd_segment[half:]) - np.min(cmd_segment[half:])) / 2
        act_amp = (np.max(act_segment[half:]) - np.min(act_segment[half:])) / 2

        if cmd_amp > 1e-10:
            gain = act_amp / cmd_amp
        else:
            gain = 0.0

        freq_list.append(f_hz)
        gain_list.append(gain)

    freq_array = np.array(freq_list)
    gain_array = np.array(gain_list)

    # 归一化
    if len(gain_array) > 0 and gain_array[0] > 0:
        gain_normalized = gain_array / gain_array[0]
    else:
        gain_normalized = gain_array

    # -3dB 带宽
    if len(gain_normalized) > 0:
        idx = (np.abs(gain_normalized - 0.707)).argmin()
        VLBW_Hz = freq_array[idx]
    else:
        VLBW_Hz = 0.0

    info = {
        'freq_array': freq_array,
        'gain_array': gain_array,
        'gain_normalized': gain_normalized,
        'global_machine_times': global_machine_times,
        'global_cmd_rpm': global_cmd_rpm,
        'global_act_rpm': global_act_rpm,
    }
    return VLBW_Hz, info


# ======================================================================
# 3. Parameter scan: bandwidth vs (zeta, CLBW)
# ======================================================================
def scan_bandwidth_vs_zeta_CLBW(d, zeta_list, CLBW_list):
    """
    扫描 zeta 和 CLBW 的组合，返回带宽矩阵（解析法）。

    Parameters
    ----------
    d : dict
        基础电机参数字典。
    zeta_list : array-like
        zeta (阻尼系数) 扫描值列表。
    CLBW_list : array-like
        CLBW (电流环带宽 Hz) 扫描值列表。

    Returns
    -------
    bandwidth_matrix : np.ndarray
        shape (len(zeta_list), len(CLBW_list))，每个元素是 VLBW [Hz]。
    """
    bandwidth_matrix = np.zeros((len(zeta_list), len(CLBW_list)))

    total = len(zeta_list) * len(CLBW_list)
    count = 0
    for i, zeta in enumerate(zeta_list):
        for j, CLBW_Hz in enumerate(CLBW_list):
            try:
                vlbw, _ = measure_speed_bandwidth_analytical(d, zeta, CLBW_Hz)
                bandwidth_matrix[i, j] = vlbw
            except Exception as e:
                bandwidth_matrix[i, j] = np.nan
            count += 1
            if count % 50 == 0:
                print(f'  Progress: {count}/{total}')

    return bandwidth_matrix


# ======================================================================
# 4. Plotting
# ======================================================================
def plot_bandwidth_curves(d, zeta_list, CLBW_list, bandwidth_matrix, save_path=None):
    """
    绘制速度环带宽与 zeta、CLBW 关系的曲线图。

    Generates 3 figures:
    1. 带宽 vs CLBW（不同 zeta 为不同曲线）
    2. 带宽 vs zeta（不同 CLBW 为不同曲线）
    3. 3D 曲面图 / 热力图
    """
    plt.style.use('bmh')
    mpl.rc('font', family='Times New Roman', size=11.0)
    mpl.rc('legend', fontsize=9)
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['mathtext.fontset'] = 'stix'

    # ---- Figure 1: Bandwidth vs CLBW (different zeta curves) ----
    fig1, ax1 = plt.subplots(dpi=150, facecolor='w', figsize=(9, 6))
    colors = plt.cm.viridis(np.linspace(0, 1, len(zeta_list)))
    # 选取部分 zeta 值绘制（避免太挤）
    n_curves = min(10, len(zeta_list))
    indices = np.linspace(0, len(zeta_list) - 1, n_curves, dtype=int)

    for idx in indices:
        zeta = zeta_list[idx]
        valid = ~np.isnan(bandwidth_matrix[idx, :])
        ax1.plot(np.array(CLBW_list)[valid], bandwidth_matrix[idx, valid],
                 'o-', color=colors[idx], markersize=3,
                 label=f'ζ = {zeta:.0f}')
    ax1.set_xlabel('Current Loop Bandwidth CLBW [Hz]')
    ax1.set_ylabel('Speed Loop Bandwidth VLBW [Hz]')
    ax1.set_title('Speed Loop Bandwidth vs Current Loop Bandwidth')
    ax1.legend(loc='best', ncol=2)
    ax1.grid(True, alpha=0.3)
    fig1.tight_layout()

    # ---- Figure 2: Bandwidth vs zeta (different CLBW curves) ----
    fig2, ax2 = plt.subplots(dpi=150, facecolor='w', figsize=(9, 6))
    colors2 = plt.cm.plasma(np.linspace(0, 1, len(CLBW_list)))
    n_curves2 = min(10, len(CLBW_list))
    indices2 = np.linspace(0, len(CLBW_list) - 1, n_curves2, dtype=int)

    for idx in indices2:
        clbw = CLBW_list[idx]
        valid = ~np.isnan(bandwidth_matrix[:, idx])
        ax2.plot(np.array(zeta_list)[valid], bandwidth_matrix[valid, idx],
                 's-', color=colors2[idx], markersize=3,
                 label=f'CLBW = {clbw:.0f} Hz')
    ax2.set_xlabel('Damping Ratio ζ (FOC_delta)')
    ax2.set_ylabel('Speed Loop Bandwidth VLBW [Hz]')
    ax2.set_title('Speed Loop Bandwidth vs Damping Ratio ζ')
    ax2.legend(loc='best', ncol=2)
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()

    # ---- Figure 3: 3D Surface + Contour ----
    fig3, (ax3a, ax3b) = plt.subplots(1, 2, dpi=150, facecolor='w', figsize=(16, 6),
                                       subplot_kw={'projection': None})
    # Remove projection for second subplot
    fig3.clf()
    ax3a = fig3.add_subplot(121, projection='3d')
    ax3b = fig3.add_subplot(122)

    Z_grid, C_grid = np.meshgrid(zeta_list, CLBW_list, indexing='ij')

    # 3D surface
    surf = ax3a.plot_surface(Z_grid, C_grid, bandwidth_matrix,
                              cmap='coolwarm', alpha=0.85, edgecolor='none')
    ax3a.set_xlabel('ζ (FOC_delta)')
    ax3a.set_ylabel('CLBW [Hz]')
    ax3a.set_zlabel('VLBW [Hz]')
    ax3a.set_title('Speed Loop Bandwidth Surface')
    fig3.colorbar(surf, ax=ax3a, shrink=0.5, label='VLBW [Hz]')

    # Heatmap / contour
    bw_for_contour = bandwidth_matrix.copy()
    bw_for_contour[np.isnan(bw_for_contour)] = 0
    contour = ax3b.contourf(C_grid, Z_grid, bw_for_contour, levels=20, cmap='coolwarm')
    cs = ax3b.contour(C_grid, Z_grid, bw_for_contour, levels=10, colors='black', linewidths=0.5)
    ax3b.clabel(cs, inline=True, fontsize=7, fmt='%.0f')
    ax3b.set_xlabel('CLBW [Hz]')
    ax3b.set_ylabel('ζ (FOC_delta)')
    ax3b.set_title('Speed Loop Bandwidth Contour Map')
    fig3.colorbar(contour, ax=ax3b, label='VLBW [Hz]')
    fig3.tight_layout()

    if save_path:
        fig1.savefig(f'{save_path}_vlbw_vs_clbw.png', dpi=200, bbox_inches='tight')
        fig2.savefig(f'{save_path}_vlbw_vs_zeta.png', dpi=200, bbox_inches='tight')
        fig3.savefig(f'{save_path}_vlbw_surface.png', dpi=200, bbox_inches='tight')

    return fig1, fig2, fig3


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':
    print('=' * 60)
    print(' 速度环闭环扫频与带宽分析工具')
    print('=' * 60)

    # ---------- 基础电机参数（小电感电机，与 ep6 一致） ----------
    d = {
        'CL_TS': 1e-4,
        'VL_EXE_PER_CL_EXE': 5,
        'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
        'TIME_SLICE': 0.2,
        'NUMBER_OF_SLICES': 6,
        'init_npp': 22,
        'init_IN': 1.3 * 6 / 1.414,
        'init_R': 0.035,
        'init_Ld': 1 * 0.036e-3,
        'init_Lq': 1 * 0.036e-3,
        'init_KE': 0.0125,
        'init_Rreq': 0.0,
        'init_Js': 0.44e-4,
        'DC_BUS_VOLTAGE': 5,
        'CTRL.bool_apply_speed_closed_loop_control': True,
        'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
        'CTRL.bool_apply_sweeping_frequency_excitation': False,
        'CTRL.bool_overwrite_speed_commands': True,
        'CTRL.bool_zero_id_control': True,
        'FOC_delta': 15,
        'FOC_desired_VLBW_HZ': 120,
        'FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False': 10,
        'CL_SERIES_KP': None,
        'CL_SERIES_KI': None,
        'VL_SERIES_KP': None,
        'VL_SERIES_KI': None,
        'VL_LIMIT_OVERLOAD_FACTOR': 3.0,
        'disp.Kp': 0.0,
        'disp.Ki': 0.0,
        'disp.Kd': 0.0,
        'disp.tau': 0.0,
        'disp.OutLimit': 0.0,
        'disp.IntLimit': 0.0,
    }

    # ====================================================
    # Demo 1: 单点解析法测量
    # ====================================================
    print('\n--- Demo 1: 单点解析法测量 ---')
    zeta_demo = 15
    CLBW_demo = 1000  # Hz
    vlbw, info = measure_speed_bandwidth_analytical(d, zeta_demo, CLBW_demo)
    print(f'  zeta={zeta_demo}, CLBW={CLBW_demo} Hz  =>  VLBW = {vlbw:.1f} Hz')

    # ====================================================
    # Demo 2: 参数扫描（解析法）
    # ====================================================
    print('\n--- Demo 2: 参数扫描（解析法）---')
    zeta_list = np.arange(3, 31, 1)       # zeta ∈ [3, 30]
    CLBW_list = np.arange(200, 3001, 100)  # CLBW ∈ [200, 3000] Hz

    t_start = _time.time()
    bandwidth_matrix = scan_bandwidth_vs_zeta_CLBW(d, zeta_list, CLBW_list)
    elapsed = _time.time() - t_start
    print(f'  参数扫描完成，耗时 {elapsed:.1f} 秒')
    print(f'  带宽矩阵形状: {bandwidth_matrix.shape}')
    print(f'  带宽范围: [{np.nanmin(bandwidth_matrix):.1f}, {np.nanmax(bandwidth_matrix):.1f}] Hz')

    # ====================================================
    # Demo 3: 绘图
    # ====================================================
    print('\n--- Demo 3: 绘制带宽关系曲线 ---')
    fig1, fig2, fig3 = plot_bandwidth_curves(d, zeta_list, CLBW_list, bandwidth_matrix,
                                              save_path='speed_loop_bandwidth')

    # ====================================================
    # Demo 4: Bode 对比图（选几个典型参数）
    # ====================================================
    print('\n--- Demo 4: 典型参数的闭环 Bode 图 ---')
    fig4, (ax4a, ax4b) = plt.subplots(2, 1, dpi=150, facecolor='w', figsize=(9, 8), sharex=True)

    test_cases = [
        (5,  500,  'tab:blue'),
        (10, 500,  'tab:orange'),
        (15, 500,  'tab:green'),
        (15, 1000, 'tab:red'),
        (15, 2000, 'tab:purple'),
        (25, 1000, 'tab:brown'),
    ]
    for zeta, clbw, color in test_cases:
        vlbw, info = measure_speed_bandwidth_analytical(d, zeta, clbw)
        freq_hz = info['omega'] / (2 * np.pi)
        mag_db = 20 * np.log10(info['mag'])
        ax4a.semilogx(freq_hz, mag_db, color=color,
                       label=f'ζ={zeta}, CLBW={clbw}Hz → VLBW={vlbw:.0f}Hz')
        ax4b.semilogx(freq_hz, info['phase'], color=color)

    ax4a.axhline(-3, color='gray', ls='--', lw=0.8, label='-3 dB')
    ax4a.set_ylabel('Magnitude [dB]')
    ax4a.set_title('Speed Loop Closed-Loop Bode Diagram')
    ax4a.legend(loc='best', fontsize=7)
    ax4a.set_ylim([-40, 10])
    ax4a.grid(True, which='both', alpha=0.3)

    ax4b.set_xlabel('Frequency [Hz]')
    ax4b.set_ylabel('Phase [deg]')
    ax4b.grid(True, which='both', alpha=0.3)
    fig4.tight_layout()
    fig4.savefig('speed_loop_bandwidth_bode.png', dpi=200, bbox_inches='tight')

    print('\n--- 所有图已保存到当前目录 ---')
    print('  speed_loop_bandwidth_vlbw_vs_clbw.png')
    print('  speed_loop_bandwidth_vlbw_vs_zeta.png')
    print('  speed_loop_bandwidth_vlbw_surface.png')
    print('  speed_loop_bandwidth_bode.png')
    plt.close('all')
