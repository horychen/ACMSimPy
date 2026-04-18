# -*- coding: utf-8 -*-
"""

TL;DR
> python .\demo_4th_order_ESO.py --time 13 --ki 0 --scenarios 4_ON --loads sweep --omega_ob 200.0 --cmd_rpm 0.0

四阶扩展状态观测器 (4th-Order ESO) 演示
========================================
基于 tutorials_ep6_svpwm.py 的仿真框架，演示:
1. 3阶 vs 4阶位置观测器在斜坡负载转矩下的扰动估计
2. 前馈补偿 (feedforward) 对速度跟踪的改善
3. 阶跃负载 vs 斜坡负载的对比

运行 4 个场景:
  (a) 3阶观测器, 无前馈 — baseline
  (b) 3阶观测器, 有前馈 — 能跟踪阶跃但不能跟踪斜坡
  (c) 4阶观测器, 无前馈 — 更好的估计，无补偿
  (d) 4阶观测器, 有前馈 — 能跟踪斜坡负载

附加: 阶跃负载场景对比
"""

# %%
############################################# PACKAGES
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
from pylab import np, plt, mpl
import copy
import time as _time

from tutorials_ep6_svpwm import (
    Simulation_Benchmark,
    ACMSimPyIncremental,
)
from collections import OrderedDict as OD
import tuner


# ======================================================================
# Base motor parameter dictionary (小电感电机, consistent with ep6)
# ======================================================================
def get_base_d():
    """返回基础参数字典。"""
    return {
        'CL_TS': 1e-4,
        'VL_EXE_PER_CL_EXE': 5,
        'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
        'TIME_SLICE': 0.01,
        'NUMBER_OF_SLICES': 150,  # 1.5 秒仿真 (0.01s × 150)
        'init_npp': 22,
        'init_IN': 1.3 * 6 / 1.414,
        'init_R': 0.035,
        'init_Ld': 0.036e-3,
        'init_Lq': 0.036e-3,
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
        'user_system_input_code': 'CTRL.cmd_rpm = 50',
    }


# ======================================================================
# Core function: run one ESO scenario
# ======================================================================
def run_scenario(d_base, observer_order=4, feedforward_on=True, load_type='step', omega_ob=100.0, cmd_rpm=50.0, verbose=True, virtual_inertia_kd=0.0):
    """
    运行单次仿真。返回时间序列和关键波形数据。

    使用 Simulation_Benchmark 框架确保 exec() 和 numba JIT 正确工作。

    Parameters
    ----------
    d_base : dict
        基础参数字典 (会被复制, 不修改原始)。
    observer_order : int
        观测器阶数: 3 或 4。
    feedforward_on : bool
        是否启用扰动前馈补偿。
    load_type : str
        'ramp' — 斜坡负载; 'step' — 阶跃负载。
    omega_ob : float
        观测器带宽 [rad/s]。
    verbose : bool
        是否打印过程信息。

    Returns
    -------
    result : dict
        包含时间、速度、负载估计等数据的字典。
    """
    dd = copy.deepcopy(d_base)

    # --- 设置负载施加代码 ---
    TIME_SLICE = dd['TIME_SLICE']
    dd['user_system_input_code'] = f"CTRL.cmd_rpm = {cmd_rpm}\n"

    # --- 先自动整定 PI 参数 ---
    tuner.tunner_wrapper(dd)
    if dd.get('override_speed_ki') is not None:
        dd['VL_SERIES_KI'] = dd['override_speed_ki']

    label = f'ESO-{observer_order}, FF={"ON" if feedforward_on else "OFF"}'
    if verbose:
        print(f'\n--- Running: {label}, load={load_type} ---')

    # --- 定义需要收集的信号 ---
    numba__scope_dict = OD([
        (r'Speed [rpm]',     ('CTRL.cmd_rpm', 'CTRL.omega_r_mech',)),
        (r'Torque [Nm]',     ('ACM.TLoad', 'CTRL.xS[2]',)),
        (r'Rate [Nm/s]',     ('CTRL.xS[3]',)),
        (r'iq [A]',          ('CTRL.cmd_idq[1]',)),
        (r'Speed err [rpm]', ('CTRL.cmd_rpm - CTRL.omega_r_mech',)),
    ])

    # --- 使用 Simulation_Benchmark 但不自动启动仿真 ---
    sim = Simulation_Benchmark(dd, bool_start_simulation=False)

    # 获取全局对象 (包含已整定的 PI 参数)
    CTRL, ACM, reg_id, reg_iq, reg_speed, reg_dispX, reg_dispY = sim.get_global_objects()
    
    if virtual_inertia_kd > 0.0:
        from tutorials_ep6_svpwm import The_PID_Regulator
        local_Kp = dd['VL_SERIES_KP']
        local_Ki = dd['VL_SERIES_KP'] * dd['VL_SERIES_KI']
        local_Kd = virtual_inertia_kd
        local_tau = 0.0005 # differentiator filtering time constant
        local_OutLimit = dd['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * dd['init_IN']
        local_IntLimit = local_OutLimit
        reg_speed = The_PID_Regulator(local_Kp, local_Ki, local_Kd, local_tau, local_OutLimit, local_IntLimit, CTRL.VL_TS)

    # --- 配置连续且高分辨率的负载模型 ---
    if load_type == 'ramp':
        ACM.custom_load_type = 2
        ACM.custom_load_params[0] = 5 * TIME_SLICE  # t_start
        ACM.custom_load_params[1] = 0.2             # slope
    elif load_type == 'step':
        ACM.custom_load_type = 1
        ACM.custom_load_params[0] = 10 * TIME_SLICE # t_step
        ACM.custom_load_params[1] = 0.15            # amplitude
    elif load_type == 'sine':
        ACM.custom_load_type = 3
        ACM.custom_load_params[0] = 5 * TIME_SLICE  # t_start
        ACM.custom_load_params[1] = 0.15            # amplitude
        ACM.custom_load_params[2] = 2.0             # freq
    elif load_type == 'parabola':
        ACM.custom_load_type = 4
        ACM.custom_load_params[0] = 5 * TIME_SLICE  # t_start
        ACM.custom_load_params[1] = 0.1             # coefficient
    elif load_type == 'sweep':
        ACM.custom_load_type = 5
        T_start = 5 * TIME_SLICE
        T_end = dd['NUMBER_OF_SLICES'] * TIME_SLICE
        T_span = T_end - T_start if (T_end - T_start) > 0.1 else 1.0
        ACM.custom_load_params[0] = T_start
        ACM.custom_load_params[1] = 0.15            # amplitude
        ACM.custom_load_params[2] = 0.2             # f0 (hz)
        ACM.custom_load_params[3] = 200.0           # f1 (hz)
        ACM.custom_load_params[4] = T_span          # span
        ACM.custom_load_params[5] = 0.2             # state: current_f
        ACM.custom_load_params[6] = T_start         # state: current_cycle_start

    # --- 配置观测器 ---
    CTRL.index_separate_speed_estimation = 1  # 使用观测器估计速度

    n_pp = dd['init_npp']
    J_s = dd['init_Js']

    # 设置观测器增益
    CTRL.ell1 = 0.0
    CTRL.ell2 = 0.0
    CTRL.ell3 = 0.0
    CTRL.ell4 = 0.0

    if observer_order == 3:
        CTRL.ell1 = 3 * omega_ob
        CTRL.ell2 = 3 * omega_ob**2
        CTRL.ell3 =     omega_ob**3 * J_s / n_pp
        CTRL.ell4 = 0.0
    elif observer_order == 4:
        CTRL.ell1 = 4 * omega_ob
        CTRL.ell2 = 6 * omega_ob**2
        CTRL.ell3 = 4 * omega_ob**3 * J_s / n_pp
        CTRL.ell4 =     omega_ob**4 * J_s / n_pp
    else:
        raise ValueError(f"observer_order must be 3 or 4, got {observer_order}")

    # --- 配置前馈补偿 ---
    if feedforward_on:
        CTRL.use_disturbance_feedforward_rejection = 1
    else:
        CTRL.use_disturbance_feedforward_rejection = 0

    if verbose:
        print(f'    ell1={CTRL.ell1:.1f}, ell2={CTRL.ell2:.1f}, '
              f'ell3={CTRL.ell3:.6f}, ell4={CTRL.ell4:.1f}')
        print(f'    feedforward={CTRL.use_disturbance_feedforward_rejection}')

    # --- 手动替换 sim 的全局对象并运行仿真循环 ---
    # (模仿 start_simulation_slices 的逻辑, 但使用我们配置好的 CTRL)
    sim.CTRL = CTRL
    sim.ACM = ACM
    sim.reg_id = reg_id
    sim.reg_iq = reg_iq
    sim.reg_speed = reg_speed

    from tutorials_ep6_svpwm import ACMSimPyWrapper

    global_trace_names = []
    max_number_of_traces = 0
    for ylabel, trace_names in numba__scope_dict.items():
        for name in trace_names:
            max_number_of_traces += 1
        for trace_index, name in enumerate(trace_names):
            global_trace_names.append(name)

    global_arrays = [None] * max_number_of_traces
    global_machine_times = None

    def save_to_global(_global, _local):
        return _local if _global is None else np.append(_global, _local)

    for ii in range(dd['NUMBER_OF_SLICES']):
        # 执行用户的系统输入代码 (与 Simulation_Benchmark 一致)
        exec(dd['user_system_input_code'])

        machine_times, numba__waveforms_dict = ACMSimPyWrapper(
            numba__scope_dict,
            t0=ii * dd['TIME_SLICE'], TIME=dd['TIME_SLICE'],
            ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed
        )

        global_machine_times = save_to_global(global_machine_times, machine_times)
        global_index = 0
        for ylabel in numba__scope_dict.keys():
            for trace_index, local_trace_data in enumerate(numba__waveforms_dict[ylabel]):
                global_arrays[global_index] = save_to_global(global_arrays[global_index], local_trace_data)
                global_index += 1

    # 构建 gdd 全局数据字典
    gdd = OD()
    for name, array in zip(global_trace_names, global_arrays):
        gdd[name] = array

    if verbose:
        print(f'    Simulation done. Keys: {list(gdd.keys())}')

    result = {
        'label': label,
        'observer_order': observer_order,
        'feedforward_on': feedforward_on,
        'load_type': load_type,
        'time': global_machine_times,
        'cmd_rpm': gdd['CTRL.cmd_rpm'],
        'act_rpm': gdd['CTRL.omega_r_mech'],
        'TLoad': gdd['ACM.TLoad'],
        'TL_est': gdd['CTRL.xS[2]'],     # xS[2]: 估计的扰动 (约等于 -TLoad)
        'pT_est': gdd['CTRL.xS[3]'],     # xS[3]: 估计的扰动变化率
        'iq_cmd': gdd['CTRL.cmd_idq[1]'],
        'speed_err': gdd['CTRL.cmd_rpm - CTRL.omega_r_mech'],
    }

    if verbose:
        half = len(result['speed_err']) // 2
        rms_err = np.sqrt(np.mean(result['speed_err'][half:]**2))
        print(f'    Speed RMS error (last half): {rms_err:.4f} rpm')

    return result


# ======================================================================
# Plotting functions
# ======================================================================
def plot_overlay(results, title_suffix='', save_path=None):
    """
    将所有场景叠加在同一组子图上 (便于直观对比)。
    """
    plt.style.use('bmh')
    mpl.rc('font', family='Times New Roman', size=10.0)
    mpl.rc('legend', fontsize=8)
    mpl.rcParams['lines.linewidth'] = 1.0
    mpl.rcParams['mathtext.fontset'] = 'stix'

    fig, axes = plt.subplots(nrows=4, ncols=1, dpi=150, facecolor='w',
                              figsize=(10, 14), sharex=True)

    style_map = {
        (3, False): ('tab:blue',   '--',  '3rd, FF OFF'),
        (3, True):  ('tab:orange', '-.',   '3rd, FF ON'),
        (4, False): ('tab:green',  '--',  '4th, FF OFF'),
        (4, True):  ('tab:red',    '-',   '4th, FF ON'),
    }

    first_load = True
    first = True
    for res in results:
        t = res['time']
        key = (res['observer_order'], res['feedforward_on'])
        color, ls, lbl = style_map.get(key, ('black', '-', ''))

        # Row 0: Speed
        ax = axes[0]
        if first:
            ax.plot(t, res['cmd_rpm'], 'k--', lw=0.8, label=r'$\omega^*$')
            first = False
        ax.plot(t, res['act_rpm'], color=color, ls=ls, lw=0.9, label=lbl)
        ax.set_ylabel('Speed [rpm]')
        ax.set_title(f'Speed Tracking Comparison {title_suffix}', fontsize=11)

        # Row 1: Speed error
        ax = axes[1]
        ax.plot(t, res['speed_err'], color=color, ls=ls, lw=0.7, label=lbl)
        ax.set_ylabel('Speed Error [rpm]')

        # Row 2: Load estimation
        ax = axes[2]
        if first_load:  # 只画一次真实负载
            ax.plot(t, res['TLoad'], 'k-', lw=1.2, label=r'$T_L$ actual')
            first_load = False
        ax.plot(t, -res['TL_est'], color=color, ls=ls, lw=0.8,
                label=f'$-\\hat{{x}}_2$ ({lbl})')
        ax.set_ylabel('Load Torque [Nm]')

        # Row 3: pT estimation
        ax = axes[3]
        ax.plot(t, -res['pT_est'], color=color, ls=ls, lw=0.7,
                label=f'$-\\hat{{x}}_3$ ({lbl})')
        ax.set_ylabel(r'$d\hat{T}_L/dt$ [Nm/s]')

    for ax in axes:
        ax.legend(loc='best', fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel('Time [s]')
    axes[1].axhline(0, color='gray', lw=0.5)

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f'  Figure saved: {save_path}')

    return fig


def plot_combined(target_results, save_path=None, cli_suffix=''):
    """
    综合对比图: 动态并排显示多个负载对比。
    """
    plt.style.use('bmh')
    mpl.rc('font', family='Times New Roman', size=10.0)
    mpl.rc('legend', fontsize=7)
    mpl.rcParams['lines.linewidth'] = 1.0
    mpl.rcParams['mathtext.fontset'] = 'stix'

    ncols = len(target_results)
    if ncols == 0: return None
    
    fig, axes = plt.subplots(4, ncols, dpi=150, facecolor='w',
                              figsize=(6*ncols + 2, 13), sharex='col')
    if ncols == 1:
        axes = np.expand_dims(axes, axis=1)

    style_map = {
        (3, False): ('tab:blue',   '--',  '3rd, FF OFF'),
        (3, True):  ('tab:orange', '-.',   '3rd, FF ON'),
        (4, False): ('tab:green',  '--',  '4th, FF OFF'),
        (4, True):  ('tab:red',    '-',   '4th, FF ON'),
    }

    for col_idx, (results_set, load_label) in enumerate(target_results):
        first_load = True
        first = True
        for res in results_set:
            t = res['time']
            key = (res['observer_order'], res['feedforward_on'])
            color, ls, lbl = style_map[key]

            # Speed
            ax = axes[0, col_idx]
            if first:
                ax.plot(t, res['cmd_rpm'], 'k--', lw=0.8, label=r'$\omega^*$')
                first = False
            ax.plot(t, res['act_rpm'], color=color, ls=ls, lw=0.8, label=lbl)
            ax.set_ylabel('Speed [rpm]')
            ax.set_title(load_label, fontsize=11)

            # Speed error
            ax = axes[1, col_idx]
            ax.plot(t, res['speed_err'], color=color, ls=ls, lw=0.7, label=lbl)
            ax.set_ylabel('Speed Err [rpm]')

            # Load estimation
            ax = axes[2, col_idx]
            if first_load:
                ax.plot(t, res['TLoad'], 'k-', lw=1.2, label=r'$T_L$ actual')
                first_load = False
            ax.plot(t, -res['TL_est'], color=color, ls=ls, lw=0.8,
                    label=f'est. ({lbl})')
            ax.set_ylabel('Torque [Nm]')

            # pT estimation (Derivative of disturbance)
            ax = axes[3, col_idx]
            ax.plot(t, -res['pT_est'], color=color, ls=ls, lw=0.7,
                    label=f'$-\\hat{{x}}_3$ ({lbl})')
            ax.set_ylabel(r'$d\hat{T}_L/dt$ [Nm/s]')

    for ax in axes.flat:
        ax.legend(loc='best', fontsize=6, ncol=2)
        ax.grid(True, alpha=0.3)
    for ax in axes[-1, :]:
        ax.set_xlabel('Time [s]')
    for i in range(ncols):
        axes[1, i].axhline(0, color='gray', lw=0.5)

    load_names = " vs ".join([label.split()[0] for _, label in target_results])
    fig.suptitle(f'3rd vs 4th Order ESO: {load_names} Load{cli_suffix}',
                 fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f'  Figure saved: {save_path}')

    return fig

def analyze_sweep_bode(res, t_start=0.05, f_end=200.0):
    """
    基于DFT提取扫频分段数据的基波幅值和相位。
    """
    f_list = []
    current_f = 0.2
    
    t_nodes = [t_start]
    while current_f <= f_end:
        f_list.append(current_f)
        T_cycle = 1.0 / current_f
        t_nodes.append(t_nodes[-1] + T_cycle)
        
        if current_f < 0.3:
            current_f = 0.5
        elif current_f < 0.7:
            current_f = 1.0
        else:
            current_f += 1.0

    freqs = []
    mag_est = []
    phase_est = []
    mag_rej = []
    phase_rej = []
    
    t = res['time']
    TLoad = res['TLoad']
    TL_est = -res['TL_est']
    speed_err = res['speed_err']
    
    for i in range(len(f_list)):
        f = f_list[i]
        t_0 = t_nodes[i]
        t_1 = t_nodes[i+1]
        
        mask = (t >= t_0) & (t < t_1)
        t_seg = t[mask]
        
        if len(t_seg) == 0:
            break
            
        dt = 1e-4
        
        basis = np.exp(-1j * 2 * np.pi * f * (t_seg - t_0))
        
        X_in = np.sum(TLoad[mask] * basis) * dt * f * 2.0
        X_est = np.sum(TL_est[mask] * basis) * dt * f * 2.0
        X_rej = np.sum(speed_err[mask] * basis) * dt * f * 2.0
        
        eps = 1e-12
        amp_in = np.abs(X_in) + eps
        
        H_est = X_est / X_in if amp_in > eps else 0j
        H_rej = X_rej / X_in if amp_in > eps else 0j
        
        freqs.append(f)
        mag_est.append(20 * np.log10(np.abs(H_est) + eps))
        phase_est.append(np.angle(H_est, deg=True))
        mag_rej.append(20 * np.log10(np.abs(H_rej) + eps))
        phase_rej.append(np.angle(H_rej, deg=True))
        
    return {
        'freqs': np.array(freqs),
        'mag_est': np.array(mag_est),
        'phase_est': np.array(phase_est),
        'mag_rej': np.array(mag_rej),
        'phase_rej': np.array(phase_rej)
    }

def plot_bode(sweep_results, omega_ob=100.0, save_path='fig_4th_order_ESO_bode.png'):
    """
    针对 Sweep 数据集绘制 Bode 图，并标注理论带宽线。
    """
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), dpi=150)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    for idx, res in enumerate(sweep_results):
        lbl = res.get('label', f'Index {idx}')
        c = colors[idx % len(colors)]
        
        bode_res = analyze_sweep_bode(res)
        f = bode_res['freqs']
        
        if len(f) == 0:
            continue
            
        axes[0].semilogx(f, bode_res['mag_est'], color=c, lw=1.5, marker='o', markersize=4, label=lbl)
        ph = np.unwrap(bode_res['phase_est'] * np.pi / 180.0) * 180.0 / np.pi
        axes[1].semilogx(f, ph, color=c, lw=1.5, marker='o', markersize=4, label=lbl)
        axes[2].semilogx(f, bode_res['mag_rej'], color=c, lw=1.5, marker='o', markersize=4, label=lbl)

    # 绘制截止频率辅助线
    f_ob = omega_ob / (2 * np.pi)
    for ax in axes:
        ax.axvline(x=f_ob, color='black', linestyle='--', lw=1.2, alpha=0.6, label=f'$f_{{ob}}$={f_ob:.1f}Hz')
        
    axes[0].axhline(y=-3, color='red', linestyle=':', lw=1.0, alpha=0.8, label='-3 dB')

    axes[0].set_title("Disturbance Estimation ($T_{L,est} / T_L$) Magnitude", fontsize=10)
    axes[0].set_ylabel("Gain [dB]")
    axes[0].grid(True, which='both', ls='--', alpha=0.5)
    axes[0].legend(loc='lower left', fontsize=8, ncol=2)
    
    axes[1].set_title("Disturbance Estimation ($T_{L,est} / T_L$) Phase", fontsize=10)
    axes[1].set_ylabel("Phase [deg]")
    axes[1].set_xlabel("Frequency [Hz]")
    axes[1].grid(True, which='both', ls='--', alpha=0.5)
    
    axes[2].set_title("Disturbance Rejection ($Speed_{err} / T_L$) Magnitude", fontsize=10)
    axes[2].set_ylabel("Gain [(rpm)/Nm dB]")
    axes[2].set_xlabel("Frequency [Hz]")
    axes[2].grid(True, which='both', ls='--', alpha=0.5)
    
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches='tight')
    print(f'  Bode figure saved: {save_path}')
    return fig

# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':
    print('=' * 70)
    print(' 四阶 ESO 扰动估计与前馈补偿演示')
    print('=' * 70)

    import argparse
    parser = argparse.ArgumentParser(description="4th-order ESO Simulation")
    parser.add_argument('-t', '--time', type=float, default=1.5,
                        help="Simulation total time in seconds (default: 1.5s)")
    parser.add_argument('--ki', type=float, default=None,
                        help="Speed controller KI coefficient (if provided, overrides auto-tuned KI)")
    parser.add_argument('--scenarios', nargs='+', default=None,
                        help="List of scenarios to run. E.g., '3_OFF' '3_ON' '4_OFF' '4_ON' (default: all)")
    parser.add_argument('--loads', nargs='+', default=None,
                        help="List of load scenarios to run. E.g., 'step' 'ramp' 'sine' 'parabola' 'sweep' (default: all except sweep)")
    parser.add_argument('--omega_ob', type=float, default=100.0,
                        help="Observer bandwidth in rad/s (default: 100.0)")
    parser.add_argument('--cmd_rpm', type=float, default=50.0,
                        help="Constant speed command in rpm (default: 50.0)")
    args = parser.parse_args()

    all_configs = {
        '3_OFF': (3, False),
        '3_ON':  (3, True),
        '4_OFF': (4, False),
        '4_ON':  (4, True)
    }
    
    if args.scenarios is None:
        target_scenarios = ['3_OFF', '3_ON', '4_OFF', '4_ON']
    else:
        target_scenarios = args.scenarios

    d = get_base_d()
    d['NUMBER_OF_SLICES'] = int(args.time / d['TIME_SLICE'])
    d['override_speed_ki'] = args.ki
    
    ki_str = f"KI={args.ki}" if args.ki is not None else "KI=Auto"
    cli_suffix = f" (Time={args.time}s, {ki_str})"
    
    omega_ob_3rd = args.omega_ob  # 3阶观测器带宽 [rad/s]
    omega_ob_4th = args.omega_ob  # 4阶观测器带宽 [rad/s]

    t_total_start = _time.time()

    # ==================================================================
    # Part 1: Run Scenarios
    # ==================================================================
    available_loads = ['step', 'ramp', 'sine', 'parabola', 'sweep']
    if args.loads is None:
        load_types = ['step', 'ramp', 'sine', 'parabola']
    else:
        load_types = [l for l in args.loads if l in available_loads]
        
    all_results = {}

    for ltype in load_types:
        print('\n' + '='*50)
        print(f' Part: {ltype.capitalize()} Load')
        print('='*50)

        results = []
        for sce in target_scenarios:
            if sce not in all_configs:
                print(f"Warning: Unknown scenario '{sce}', skipping.")
                continue
            obs_order, ff_on = all_configs[sce]
            w_ob = omega_ob_3rd if obs_order == 3 else omega_ob_4th
            
            res = run_scenario(d, observer_order=obs_order, feedforward_on=ff_on,
                               load_type=ltype, omega_ob=w_ob, cmd_rpm=args.cmd_rpm)
            results.append(res)
        all_results[ltype] = results

        # Overlay plot for this load type
        plot_overlay(results, title_suffix=f' ({ltype.capitalize()} Load){cli_suffix}',
                     save_path=f'fig_4th_order_ESO_{ltype}_overlay.png')

    # ==================================================================
    # Part 2: 综合对比图
    # ==================================================================
    print('\n' + '='*50)
    print(' Part 2: Combined comparison')
    print('='*50)

    target_results_for_plot = []
    for ltype in load_types:
        if ltype == 'sweep':
            # Cut at t = 9.98 seconds, which is when f finishes 10 Hz cycle and steps to 11 Hz
            t_cut = 9.98 
            res_low = []
            res_high = []
            
            for r in all_results[ltype]:
                t = r['time']
                mask_low = t < t_cut
                mask_high = t >= t_cut
                
                r_low = r.copy()
                r_high = r.copy()
                for k in ['time', 'cmd_rpm', 'act_rpm', 'speed_err', 'TLoad', 'TL_est', 'pT_est']:
                    r_low[k] = r[k][mask_low]
                    r_high[k] = r[k][mask_high]
                
                # Check if we actually have data for these sweeps
                if len(r_low['time']) > 0: res_low.append(r_low)
                if len(r_high['time']) > 0: res_high.append(r_high)
            
            if res_low: target_results_for_plot.append((res_low, 'Sweep (≤10Hz)'))
            if res_high: target_results_for_plot.append((res_high, 'Sweep (>10Hz)'))
        else:
            target_results_for_plot.append((all_results[ltype], f'{ltype.capitalize()} Load'))

    if target_results_for_plot:
        fig_combined = plot_combined(target_results_for_plot,
                                     save_path='fig_4th_order_ESO_comparison.png',
                                     cli_suffix=cli_suffix)

    if 'sweep' in load_types:
        plot_bode(all_results['sweep'], omega_ob=args.omega_ob, save_path='fig_4th_order_ESO_bode.png')

    elapsed = _time.time() - t_total_start
    print(f'\n总耗时: {elapsed:.1f} 秒')

    # ==================================================================
    # 结果总结
    # ==================================================================
    print('\n' + '='*70)
    print(' Results Summary')
    print('='*70)
    print(f'  Observer bandwidth: 3rd-order omega_ob = {omega_ob_3rd} rad/s, 4th-order omega_ob = {omega_ob_4th} rad/s\n')

    for ltype in load_types:
        print(f'  [{ltype.capitalize()} Load]')
        for res in all_results[ltype]:
            half = len(res['speed_err']) // 2
            rms = np.sqrt(np.mean(res['speed_err'][half:]**2))
            max_err = np.max(np.abs(res['speed_err'][half:]))
            TL_est_err = res['TLoad'][half:] - (-res['TL_est'][half:])
            rms_TL = np.sqrt(np.mean(TL_est_err**2))
            print(f'    {res["label"]:30s}  '
                  f'speed_RMS={rms:8.4f} rpm, speed_MAX={max_err:8.4f} rpm, '
                  f'TL_est_RMS={rms_TL:.6f} Nm')
        print()

    print('Output files:')
    for ltype in load_types:
        print(f'  fig_4th_order_ESO_{ltype}_overlay.png')
    print('  fig_4th_order_ESO_comparison.png')

    plt.close('all')
