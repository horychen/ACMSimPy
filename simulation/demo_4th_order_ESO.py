# -*- coding: utf-8 -*-
"""
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
        'NUMBER_OF_SLICES': 300,  # 3 秒仿真 (0.01s × 300)
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
def run_scenario(d_base, observer_order=3, feedforward_on=False,
                 load_type='ramp', omega_ob=100, verbose=True):
    """
    运行一个仿真场景, 返回时间序列和关键波形数据。

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
    if load_type == 'ramp':
        # 线性增加的负载转矩: 从 slice 5 (0.5s) 开始, 斜率 0.2 Nm/s
        dd['user_system_input_code'] = (
            f"CTRL.cmd_rpm = 50\n"
            f"if ii >= 5: ACM.TLoad = 0.2 * (ii - 5) * {TIME_SLICE}"
        )
    elif load_type == 'step':
        # 阶跃负载: 在 slice 10 (1.0s) 时突加 0.15 Nm
        dd['user_system_input_code'] = (
            "CTRL.cmd_rpm = 50\n"
            "if ii >= 10: ACM.TLoad = 0.15"
        )
    elif load_type == 'sine':
        # 正弦负载: 2Hz, 幅值 0.15 Nm
        dd['user_system_input_code'] = (
            "CTRL.cmd_rpm = 50\n"
            f"if ii >= 5: ACM.TLoad = 0.15 * np.sin(2 * np.pi * 2.0 * (ii - 5) * {TIME_SLICE})"
        )
    elif load_type == 'parabola':
        # 抛物线负载: 从 slice 5 (0.5s) 开始, T = 0.1 * t^2
        dd['user_system_input_code'] = (
            "CTRL.cmd_rpm = 50\n"
            f"if ii >= 5: ACM.TLoad = 0.1 * ((ii - 5) * {TIME_SLICE})**2"
        )
    else:
        raise ValueError(f"Unknown load_type: {load_type}")

    # --- 先自动整定 PI 参数 ---
    tuner.tunner_wrapper(dd)

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
        if key == (3, False):  # 只画一次真实负载
            ax.plot(t, res['TLoad'], 'k-', lw=1.2, label=r'$T_L$ actual')
        ax.plot(t, -res['TL_est'], color=color, ls=ls, lw=0.8,
                label=f'$-\\hat{{x}}_2$ ({lbl})')
        ax.set_ylabel('Load Torque [Nm]')

        # Row 3: pT estimation
        ax = axes[3]
        ax.plot(t, res['pT_est'], color=color, ls=ls, lw=0.7,
                label=f'$\\hat{{x}}_3$ ({lbl})')
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


def plot_combined(ramp_results, step_results, save_path=None):
    """
    综合对比图: Ramp + Step 并排显示。
    """
    plt.style.use('bmh')
    mpl.rc('font', family='Times New Roman', size=10.0)
    mpl.rc('legend', fontsize=7)
    mpl.rcParams['lines.linewidth'] = 1.0
    mpl.rcParams['mathtext.fontset'] = 'stix'

    fig, axes = plt.subplots(3, 2, dpi=150, facecolor='w',
                              figsize=(14, 10), sharex='col')

    style_map = {
        (3, False): ('tab:blue',   '--',  '3rd, FF OFF'),
        (3, True):  ('tab:orange', '-.',   '3rd, FF ON'),
        (4, False): ('tab:green',  '--',  '4th, FF OFF'),
        (4, True):  ('tab:red',    '-',   '4th, FF ON'),
    }

    for col_idx, (results_set, load_label) in enumerate([
        (ramp_results, 'Ramp Load'),
        (step_results, 'Step Load'),
    ]):
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
            if key == (3, False):
                ax.plot(t, res['TLoad'], 'k-', lw=1.2, label=r'$T_L$ actual')
            ax.plot(t, -res['TL_est'], color=color, ls=ls, lw=0.8,
                    label=f'est. ({lbl})')
            ax.set_ylabel('Torque [Nm]')

    for ax in axes.flat:
        ax.legend(loc='best', fontsize=6, ncol=2)
        ax.grid(True, alpha=0.3)
    for ax in axes[-1, :]:
        ax.set_xlabel('Time [s]')
    axes[1, 0].axhline(0, color='gray', lw=0.5)
    axes[1, 1].axhline(0, color='gray', lw=0.5)

    fig.suptitle('3rd vs 4th Order ESO: Ramp Load vs Step Load',
                 fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches='tight')
        print(f'  Figure saved: {save_path}')

    return fig


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':
    print('=' * 70)
    print(' 四阶 ESO 扰动估计与前馈补偿演示')
    print('=' * 70)

    d = get_base_d()
    omega_ob_3rd = 100  # 3阶观测器带宽 [rad/s]
    omega_ob_4th = 100  # 4阶观测器带宽 [rad/s] (ell4 = 100^4 * Js/npp = 200)

    t_total_start = _time.time()

    # ==================================================================
    # Part 1: Run all 4 Scenarios (Step, Ramp, Sine, Parabola)
    # ==================================================================
    load_types = ['step', 'ramp', 'sine', 'parabola']
    all_results = {}

    for ltype in load_types:
        print('\n' + '='*50)
        print(f' Part: {ltype.capitalize()} Load')
        print('='*50)

        results = []
        for obs_order in [3, 4]:
            w_ob = omega_ob_3rd if obs_order == 3 else omega_ob_4th
            for ff_on in [False, True]:
                res = run_scenario(d, observer_order=obs_order, feedforward_on=ff_on,
                                   load_type=ltype, omega_ob=w_ob)
                results.append(res)
        all_results[ltype] = results

        # Overlay plot for this load type
        plot_overlay(results, title_suffix=f'({ltype.capitalize()} Load)',
                     save_path=f'fig_4th_order_ESO_{ltype}_overlay.png')

    # ==================================================================
    # Part 2: 综合对比图 (Ramp vs Step for main paper/report)
    # ==================================================================
    print('\n' + '='*50)
    print(' Part 2: Combined comparison (Ramp vs Step)')
    print('='*50)

    fig_combined = plot_combined(all_results['ramp'], all_results['step'],
                                 save_path='fig_4th_order_ESO_comparison.png')

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
