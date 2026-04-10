# -*- coding: utf-8 -*-
"""
斜坡负载扰动仿真对比 (无ESO版本, 仅PI速度环)
=============================================
仿真工况：
  - 电机空载运行至 500 rpm
  - t = 1 s 时突然施加斜坡负载（线性增长至额定负载）
  - 对比三组不同 (zeta, CLBW) 参数下的速度与转矩响应

本脚本展示了纯 PI 速度环在面对斜坡负载扰动时的局限性：
  - PI 积分器只能抗阶跃扰动
  - 对斜坡型扰动存在稳态跟踪误差
  - 不同的 (zeta, CLBW) 参数组合影响瞬态响应和稳态误差
"""

from pylab import np, plt, mpl
import copy

from tuner import (
    get_coeffs_dc_motor_current_regulator,
    get_coeffs_dc_motor_SPEED_regulator,
)
from tutorials_ep6_svpwm import (
    The_Motor_Controller, The_AC_Machine,
    The_PID_Regulator, ACMSimPyIncremental,
)


def run_ramp_load_simulation(d, zeta, CLBW_Hz, cmd_rpm=500,
                              ramp_start_time=1.0, ramp_duration=0.5,
                              ramp_load_max=0.2, total_time=2.5):
    """
    运行一次完整仿真：空载加速到 cmd_rpm，然后在 ramp_start_time 施加斜坡负载。

    Parameters
    ----------
    d : dict             基础电机参数
    zeta : float         阻尼系数 (FOC_delta)
    CLBW_Hz : float      电流环带宽 [Hz]
    cmd_rpm : float      目标转速 [rpm]
    ramp_start_time : float  斜坡负载开始时间 [s]
    ramp_duration : float    斜坡负载上升时间 [s]
    ramp_load_max : float    最终负载转矩 [Nm]
    total_time : float       仿真总时间 [s]

    Returns
    -------
    times : ndarray      时间数组
    speed_rpm : ndarray  实际转速 [rpm]
    cmd_rpm_arr : ndarray 指令转速 [rpm]
    tem : ndarray        电磁转矩 [Nm]
    tload : ndarray      负载转矩 [Nm]
    iq : ndarray         q轴电流 [A]
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
    CTRL.bool_overwrite_speed_commands = True  # 我们手动控制速度指令
    CTRL.bool_zero_id_control = dd['CTRL.bool_zero_id_control']
    CTRL.bool_apply_speed_closed_loop_control = True

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=dd['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])

    # PI 调节器
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

    # === 分段仿真 ===
    MACHINE_TS = CL_TS / ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    # 小切片时间步长（与负载施加精度相关）
    slice_dt = 0.01  # 每 10 ms 更新一次负载
    n_slices = int(total_time / slice_dt)

    all_times = None
    all_speed = None
    all_cmd   = None
    all_tem   = None
    all_tload = None
    all_iq    = None

    def append_global(g, l):
        return l if g is None else np.append(g, l)

    current_time = 0.0
    for si in range(n_slices):
        t_now = si * slice_dt

        # 设置速度指令
        CTRL.cmd_rpm = cmd_rpm

        # 设置斜坡负载
        if t_now < ramp_start_time:
            ACM.TLoad = 0.0
        elif t_now < ramp_start_time + ramp_duration:
            # 线性斜坡
            progress = (t_now - ramp_start_time) / ramp_duration
            ACM.TLoad = ramp_load_max * progress
        else:
            ACM.TLoad = ramp_load_max

        # 运行一个切片
        machine_times, watch_data = ACMSimPyIncremental(
            t0=t_now, TIME=slice_dt,
            ACM=ACM, CTRL=CTRL,
            reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed
        )

        # watch_data 索引参考：
        #  [1]  = ACM.omega_r_mech (rpm, 已转换)
        #  [5]  = ACM.Tem
        #  [12] = CTRL.cmd_rpm
        #  [41] = ACM.TLoad
        #  [9]  = CTRL.idq[1] (iq)
        all_times = append_global(all_times, machine_times)
        all_speed = append_global(all_speed, watch_data[1])    # 实际转速 [rpm]
        all_cmd   = append_global(all_cmd,   watch_data[12])   # 指令转速
        all_tem   = append_global(all_tem,   watch_data[5])    # 电磁转矩
        all_tload = append_global(all_tload, watch_data[41])   # 负载转矩
        all_iq    = append_global(all_iq,    watch_data[9])    # iq 电流

    return all_times, all_speed, all_cmd, all_tem, all_tload, all_iq


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':
    print('=' * 60)
    print(' 斜坡负载扰动仿真对比 (不同 zeta & CLBW)')
    print('=' * 60)

    # ---------- 基础电机参数（小电感电机，与 ep6 一致） ----------
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
        'DC_BUS_VOLTAGE': 48,  # 提高母线电压以支持 500 rpm 运行 (back-EMF ≈ 14.4 V)
        'CTRL.bool_apply_speed_closed_loop_control': True,
        'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
        'CTRL.bool_apply_sweeping_frequency_excitation': False,
        'CTRL.bool_overwrite_speed_commands': True,
        'CTRL.bool_zero_id_control': True,
        'FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False': 10,
        'VL_LIMIT_OVERLOAD_FACTOR': 3.0,
    }

    # ---------- 三组参数 ----------
    test_cases = [
        # (zeta, CLBW_Hz, label, color, linestyle)
        (5,   500,  r'$\zeta$=5,  CLBW=500 Hz  (Low damping, Low BW)',   '#e74c3c', '-'),
        (15, 1000,  r'$\zeta$=15, CLBW=1000 Hz (Mid damping, Mid BW)',  '#2ecc71', '--'),
        (25, 2000,  r'$\zeta$=25, CLBW=2000 Hz (High damping, High BW)',  '#3498db', '-.'),
    ]

    # ---------- 仿真参数 ----------
    CMD_RPM       = 500     # 目标转速
    RAMP_START    = 1.0     # 斜坡负载开始时间
    RAMP_DURATION = 0.5     # 斜坡上升持续时间
    RAMP_LOAD_MAX = 5*0.4   # 最终负载 [Nm] (增大以凸显差异)
    TOTAL_TIME    = 2.5     # 仿真总时间 [s]

    # ---------- 运行仿真 ----------
    results = []
    for zeta, clbw, label, color, ls in test_cases:
        print(f'\n  正在仿真: {label} ...')
        times, speed, cmd, tem, tload, iq = run_ramp_load_simulation(
            d, zeta, clbw,
            cmd_rpm=CMD_RPM,
            ramp_start_time=RAMP_START,
            ramp_duration=RAMP_DURATION,
            ramp_load_max=RAMP_LOAD_MAX,
            total_time=TOTAL_TIME,
        )
        results.append({
            'zeta': zeta, 'clbw': clbw,
            'label': label, 'color': color, 'ls': ls,
            'times': times, 'speed': speed, 'cmd': cmd,
            'tem': tem, 'tload': tload, 'iq': iq,
        })
        print(f'    完成! 数据点数: {len(times)}')

    # ---------- 绘图 ----------
    print('\n  正在绘图...')
    plt.style.use('bmh')
    mpl.rc('font', family='Times New Roman', size=11.0)
    mpl.rc('legend', fontsize=9)
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.unicode_minus'] = False

    fig, axes = plt.subplots(4, 1, dpi=150, facecolor='w', figsize=(13, 16), sharex=True)
    fig.suptitle('PMSM Ramp Load Disturbance Response Comparison\n'
                 r'(No-load $\rightarrow$ 500 rpm, ramp load at $t$ = 1 s)',
                 fontsize=14, fontweight='bold', y=0.98, fontfamily='Times New Roman')

    # --- 子图1: 转速响应 (全局 + 放大 inset) ---
    ax = axes[0]
    r0 = results[0]
    ax.plot(r0['times'], r0['cmd'], color='gray', lw=1.0, alpha=0.6, label='Speed command')
    for r in results:
        ax.plot(r['times'], r['speed'], color=r['color'], ls=r['ls'], lw=1.5, label=r['label'])
    ax.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.5)
    ax.axvspan(RAMP_START, RAMP_START + RAMP_DURATION, alpha=0.06, color='red')
    ax.annotate('Ramp load', xy=(RAMP_START + 0.02, CMD_RPM * 0.55),
                fontsize=8, color='dimgray', fontweight='bold')
    ax.set_ylabel('Speed [rpm]')
    ax.set_title('(a) Speed Response', fontweight='bold', loc='left')
    ax.legend(loc='lower right', fontsize=7, ncol=2, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # Inset: zoom in on load disturbance region
    axins1 = ax.inset_axes([0.15, 0.08, 0.35, 0.40])  # [x, y, width, height] in axes coords
    for r in results:
        axins1.plot(r['times'], r['speed'], color=r['color'], ls=r['ls'], lw=1.5)
    axins1.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.5)
    axins1.set_xlim(0.9, 2.0)
    # auto-determine y range around disturbance
    speed_min = min(np.min(r['speed'][(r['times'] > 0.9) & (r['times'] < 2.0)]) for r in results)
    speed_max = max(np.max(r['speed'][(r['times'] > 0.9) & (r['times'] < 2.0)]) for r in results)
    margin = max((speed_max - speed_min) * 0.3, 2.0)
    axins1.set_ylim(speed_min - margin, speed_max + margin)
    axins1.set_title('Zoom in', fontsize=7, color='gray')
    axins1.tick_params(labelsize=7)
    axins1.grid(True, alpha=0.3)
    ax.indicate_inset_zoom(axins1, edgecolor='gray', alpha=0.5)

    # --- 子图2: 转速误差 (仅显示负载施加后的区域) ---
    ax = axes[1]
    for r in results:
        speed_error = r['cmd'] - r['speed']
        ax.plot(r['times'], speed_error, color=r['color'], ls=r['ls'], lw=1.5, label=r['label'])
    ax.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.5)
    ax.axvspan(RAMP_START, RAMP_START + RAMP_DURATION, alpha=0.06, color='red')
    ax.axhline(y=0, color='gray', ls='-', lw=0.5)
    ax.set_ylabel('Speed Error [rpm]')
    ax.set_title('(b) Speed Tracking Error', fontweight='bold', loc='left')
    ax.legend(loc='best', fontsize=7, ncol=1, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # Inset: focus on disturbance region
    axins2 = ax.inset_axes([0.55, 0.45, 0.40, 0.50])
    for r in results:
        speed_error = r['cmd'] - r['speed']
        axins2.plot(r['times'], speed_error, color=r['color'], ls=r['ls'], lw=1.5)
    axins2.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.5)
    axins2.axhline(y=0, color='gray', ls='-', lw=0.5)
    axins2.set_xlim(0.9, 2.0)
    err_vals = []
    for r in results:
        mask = (r['times'] > 0.9) & (r['times'] < 2.0)
        err_vals.extend((r['cmd'][mask] - r['speed'][mask]).tolist())
    err_min, err_max = min(err_vals), max(err_vals)
    err_margin = max((err_max - err_min) * 0.2, 0.5)
    axins2.set_ylim(err_min - err_margin, err_max + err_margin)
    axins2.set_title('Zoom in (disturbance region)', fontsize=7, color='gray')
    axins2.tick_params(labelsize=7)
    axins2.grid(True, alpha=0.3)
    ax.indicate_inset_zoom(axins2, edgecolor='gray', alpha=0.5)

    # --- 子图3: 电磁转矩 & 负载转矩 ---
    ax = axes[2]
    ax.fill_between(r0['times'], 0, r0['tload'], alpha=0.20, color='orange',
                     label=f'Load torque (ramp to {RAMP_LOAD_MAX} Nm)')
    for r in results:
        short_label = f'$T_{{em}}$ ({r["label"][:r["label"].find("(")].strip()})'
        ax.plot(r['times'], r['tem'], color=r['color'], ls=r['ls'], lw=1.5, label=short_label)
    ax.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.5)
    ax.axvspan(RAMP_START, RAMP_START + RAMP_DURATION, alpha=0.06, color='red')
    ax.set_ylabel('Torque [Nm]')
    ax.set_title('(c) Electromagnetic Torque vs Load Torque', fontweight='bold', loc='left')
    ax.legend(loc='best', fontsize=7, ncol=2, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # --- 子图4: q轴电流 ---
    ax = axes[3]
    for r in results:
        ax.plot(r['times'], r['iq'], color=r['color'], ls=r['ls'], lw=1.5, label=r['label'])
    ax.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.5)
    ax.axvspan(RAMP_START, RAMP_START + RAMP_DURATION, alpha=0.06, color='red')
    ax.set_ylabel('$i_q$ [A]')
    ax.set_title('(d) $q$-axis Current Response', fontweight='bold', loc='left')
    ax.legend(loc='best', fontsize=7, ncol=1, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('Time [s]')

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    save_path = 'sim_ramp_load_comparison_no_eso.png'
    fig.savefig(save_path, dpi=200, bbox_inches='tight')
    print(f'\n  图已保存: {save_path}')
    plt.show()
    plt.close('all')
    print('  完成!')
