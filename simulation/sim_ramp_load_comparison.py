# -*- coding: utf-8 -*-
"""
斜坡负载扰动仿真对比 — 4阶ESO前馈补偿
========================================
仿真工况：
  - 电机空载运行至 500 rpm
  - t = 1 s 时突然施加斜坡负载
  - 对比三组不同 (zeta, CLBW) 参数:  无ESO vs 4阶ESO前馈补偿

思路：
  传统PI速度环的内模原理只能抗阶跃扰动（积分器），对斜坡扰动有稳态误差。
  4阶ESO把扰动模型拓展为 d²T_L/dt² = 0（即T_L可以是斜坡），
  多估一个状态 x[3] = dT_L/dt，从而使得 x[2] ≈ -T_L 在斜坡下也无稳态误差。
  将 -x[2] 前馈到 iq* 上，等效于在速度环外加了一个扰动补偿通道。
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
                              ramp_load_max=0.2, total_time=2.5,
                              enable_4th_ESO=False, omega_ob=200):
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
    enable_4th_ESO : bool    是否启用4阶ESO + 前馈补偿
    omega_ob : float         观测器带宽 [rad/s]

    Returns
    -------
    times, speed_rpm, cmd_rpm_arr, tem, tload, iq, tl_est
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
        CTRL.index_separate_speed_estimation = 1   # 启用观测器
        CTRL.use_disturbance_feedforward_rejection = 1  # xS[2] 前馈
        # 4阶位置观测器增益 (特征多项式 = (s + omega_ob)^4)
        CTRL.ell1 = 4 * omega_ob
        CTRL.ell2 = 6 * omega_ob**2
        CTRL.ell3 = 4 * omega_ob**3 * J_s / n_pp
        CTRL.ell4 =     omega_ob**4 * J_s / n_pp
    # else: 默认 index_separate_speed_estimation=0, 不跑观测器, 不前馈

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
    slice_dt = 0.01  # 每 10 ms 更新一次负载
    n_slices = int(total_time / slice_dt)

    all_times  = None
    all_speed  = None
    all_cmd    = None
    all_tem    = None
    all_tload  = None
    all_iq     = None
    all_tl_est = None  # ESO 估计的负载转矩

    def append_global(g, l):
        return l if g is None else np.append(g, l)

    for si in range(n_slices):
        t_now = si * slice_dt

        # 设置速度指令
        CTRL.cmd_rpm = cmd_rpm

        # 设置斜坡负载
        if t_now < ramp_start_time:
            ACM.TLoad = 0.0
        elif t_now < ramp_start_time + ramp_duration:
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
        #  [1]  = ACM.omega_r_mech (rpm)
        #  [5]  = ACM.Tem
        #  [12] = CTRL.cmd_rpm
        #  [17] = CTRL.xS[2]  (ESO 估计的扰动, ≈ -TLoad)
        #  [41] = ACM.TLoad
        #  [9]  = CTRL.idq[1] (iq)
        all_times  = append_global(all_times,  machine_times)
        all_speed  = append_global(all_speed,  watch_data[1])
        all_cmd    = append_global(all_cmd,    watch_data[12])
        all_tem    = append_global(all_tem,    watch_data[5])
        all_tload  = append_global(all_tload,  watch_data[41])
        all_iq     = append_global(all_iq,     watch_data[9])
        all_tl_est = append_global(all_tl_est, watch_data[17])  # xS[2]

    return all_times, all_speed, all_cmd, all_tem, all_tload, all_iq, all_tl_est


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':
    print('=' * 60)
    print(' 斜坡负载扰动: PI only vs PI + 4阶ESO前馈补偿')
    print('=' * 60)

    # ---------- 基础电机参数 ----------
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

    # ---------- 三组控制器参数 ----------
    test_cases = [
        (5,   500,  r'$\zeta$=5,  CLBW=500 Hz',   '#e74c3c', '-'),
        (15, 1000,  r'$\zeta$=15, CLBW=1000 Hz',  '#2ecc71', '--'),
        (25, 2000,  r'$\zeta$=25, CLBW=2000 Hz',  '#3498db', '-.'),
    ]

    # ---------- 仿真参数 ----------
    CMD_RPM       = 500
    RAMP_START    = 1.0
    RAMP_DURATION = 0.5
    RAMP_LOAD_MAX = 5 * 0.4   # 2.0 Nm
    TOTAL_TIME    = 2.5
    OMEGA_OB      = 200       # ESO 带宽 [rad/s]

    # ---------- 运行 6 组仿真 (3 param × 2 ESO modes) ----------
    results_no_eso  = []
    results_with_eso = []

    for zeta, clbw, label, color, ls in test_cases:
        # --- 无 ESO ---
        print(f'\n  [No ESO]  {label} ...')
        out = run_ramp_load_simulation(
            d, zeta, clbw, cmd_rpm=CMD_RPM,
            ramp_start_time=RAMP_START, ramp_duration=RAMP_DURATION,
            ramp_load_max=RAMP_LOAD_MAX, total_time=TOTAL_TIME,
            enable_4th_ESO=False,
        )
        results_no_eso.append({
            'label': label, 'color': color, 'ls': ls,
            'times': out[0], 'speed': out[1], 'cmd': out[2],
            'tem': out[3], 'tload': out[4], 'iq': out[5], 'tl_est': out[6],
        })

        # --- 4阶 ESO + 前馈 ---
        print(f'  [4th ESO] {label} ...')
        out = run_ramp_load_simulation(
            d, zeta, clbw, cmd_rpm=CMD_RPM,
            ramp_start_time=RAMP_START, ramp_duration=RAMP_DURATION,
            ramp_load_max=RAMP_LOAD_MAX, total_time=TOTAL_TIME,
            enable_4th_ESO=True, omega_ob=OMEGA_OB,
        )
        results_with_eso.append({
            'label': label, 'color': color, 'ls': ls,
            'times': out[0], 'speed': out[1], 'cmd': out[2],
            'tem': out[3], 'tload': out[4], 'iq': out[5], 'tl_est': out[6],
        })

    # ===================================================================
    # 绘图：左列 = 无ESO, 右列 = 4阶ESO前馈
    # ===================================================================
    print('\n  正在绘图...')
    plt.style.use('bmh')
    mpl.rc('font', family='Times New Roman', size=10.0)
    mpl.rc('legend', fontsize=8)
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.unicode_minus'] = False

    fig, axes = plt.subplots(3, 2, dpi=150, facecolor='w', figsize=(16, 13),
                              sharex='col', sharey='row')

    fig.suptitle(r'Ramp Load Rejection: PI Only  vs  PI + 4th-Order ESO Feedforward'
                 f'\n(500 rpm, ramp load 0→{RAMP_LOAD_MAX:.1f} Nm at t=1 s, '
                 r'ESO $\omega_{ob}$' + f'={OMEGA_OB} rad/s)',
                 fontsize=13, fontweight='bold', y=0.99, fontfamily='Times New Roman')

    col_titles = ['PI Speed Loop Only (No ESO)',
                  'PI + 4th-Order ESO Feedforward']

    for col_idx, (results_set, col_title) in enumerate([
        (results_no_eso, col_titles[0]),
        (results_with_eso, col_titles[1]),
    ]):
        r0 = results_set[0]

        # ---- Row 0: Speed ----
        ax = axes[0, col_idx]
        ax.plot(r0['times'], r0['cmd'], color='gray', lw=1.0, alpha=0.5, label='Speed command')
        for r in results_set:
            ax.plot(r['times'], r['speed'], color=r['color'], ls=r['ls'], lw=1.5, label=r['label'])
        ax.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.4)
        ax.axvspan(RAMP_START, RAMP_START + RAMP_DURATION, alpha=0.05, color='red')
        ax.set_ylabel('Speed [rpm]')
        ax.set_title(col_title, fontweight='bold', fontsize=11)
        ax.legend(loc='lower right', fontsize=6.5, ncol=1, framealpha=0.9)
        ax.grid(True, alpha=0.3)

        # Inset zoom on disturbance region
        axins = ax.inset_axes([0.12, 0.08, 0.38, 0.42])
        for r in results_set:
            axins.plot(r['times'], r['speed'], color=r['color'], ls=r['ls'], lw=1.5)
        axins.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.4)
        axins.set_xlim(0.9, 2.2)
        speed_vals_region = []
        for r in results_set:
            mask = (r['times'] > 0.9) & (r['times'] < 2.2)
            speed_vals_region.extend(r['speed'][mask].tolist())
        if speed_vals_region:
            sp_min, sp_max = min(speed_vals_region), max(speed_vals_region)
            margin = max((sp_max - sp_min) * 0.15, 5.0)
            axins.set_ylim(sp_min - margin, sp_max + margin)
        axins.tick_params(labelsize=6)
        axins.grid(True, alpha=0.3)
        axins.set_title('Zoom', fontsize=6, color='gray')
        ax.indicate_inset_zoom(axins, edgecolor='gray', alpha=0.4)

        # ---- Row 1: Speed error ----
        ax = axes[1, col_idx]
        for r in results_set:
            err = r['cmd'] - r['speed']
            ax.plot(r['times'], err, color=r['color'], ls=r['ls'], lw=1.5, label=r['label'])
        ax.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.4)
        ax.axvspan(RAMP_START, RAMP_START + RAMP_DURATION, alpha=0.05, color='red')
        ax.axhline(y=0, color='gray', ls='-', lw=0.5)
        ax.set_ylabel('Speed Error [rpm]')
        ax.legend(loc='best', fontsize=6.5, ncol=1, framealpha=0.9)
        ax.grid(True, alpha=0.3)

        # ---- Row 2: Load torque & ESO estimation ----
        ax = axes[2, col_idx]
        ax.fill_between(r0['times'], 0, r0['tload'], alpha=0.20, color='orange',
                         label=f'Actual load torque')
        for r in results_set:
            # xS[2] ≈ -TLoad (观测器约定), 所以画 -xS[2]
            ax.plot(r['times'], -r['tl_est'], color=r['color'], ls=r['ls'], lw=1.2,
                    label=r'$-\hat{x}_2$ (' + r['label'] + ')')
        ax.axvline(x=RAMP_START, color='k', ls=':', lw=0.8, alpha=0.4)
        ax.axvspan(RAMP_START, RAMP_START + RAMP_DURATION, alpha=0.05, color='red')
        ax.set_ylabel('Torque [Nm]')
        ax.set_xlabel('Time [s]')
        ax.legend(loc='best', fontsize=6, ncol=1, framealpha=0.9)
        ax.grid(True, alpha=0.3)

    # 行标题
    row_labels = ['(a) Speed Response', '(b) Speed Tracking Error',
                  '(c) Load Torque Estimation']
    for i, lbl in enumerate(row_labels):
        axes[i, 0].annotate(lbl, xy=(0.02, 0.95), xycoords='axes fraction',
                            fontsize=9, fontweight='bold', va='top',
                            bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8))

    fig.tight_layout(rect=[0, 0, 1, 0.94])

    save_path = 'sim_ramp_load_comparison.png'
    fig.savefig(save_path, dpi=200, bbox_inches='tight')
    print(f'\n  图已保存: {save_path}')
    plt.show()
    plt.close('all')
    print('  完成!')
