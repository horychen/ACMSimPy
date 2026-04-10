# -*- coding: utf-8 -*-
"""
ESO 观测器极点配置图 (Pole-Zero Map)
=====================================
展示不同阶数和不同带宽下的观测器误差动态系统极点，
以及三组控制器参数对应的速度环闭环极点。

观测器误差动态系统 (以角度估计误差 e_theta = theta - x[0] 为输出误差):
  de_theta/dt = -ell1 * e_theta + e_omega
  de_omega/dt = -ell2 * e_theta + (npp/Js) * e_d
  de_d/dt     = -ell3 * e_theta + e_p
  de_p/dt     = -ell4 * e_theta

其中：
  e_omega = omega - x[1],  e_d = -(T_L + x[2]),  e_p = -(pT_L + x[3])

采用二项式极点配置法：
  3阶: (s + omega_ob)^3, 增益 ell = [3, 3, 1] * omega_ob^[1,2,3]
  4阶: (s + omega_ob)^4, 增益 ell = [4, 6, 4, 1] * omega_ob^[1,2,3,4]
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from tuner import (
    get_coeffs_dc_motor_current_regulator,
    get_coeffs_dc_motor_SPEED_regulator,
)


# ======================================================================
# Motor parameters (与仿真脚本一致)
# ======================================================================
n_pp = 22
J_s = 0.44e-4
KE = 0.0125
KA = KE
R = 0.035
L = 0.036e-3


def observer_A_matrix(order, omega_ob, n_pp=n_pp, J_s=J_s):
    """
    构建观测器误差动态系统的 A 矩阵。

    3阶 ESO (状态: e_theta, e_omega, e_d):
      扰动模型: dT_L/dt = 0 (阶跃扰动)

    4阶 ESO (状态: e_theta, e_omega, e_d, e_p):
      扰动模型: d²T_L/dt² = 0 (斜坡扰动)
    """
    ratio = n_pp / J_s

    if order == 3:
        ell1 = 3 * omega_ob
        ell2 = 3 * omega_ob**2
        ell3 = omega_ob**3 * J_s / n_pp

        A = np.array([
            [-ell1,    1,       0     ],
            [-ell2,    0,       ratio ],
            [-ell3,    0,       0     ],
        ])
    elif order == 4:
        ell1 = 4 * omega_ob
        ell2 = 6 * omega_ob**2
        ell3 = 4 * omega_ob**3 * J_s / n_pp
        ell4 = omega_ob**4 * J_s / n_pp

        A = np.array([
            [-ell1,    1,       0,       0],
            [-ell2,    0,       ratio,   0],
            [-ell3,    0,       0,       1],
            [-ell4,    0,       0,       0],
        ])
    else:
        raise ValueError(f"order must be 3 or 4, got {order}")

    return A


def speed_loop_closed_loop_poles(zeta, CLBW_Hz):
    """
    计算速度环闭环极点 (简化模型)。

    速度环开环传递函数:
      L(s) = C(s) * G_cl(s) * G_plant(s)

    其中:
      G_plant(s)  = 1.5*npp*KA / (Js*s)        (从 iq 到 omega_elec)
      G_cl(s)     = omega_cl / (s + omega_cl)   (电流环近似)
      C(s)        = Kp_speed * (s + Ki) / s     (PI 控制器)

    闭环特征多项式: s^3 + a2*s^2 + a1*s + a0 = 0
    """
    currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R, L, CLBW_Hz)
    omega_cl = currentKp / L  # 电流环带宽 [rad/s]

    speedKp, speedKi = get_coeffs_dc_motor_SPEED_regulator(
        J_s, n_pp, KA, zeta, omega_cl
    )

    # 开环: Kp_speed * (s + Ki) / s * omega_cl / (s + omega_cl) * 1.5*npp*KA / (Js*s)
    # = Kp_speed * omega_cl * 1.5*npp*KA / Js * (s + Ki) / (s^2 * (s + omega_cl))
    #
    # 闭环特征方程: s^2 * (s + omega_cl) + K * (s + Ki) = 0
    # 其中 K = Kp_speed * omega_cl * 1.5*npp*KA / Js
    K = speedKp * omega_cl * 1.5 * n_pp * KA / J_s

    # s^3 + omega_cl*s^2 + K*s + K*speedKi = 0
    coeffs = [1.0, omega_cl, K, K * speedKi]
    poles = np.roots(coeffs)
    return poles, omega_cl, speedKp, speedKi


# ======================================================================
# MAIN
# ======================================================================
if __name__ == '__main__':

    # --- 三组控制器参数 ---
    test_cases = [
        (5,   500,  r'$\zeta$=5,  CLBW=500 Hz',    '#e74c3c', 'o'),
        (15, 1000,  r'$\zeta$=15, CLBW=1000 Hz',   '#2ecc71', 's'),
        (25, 2000,  r'$\zeta$=25, CLBW=2000 Hz',   '#3498db', 'D'),
    ]

    omega_ob_values = [100, 200, 500]  # 三种观测器带宽

    # ===========================================================
    # Figure 1: 观测器极点 — 3阶 vs 4阶, 不同 omega_ob
    # ===========================================================
    plt.style.use('bmh')
    mpl.rc('font', family='Times New Roman', size=10.0)
    mpl.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.unicode_minus'] = False

    fig, axes = plt.subplots(1, 2, dpi=150, facecolor='w', figsize=(14, 6))

    colors_ob = ['#e74c3c', '#2ecc71', '#3498db']
    markers_ob = ['o', 's', 'D']
    marker_sizes = [160, 120, 80]

    for ax_idx, order in enumerate([3, 4]):
        ax = axes[ax_idx]
        title = f'{order}rd-Order ESO' if order == 3 else f'{order}th-Order ESO'
        ax.set_title(title + ' Pole Placement', fontsize=13, fontweight='bold')

        for i, omega_ob in enumerate(omega_ob_values):
            A = observer_A_matrix(order, omega_ob)
            eigvals = np.linalg.eigvals(A)

            # 极点在实轴上 (imag ≈ 0)
            ax.plot(eigvals[0].real, 0, markers_ob[i],
                    color=colors_ob[i], markersize=np.sqrt(marker_sizes[i]),
                    markeredgecolor='k', markeredgewidth=0.5, zorder=5 + i,
                    label=rf'$\omega_{{ob}}$={omega_ob}: poles at $s$={eigvals[0].real:.0f} (×{order})')

        # 画实轴和虚轴
        ax.axhline(y=0, color='k', lw=1.2, alpha=0.7)
        ax.axvline(x=0, color='k', lw=1.0, alpha=0.5, ls='--')

        # 稳定/不稳定域着色
        ax.axvspan(-700, 0, alpha=0.05, color='green')
        ax.axvspan(0, 150, alpha=0.08, color='red')

        ax.set_xlabel(r'Real Axis $\sigma$ [rad/s]', fontsize=11)
        ax.set_ylabel(r'Imaginary Axis $j\omega$ [rad/s]', fontsize=11)
        ax.legend(loc='upper left', fontsize=8.5, framealpha=0.9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-650, 120)
        ax.set_ylim(-100, 100)

        # 注意：极点全在实轴上，所以 y 方向只显示一个窄范围来强调这一点
        # 但保留足够空间让标注不飞出去

    fig.suptitle('ESO Observer Error Dynamics — Pole Placement Map\n'
                 r'Binomial design $(s + \omega_{ob})^n$:  all poles are repeated real poles at $s = -\omega_{ob}$',
                 fontsize=12, fontweight='bold', fontfamily='Times New Roman')
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig('fig_ESO_pzmap_observer.png', dpi=200, bbox_inches='tight')
    print('  Saved: fig_ESO_pzmap_observer.png')

    # ===========================================================
    # Figure 2: 综合极点图 — 观测器 + 速度环闭环极点
    # ===========================================================
    fig2, ax2 = plt.subplots(1, 1, dpi=150, facecolor='w', figsize=(14, 7))

    ax2.set_title('Combined Pole Map: Speed Loop + 4th-Order ESO\n'
                  r'(ESO $\omega_{ob}$ = 200 rad/s, binomial pole placement)',
                  fontsize=13, fontweight='bold', fontfamily='Times New Roman')

    omega_ob_fixed = 200

    # --- 画观测器极点 ---
    A4 = observer_A_matrix(4, omega_ob_fixed)
    obs_poles = np.linalg.eigvals(A4)
    ax2.scatter(obs_poles.real, obs_poles.imag, s=250, c='#FFD700', marker='*',
                edgecolors='k', linewidths=0.5, zorder=6,
                label=rf'4th ESO poles ($\omega_{{ob}}$={omega_ob_fixed}): all at $s$={obs_poles[0].real:.0f}')
    ax2.annotate(f'ESO: s = {obs_poles[0].real:.0f}\n(4 repeated poles)',
                 xy=(obs_poles[0].real, 0),
                 xytext=(-50, 60), textcoords='offset points',
                 fontsize=9, color='#996600', fontweight='bold',
                 arrowprops=dict(arrowstyle='->', color='#996600', lw=1.5),
                 bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', ec='#996600', alpha=0.9))

    # --- 画速度环极点 ---
    # 每组极点的标注偏移（手动调整避免重叠）
    annot_offsets = [
        [(20, 25), (25, 20)],   # zeta=5:  fast pole, dominant pair
        [(20, 25), (25, -30)],  # zeta=15
        [(20, 25), (25, -30)],  # zeta=25
    ]

    for idx, (zeta, clbw, label, color, marker) in enumerate(test_cases):
        poles, omega_cl, Kp_sp, Ki_sp = speed_loop_closed_loop_poles(zeta, clbw)

        ax2.scatter(poles.real, poles.imag, s=150, c=color, marker=marker,
                    edgecolors='k', linewidths=0.5, zorder=5,
                    label=f'Speed loop: {label}')

        # 标注极点值
        annot_idx = 0
        for j, p in enumerate(poles):
            if abs(p.imag) < 1:
                txt = f's = {p.real:.0f}'
            else:
                txt = f'{p.real:.0f}±j{abs(p.imag):.0f}'
                if p.imag < 0:
                    continue  # 只标注正虚部
            ofs = annot_offsets[idx][min(annot_idx, len(annot_offsets[idx])-1)]
            ax2.annotate(txt, xy=(p.real, p.imag),
                         xytext=ofs, textcoords='offset points',
                         fontsize=7, color=color,
                         arrowprops=dict(arrowstyle='->', color=color, lw=0.6),
                         bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=color, alpha=0.7))
            annot_idx += 1

    # 画坐标轴
    ax2.axhline(y=0, color='k', lw=0.8, alpha=0.4)
    ax2.axvline(x=0, color='k', lw=0.8, alpha=0.5, ls='--')

    # 稳定域着色
    ax2.axvspan(-15000, 0, alpha=0.03, color='green')
    ax2.axvspan(0, 1000, alpha=0.03, color='red')

    ax2.set_xlabel(r'Real Part $\sigma$ [rad/s]', fontsize=11)
    ax2.set_ylabel(r'Imaginary Part $j\omega$ [rad/s]', fontsize=11)
    ax2.legend(loc='upper left', fontsize=8.5, framealpha=0.9, ncol=1)
    ax2.grid(True, alpha=0.3)

    # 自动适配范围
    all_reals = list(obs_poles.real)
    all_imags = list(obs_poles.imag)
    for zeta, clbw, *_ in test_cases:
        poles, *_ = speed_loop_closed_loop_poles(zeta, clbw)
        all_reals.extend(poles.real.tolist())
        all_imags.extend(poles.imag.tolist())
    r_min, r_max = min(all_reals), max(all_reals)
    i_min, i_max = min(all_imags), max(all_imags)
    r_margin = max(abs(r_max - r_min) * 0.15, 500)
    i_margin = max(abs(i_max - i_min) * 0.4, 100)
    ax2.set_xlim(r_min - r_margin, max(r_max + r_margin, 300))
    ax2.set_ylim(i_min - i_margin, i_max + i_margin)

    fig2.tight_layout()
    fig2.savefig('fig_ESO_pzmap_combined.png', dpi=200, bbox_inches='tight')
    print('  Saved: fig_ESO_pzmap_combined.png')

    # ===========================================================
    # 打印数值汇总
    # ===========================================================
    print('\n' + '=' * 70)
    print(' 极点数值汇总')
    print('=' * 70)

    print('\n--- 观测器误差动态极点 ---')
    for order in [3, 4]:
        for omega_ob in omega_ob_values:
            A = observer_A_matrix(order, omega_ob)
            eigvals = np.linalg.eigvals(A)
            poles_str = ', '.join([f'{p.real:.2f}+j{p.imag:.2f}' if abs(p.imag) > 0.01
                                   else f'{p.real:.2f}' for p in eigvals])
            print(f'  {order}阶 ESO, omega_ob={omega_ob:4d} rad/s: poles = [{poles_str}]')

    print('\n--- 速度环闭环极点 ---')
    for zeta, clbw, label, *_ in test_cases:
        poles, omega_cl, Kp_sp, Ki_sp = speed_loop_closed_loop_poles(zeta, clbw)
        poles_str = ', '.join([f'{p.real:.1f}+j{p.imag:.1f}' if abs(p.imag) > 0.01
                               else f'{p.real:.1f}' for p in poles])
        print(f'  {label}: poles = [{poles_str}]')
        print(f'    Current loop BW: {omega_cl:.0f} rad/s, '
              f'Speed Kp={Kp_sp:.4f}, Ki={Ki_sp:.4f}')

    plt.show()
    plt.close('all')
    print('\nDone!')
