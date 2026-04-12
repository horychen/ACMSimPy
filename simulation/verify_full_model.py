# -*- coding: utf-8 -*-
"""
验证2：构建完整的连续时间闭环模型（含 ESO 速度反馈 + 扰动前馈）
=================================================================
当前解析模型的错误假设：
  - 速度反馈 = 真实速度
  - ESO 只影响前馈通道

真实系统：
  - 速度反馈 = ESO 估计速度 (xS[1])
  - ESO 同时影响速度反馈和扰动前馈

本脚本验证：凹坑是否在正确的连续时间模型中就存在？
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import control
import copy, time

from tuner import (
    get_coeffs_dc_motor_current_regulator,
    get_coeffs_dc_motor_SPEED_regulator,
)
from sim_bode_sweep_v2 import MOTOR_PARAMS, run_single_frequency

d = copy.deepcopy(MOTOR_PARAMS)
R    = d['init_R']
L    = d['init_Lq']
J_s  = d['init_Js']
n_pp = d['init_npp']
KE   = d['init_KE']
KA   = KE

zeta = 15
CLBW_Hz = 1000

currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R, L, CLBW_Hz)
alpha_c = currentKp / L  # 电流环带宽 [rad/s]
speedKp, speedKi = get_coeffs_dc_motor_SPEED_regulator(J_s, n_pp, KA, zeta, alpha_c)

s = control.tf('s')
KT = 1.5 * n_pp * KA
rpm_per_radps = 30.0 / np.pi

# ======================================================================
# 基础传递函数
# ======================================================================
G_CL = alpha_c / (s + alpha_c)               # 电流闭环
G_plant_elec = n_pp / (J_s * s)               # iq → elec.rad/s
G_plant_rpm = G_plant_elec * rpm_per_radps / n_pp  # iq → RPM (mech)
G_PI = speedKp + speedKp * speedKi / s        # Speed PI (RPM误差 → iq*)

# 不含 ESO 的开环和敏感度
L_open = G_plant_rpm * G_CL * KT * G_PI       # 开环增益: RPM误差 → RPM
S_noESO = 1 / (1 + L_open)                    # 灵敏度函数

# TLoad → RPM (无 ESO)
G_TL_to_RPM = rpm_per_radps / (J_s / n_pp) / n_pp / s  # -1/(Js) * rpm_conv
G_dist_noESO = G_TL_to_RPM * S_noESO


def build_full_system_with_ESO(omega_ob):
    """
    构建完整的连续时间闭环系统：
      Plant + 4th-order ESO (提供速度反馈 + 扰动前馈) + Speed PI + Current Loop

    返回: 从 TLoad 到 speed(RPM) 的闭环传递函数
    """
    # ------------------------------------------------------------------
    # ESO 状态空间模型
    # ------------------------------------------------------------------
    # 增广系统状态: [theta_elec, omega_elec, T_d, T_d_dot]
    # 测量: y = theta_elec
    # 
    # 增广系统动力学:
    #   d(theta)   /dt = omega
    #   d(omega)   /dt = (npp/Js) * (Tem + T_d)
    #   d(T_d)     /dt = T_d_dot
    #   d(T_d_dot) /dt = 0  (assumed)
    #
    # ESO 观测器增益 (from code):
    #   ell1 = 4*wob
    #   ell2 = 6*wob^2
    #   ell3 = 4*wob^3 * Js/npp
    #   ell4 =   wob^4 * Js/npp
    #
    # 观测器:
    #   d(x0_hat)/dt = ell1*(theta - x0_hat) + x1_hat
    #   d(x1_hat)/dt = ell2*(theta - x0_hat) + (npp/Js)*(Tem + x2_hat)
    #   d(x2_hat)/dt = ell3*(theta - x0_hat) + x3_hat
    #   d(x3_hat)/dt = ell4*(theta - x0_hat)
    #
    # 误差动力学 (e = x - x_hat):
    #   de0/dt = -ell1*e0 + e1
    #   de1/dt = -ell2*e0 + (npp/Js)*e2
    #   de2/dt = -ell3*e0 + e3
    #   de3/dt = -ell4*e0 + T_d_ddot  (= 0 for constant T_d_dot)
    #
    # 特征方程: (s + wob)^4

    w = omega_ob
    ell1 = 4 * w
    ell2 = 6 * w**2
    ell3 = 4 * w**3 * J_s / n_pp
    ell4 =     w**4 * J_s / n_pp

    # ------------------------------------------------------------------
    # 方法: 直接推导 ESO 估计的传递函数
    # ------------------------------------------------------------------
    # 
    # ESO 的输入: theta (测量), Tem (已知)
    # ESO 的输出: x1_hat (速度估计), x2_hat (扰动估计)
    #
    # 关键关系:
    #   x1_hat = omega - e1
    #   x2_hat = T_d - e2
    #
    # 误差 e = [e0, e1, e2, e3] 的传递函数矩阵:
    #   从 T_d 到 e_i 可以通过解误差动力学得到
    #
    # 误差系统是自治的（如果 T_d 是其中一个状态），但作为外部输入处理时：
    #   T_d 通过增广系统间接影响 theta（测量），
    #    而测量误差驱动 ESO 校正

    # ------------------------------------------------------------------
    # 用 transfer function 直接计算
    # ------------------------------------------------------------------
    # 
    # 1. ESO 残差: H_T(s) = e2/T_d = s^4 / (s+w)^4
    # 2. ESO 速度误差: e1/T_d = ?
    #
    # 从误差矩阵 A_err = [[-ell1, 1, 0, 0],
    #                      [-ell2, 0, npp/Js, 0],
    #                      [-ell3, 0, 0, 1],
    #                      [-ell4, 0, 0, 0]]
    # 
    # 输入矩阵 B_err = ... (T_d 作为输入进入 e2 的导数方程)
    # 
    # 实际上，更精确的做法是用状态空间在 python-control 中计算

    # -- 构建误差动力学的状态空间 --
    # 
    # 真实系统 (只看机械部分):
    #   d(theta)/dt = omega_elec
    #   d(omega_elec)/dt = (npp/Js) * (Tem + T_d)
    #
    # ESO:
    #   d(x0)/dt = ell1*(theta - x0) + x1
    #   d(x1)/dt = ell2*(theta - x0) + (npp/Js)*(Tem + x2)
    #   d(x2)/dt = ell3*(theta - x0) + x3
    #   d(x3)/dt = ell4*(theta - x0)
    #
    # 选择总状态: [theta, omega, x0, x1, x2, x3]
    # 输入: [Tem, T_d]
    # 输出: [omega_mech_rpm, x1 (speed_estimate), x2 (disturbance_estimate)]

    kk = n_pp / J_s  # 转换系数

    A = np.array([
        # theta  omega    x0       x1      x2     x3
        [  0,     1,       0,       0,      0,      0],   # d(theta)/dt = omega
        [  0,     0,       0,       0,      0,      0],   # d(omega)/dt = kk*(Tem + T_d), 通过 B 输入
        [  ell1,  0,      -ell1,    1,      0,      0],   # d(x0)/dt = ell1*(theta-x0) + x1
        [  ell2,  0,      -ell2,    0,      kk,     0],   # d(x1)/dt = ell2*(theta-x0) + kk*(Tem + x2)
        [  ell3,  0,      -ell3,    0,      0,      1],   # d(x2)/dt = ell3*(theta-x0) + x3
        [  ell4,  0,      -ell4,    0,      0,      0],   # d(x3)/dt = ell4*(theta-x0)
    ], dtype=float)

    # B 矩阵: 输入 = [Tem, T_d]
    B = np.array([
        [0,    0  ],   # theta
        [kk,   kk ],   # omega: d(omega)/dt = kk*(Tem + T_d)
        [0,    0  ],   # x0
        [kk,   0  ],   # x1: d(x1)/dt 里的 kk*Tem (Tem 是已知输入)
        [0,    0  ],   # x2
        [0,    0  ],   # x3
    ], dtype=float)

    # C 矩阵: 输出 = [omega_mech_rpm, speed_estimate_elec, disturbance_estimate]
    C = np.array([
        # theta  omega  x0  x1  x2  x3
        [0, rpm_per_radps/n_pp, 0, 0, 0, 0],   # 实际 RPM = omega_elec * rpm/(npp)
        [0, 0, 0, 1, 0, 0],                      # ESO 速度估计 (elec.rad/s)
        [0, 0, 0, 0, 1, 0],                      # ESO 扰动估计
    ], dtype=float)

    D = np.zeros((3, 2))

    plant_eso = control.ss(A, B, C, D)

    # 从这个状态空间提取传递函数矩阵
    # G[i,j] = C[i] * (sI-A)^-1 * B[:,j]
    # (We compute frequency response directly via matrix inversion below)

    # 输出0 = omega_rpm, 输出1 = speed_est_elec, 输出2 = Td_est
    # 输入0 = Tem, 输入1 = T_d

    # T_d → omega_rpm (开环, 无PI)
    # T_d → speed_est_elec (ESO速度估计)
    # T_d → Td_est (ESO扰动估计)

    # 闭环构建:
    # 控制律: iq* = PI(s) * (omega_ref_elec - speed_est_elec) - Td_est / KT
    # Tem = KT * G_CL(s) * iq*
    # omega_ref_elec = omega_ref_rpm * 2*pi/60 * npp

    # ========== 用频率点逐一计算闭环 TF ==========
    # 因为 python-control 对复杂互联的传递函数处理有时不够稳健，
    # 我们直接在每个频率点上数值计算闭环传递函数

    freqs_dense = np.logspace(np.log10(1), np.log10(300), 500)
    mag_correct = np.zeros(len(freqs_dense))
    mag_simple = np.zeros(len(freqs_dense))
    mag_noESO = np.zeros(len(freqs_dense))

    for i, f in enumerate(freqs_dense):
        jw = 1j * 2 * np.pi * f

        # 基础传递函数在 jw 处的值
        G_CL_val = alpha_c / (jw + alpha_c)
        G_PI_val = speedKp + speedKp * speedKi / jw  # RPM误差 → iq*

        # 从状态空间计算传递函数矩阵的频率响应
        # plant_eso = ss(A, B, C, D)
        # H(jw) = C * (jwI - A)^(-1) * B + D
        H = C @ np.linalg.solve(jw * np.eye(6) - A, B) + D

        # H[0,0] = Tem → omega_rpm
        # H[0,1] = T_d → omega_rpm (开环植物响应)
        # H[1,0] = Tem → speed_est_elec
        # H[1,1] = T_d → speed_est_elec
        # H[2,0] = Tem → Td_est
        # H[2,1] = T_d → Td_est

        # --- 正确的闭环 (ESO提供速度反馈+扰动前馈) ---
        # iq* = G_PI * (omega_ref_elec - speed_est_elec) - Td_est / KT
        # Tem = KT * G_CL * iq*
        #
        # 对于扰动通道: omega_ref = 0
        # 令 Td 为输入，解 omega_rpm 为输出
        #
        # omega_rpm = H[0,0]*Tem + H[0,1]*Td
        # speed_est = H[1,0]*Tem + H[1,1]*Td
        # Td_est    = H[2,0]*Tem + H[2,1]*Td
        # iq*  = G_PI * (0 - speed_est * rpm_conv) - Td_est/KT
        #       其中 speed_est 的单位是 elec.rad/s
        #       PI 的输入是 RPM，所以要转换: speed_est_rpm = speed_est * rpm_per_radps / n_pp
        # Tem = KT * G_CL * iq*

        rpm_conv = rpm_per_radps / n_pp  # elec.rad/s → RPM

        # 把方程化成 Tem = f(Td):
        # speed_est_rpm = (H[1,0]*Tem + H[1,1]*Td) * rpm_conv
        # Td_est = H[2,0]*Tem + H[2,1]*Td
        # iq* = -G_PI * (H[1,0]*Tem + H[1,1]*Td)*rpm_conv - (H[2,0]*Tem + H[2,1]*Td)/KT
        # Tem = KT * G_CL * iq*
        #     = KT * G_CL * [-G_PI*rpm_conv*(H[1,0]*Tem + H[1,1]*Td) - (H[2,0]*Tem + H[2,1]*Td)/KT]
        #     = -KT*G_CL*G_PI*rpm_conv*H[1,0]*Tem - KT*G_CL*G_PI*rpm_conv*H[1,1]*Td
        #       - G_CL*H[2,0]*Tem - G_CL*H[2,1]*Td
        # Tem * (1 + KT*G_CL*G_PI*rpm_conv*H[1,0] + G_CL*H[2,0])
        #     = -Td * (KT*G_CL*G_PI*rpm_conv*H[1,1] + G_CL*H[2,1])

        denom = 1 + KT * G_CL_val * G_PI_val * rpm_conv * H[1,0] + G_CL_val * H[2,0]
        numer = -(KT * G_CL_val * G_PI_val * rpm_conv * H[1,1] + G_CL_val * H[2,1])
        Tem_over_Td = numer / denom

        # omega_rpm = H[0,0]*Tem + H[0,1]*Td
        # omega_rpm/Td = H[0,0]*(Tem/Td) + H[0,1]
        G_dist_correct = H[0,0] * Tem_over_Td + H[0,1]
        mag_correct[i] = 20 * np.log10(max(abs(G_dist_correct), 1e-30))

        # --- 简化模型 (当前解析: ESO只影响前馈, 速度反馈=真实速度) ---
        # 这就是 get_analytical_bode() 里的模型
        # G_dist_simple = G_TL_to_RPM / (1 + L_open) * s^4/(s+w)^4
        P_val = rpm_per_radps / (J_s / n_pp) / n_pp / jw  # TL → RPM
        L_val = P_val * KT * G_CL_val * G_PI_val  # 开环增益
        eso_residual = jw**4 / (jw + w)**4
        G_dist_simple_val = P_val / (1 + L_val) * eso_residual
        mag_simple[i] = 20 * np.log10(max(abs(G_dist_simple_val), 1e-30))

        # --- 无 ESO ---
        G_dist_noESO_val = P_val / (1 + L_val)
        mag_noESO[i] = 20 * np.log10(max(abs(G_dist_noESO_val), 1e-30))

    return freqs_dense, mag_correct, mag_simple, mag_noESO


# ======================================================================
# 多个 omega_ob 并行计算
# ======================================================================
omega_obs = [100, 200, 400]
colors_eso = ['#e74c3c', '#3498db', '#2ecc71']

plt.style.use('bmh')
mpl.rc('font', family='Times New Roman', size=10)

fig, axes = plt.subplots(3, 1, figsize=(12, 14), sharex=True)

for omega_ob, color in zip(omega_obs, colors_eso):
    freqs_dense, mag_correct, mag_simple, mag_noESO = build_full_system_with_ESO(omega_ob)
    bw_hz = omega_ob / (2 * np.pi)

    # -- 上图: 三条解析曲线对比 --
    ax = axes[0]
    ax.plot(freqs_dense, mag_correct, '-', color=color, lw=1.5,
            label=f'Full Model (ESO fb+ff) $\\omega_{{ob}}$={omega_ob}')
    ax.plot(freqs_dense, mag_simple, '--', color=color, lw=1, alpha=0.5,
            label=f'_Simplified (ESO ff only) $\\omega_{{ob}}$={omega_ob}')
    ax.axvline(bw_hz, color=color, ls=':', lw=0.8, alpha=0.5)

ax = axes[0]
ax.plot(freqs_dense, mag_noESO, 'k--', lw=1, alpha=0.4, label='No ESO')
ax.set_ylabel('Magnitude [dB]')
ax.set_title('Continuous-Time Analytical: Full Model vs Simplified Model')
ax.legend(fontsize=7, ncol=2, loc='lower right')
ax.grid(True, alpha=0.3)

# -- 中图: Full - Simplified (体现 ESO 速度反馈的额外效应) --
ax = axes[1]
for omega_ob, color in zip(omega_obs, colors_eso):
    freqs_dense, mag_correct, mag_simple, mag_noESO = build_full_system_with_ESO(omega_ob)
    bw_hz = omega_ob / (2 * np.pi)
    ax.plot(freqs_dense, mag_correct - mag_simple, '-', color=color, lw=1.5,
            label=f'$\\omega_{{ob}}$={omega_ob} ({bw_hz:.1f} Hz)')
    ax.axvline(bw_hz, color=color, ls=':', lw=0.8, alpha=0.5)
ax.axhline(0, color='k', ls='-', lw=0.5)
ax.set_ylabel('Full − Simplified [dB]')
ax.set_title('Impact of ESO Speed Feedback (difference from simplified model)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# -- 下图: 与数值仿真对比 --
ax = axes[2]
for omega_ob, color in zip(omega_obs, colors_eso):
    freqs_dense, mag_correct, mag_simple, mag_noESO = build_full_system_with_ESO(omega_ob)
    bw_hz = omega_ob / (2 * np.pi)

    # 正确的连续时间模型
    ax.plot(freqs_dense, mag_correct, '-', color=color, lw=1.5, alpha=0.7,
            label=f'Analytical (Full) $\\omega_{{ob}}$={omega_ob}')
    ax.axvline(bw_hz, color=color, ls=':', lw=0.8, alpha=0.5)

# 数值仿真
freqs_sim = np.linspace(5, 100, 50)
for omega_ob, color in zip(omega_obs, colors_eso):
    print(f'\nSim sweep: omega_ob={omega_ob} ...')
    sim_mag = []
    t0 = time.time()
    for f in freqs_sim:
        m, _ = run_single_frequency(d, f, zeta, CLBW_Hz, True, omega_ob,
                                     mode='disturbance', rpm_0=500.0)
        sim_mag.append(m)
    sim_mag = np.array(sim_mag)
    print(f'  Done in {time.time()-t0:.1f}s')
    ax.plot(freqs_sim, sim_mag, 'o', color=color, ms=3, alpha=0.8,
            label=f'Sim $\\omega_{{ob}}$={omega_ob}')

ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('Magnitude [dB]')
ax.set_title('Full Continuous-Time Model vs Numerical Simulation')
ax.legend(fontsize=7, ncol=2, loc='lower right')
ax.grid(True, alpha=0.3)

fig.suptitle(
    'ESO Disturbance Rejection: The Missing Speed-Feedback Effect\n'
    f'(ζ={zeta}, CLBW={CLBW_Hz}Hz, KT={KT:.4f}, Js={J_s:.2e})',
    fontsize=12, fontweight='bold'
)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('fig_full_vs_simplified_model.png', dpi=150, bbox_inches='tight')
print('\nSaved: fig_full_vs_simplified_model.png')
