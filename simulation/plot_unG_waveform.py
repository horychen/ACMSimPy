"""
SVPWM u_nG (中性点对地电压/共模电压) 波形可视化
===============================================
本脚本基于 tutorials_ep3_svpwm.py 的仿真框架，
修改 watch_data 以记录并绘制 u_nG 波形。

u_nG = (u_AG + u_BG + u_CG) / 3
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for saving figures
import matplotlib.pyplot as plt

# 配置中文字体
plt.rcParams['font.sans-serif'] = ['Heiti TC', 'STHeiti', 'Hiragino Sans', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# Import the simulation module
sys.path.insert(0, os.path.dirname(__file__))

# We need to reimplement a simplified version that captures u_nG
# Let's import the building blocks from tutorials_ep3_svpwm
from tutorials_ep3_svpwm import (
    The_Motor_Controller,
    The_AC_Machine,
    The_PI_Regulator,
    SVgen_Object,
    SVGEN_DQ,
    gate_signal_generator,
    RK4_MACHINE,
    DSP,
)
from numba import njit
import numba

@njit(nogil=True)
def simulate_with_unG(
    t0, TIME,
    ACM, CTRL,
    reg_id, reg_iq, reg_speed,
):
    """Modified simulation loop that records u_nG."""
    MACHINE_TS = CTRL.CL_TS / ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    controller_down_sampling_ceiling = int(CTRL.CL_TS / MACHINE_TS)

    CPU_TICK_PER_SAMPLING_PERIOD = ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    DEAD_TIME_AS_COUNT = int(200 * 0.5e-4 * CPU_TICK_PER_SAMPLING_PERIOD)
    Vdc = 150
    one_over_Vdc = 1 / Vdc
    svgen1 = SVgen_Object(CPU_TICK_PER_SAMPLING_PERIOD)

    machine_times = np.arange(t0, t0 + TIME, MACHINE_TS)
    # Channels: 0=u_AG, 1=u_BG, 2=u_CG, 3=u_nG, 4=u_An, 5=u_Bn, 6=u_Cn,
    #           7=Ta, 8=Tb, 9=Tc, 10=S1, 11=S2, 12=S3, 13=carrier
    #           14=line_AB, 15=line_BC, 16=line_AC
    #           17=ACM.iAlfa, 18=ACM.iBeta
    watch_data = np.zeros((20, len(machine_times)))

    jj = controller_down_sampling_ceiling
    watch_index = 0
    for ii in range(len(machine_times)):
        t = machine_times[ii]

        RK4_MACHINE(t, ACM, hs=MACHINE_TS)

        ACM.theta_d_mech = ACM.x[0]
        ACM.omega_r_mech = ACM.x[1]
        ACM.KA           = ACM.x[2]
        ACM.iD           = ACM.x[3]
        ACM.iQ           = ACM.x[4]
        ACM.theta_d      = ACM.theta_d_mech * ACM.npp
        ACM.omega_r_elec = ACM.omega_r_mech * ACM.npp
        ACM.omega_syn    = ACM.omega_r_elec + ACM.omega_slip

        ACM.cosT = np.cos(ACM.theta_d)
        ACM.sinT = np.sin(ACM.theta_d)
        ACM.iAlfa = ACM.iD * ACM.cosT + ACM.iQ * -ACM.sinT
        ACM.iBeta = ACM.iD * ACM.sinT + ACM.iQ *  ACM.cosT

        jj += 1
        if jj >= controller_down_sampling_ceiling:
            jj = 0

            if CTRL.bool_overwrite_speed_commands == False:
                if t < 0.5:
                    CTRL.cmd_rpm = 200
                pass

            DSP(ACM=ACM, CTRL=CTRL, reg_speed=reg_speed, reg_id=reg_id, reg_iq=reg_iq)

            svgen1.Ualfa = CTRL.cmd_uab[0]
            svgen1.Ubeta = CTRL.cmd_uab[1]
            SVGEN_DQ(svgen1, one_over_Vdc)
            # 仿真中反转回来
            svgen1.Ta = 1 - svgen1.Ta
            svgen1.Tb = 1 - svgen1.Tb
            svgen1.Tc = 1 - svgen1.Tc
            svgen1.EPwm1Regs_CMPA_bit_CMPA = int(svgen1.Ta * CPU_TICK_PER_SAMPLING_PERIOD * 0.5)
            svgen1.EPwm2Regs_CMPA_bit_CMPA = int(svgen1.Tb * CPU_TICK_PER_SAMPLING_PERIOD * 0.5)
            svgen1.EPwm3Regs_CMPA_bit_CMPA = int(svgen1.Tc * CPU_TICK_PER_SAMPLING_PERIOD * 0.5)
            svgen1.bool_interupt_event = True

        if CPU_TICK_PER_SAMPLING_PERIOD >= 20:
            ACM.ia = ACM.iAlfa
            ACM.ib = ACM.iAlfa * -0.5 + ACM.iBeta * 0.8660254
            ACM.ic = ACM.iAlfa * -0.5 + ACM.iBeta * -0.8660254

            gate_signal_generator(ii, svgen1,
                                  CPU_TICK_PER_SAMPLING_PERIOD=CPU_TICK_PER_SAMPLING_PERIOD,
                                  DEAD_TIME_AS_COUNT=DEAD_TIME_AS_COUNT)

            if svgen1.S1 == True:
                svgen1.voltage_potential_at_terminal[0] = Vdc
            elif svgen1.S4 == True:
                svgen1.voltage_potential_at_terminal[0] = 0
            else:
                svgen1.voltage_potential_at_terminal[0] = Vdc if ACM.ia < 0 else 0

            if svgen1.S2 == True:
                svgen1.voltage_potential_at_terminal[1] = Vdc
            elif svgen1.S5 == True:
                svgen1.voltage_potential_at_terminal[1] = 0
            else:
                svgen1.voltage_potential_at_terminal[1] = Vdc if ACM.ib < 0 else 0

            if svgen1.S3 == True:
                svgen1.voltage_potential_at_terminal[2] = Vdc
            elif svgen1.S6 == True:
                svgen1.voltage_potential_at_terminal[2] = 0
            else:
                svgen1.voltage_potential_at_terminal[2] = Vdc if ACM.ic < 0 else 0

            svgen1.line_to_line_voltage_AC = svgen1.voltage_potential_at_terminal[0] - svgen1.voltage_potential_at_terminal[2]
            svgen1.line_to_line_voltage_BC = svgen1.voltage_potential_at_terminal[1] - svgen1.voltage_potential_at_terminal[2]
            svgen1.line_to_line_voltage_AB = svgen1.voltage_potential_at_terminal[0] - svgen1.voltage_potential_at_terminal[1]

            ACM.uab[0] = svgen1.line_to_line_voltage_AC * 0.6666667 - (svgen1.line_to_line_voltage_BC + 0) * 0.3333333
            ACM.uab[1] = 0.577350269 * (svgen1.line_to_line_voltage_BC - 0)
        else:
            ACM.uab[0] = CTRL.cmd_uab[0]
            ACM.uab[1] = CTRL.cmd_uab[1]

        ACM.udq[0] = ACM.uab[0] *  ACM.cosT + ACM.uab[1] * ACM.sinT
        ACM.udq[1] = ACM.uab[0] * -ACM.sinT + ACM.uab[1] * ACM.cosT

        # Record data
        u_AG = svgen1.voltage_potential_at_terminal[0]
        u_BG = svgen1.voltage_potential_at_terminal[1]
        u_CG = svgen1.voltage_potential_at_terminal[2]
        u_nG = (u_AG + u_BG + u_CG) / 3.0

        watch_data[0][watch_index] = u_AG
        watch_data[1][watch_index] = u_BG
        watch_data[2][watch_index] = u_CG
        watch_data[3][watch_index] = u_nG                # ★ 核心: 共模电压
        watch_data[4][watch_index] = u_AG - u_nG          # u_An (A相电压)
        watch_data[5][watch_index] = u_BG - u_nG          # u_Bn
        watch_data[6][watch_index] = u_CG - u_nG          # u_Cn
        watch_data[7][watch_index] = svgen1.Ta
        watch_data[8][watch_index] = svgen1.Tb
        watch_data[9][watch_index] = svgen1.Tc
        watch_data[10][watch_index] = svgen1.S1
        watch_data[11][watch_index] = svgen1.S2
        watch_data[12][watch_index] = svgen1.S3
        watch_data[13][watch_index] = svgen1.carrier_counter
        watch_data[14][watch_index] = svgen1.line_to_line_voltage_AB
        watch_data[15][watch_index] = svgen1.line_to_line_voltage_BC
        watch_data[16][watch_index] = svgen1.line_to_line_voltage_AC
        watch_data[17][watch_index] = ACM.iAlfa
        watch_data[18][watch_index] = ACM.iBeta
        watch_data[19][watch_index] = ACM.omega_r_mech / (2*np.pi) * 60  # rpm
        watch_index += 1

    return machine_times, watch_data


def main():
    print("Initializing simulation...")

    CL_TS = 1e-4
    TIME_SLICE = 0.1  # 100ms — enough to see several electrical cycles

    CTRL = The_Motor_Controller(CL_TS, 5*CL_TS,
        init_npp=21,
        init_IN=72/1.414,
        init_R=0.1222,
        init_Ld=0.000502,
        init_Lq=0.000571,
        init_KE=0.188492,
        init_Rreq=0.0,
        init_Js=0.203)

    ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=100)

    reg_id    = The_PI_Regulator(0.737168, 0.737168*214.011*CTRL.CL_TS, 150/1.732)
    reg_iq    = The_PI_Regulator(0.737168, 0.737168*214.011*CTRL.CL_TS, 150/1.732)
    reg_speed = The_PI_Regulator(0.323363, 0.323363*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)

    print("Running simulation (with JIT compilation on first run)...")
    machine_times, watch_data = simulate_with_unG(
        t0=0, TIME=TIME_SLICE,
        ACM=ACM, CTRL=CTRL,
        reg_id=reg_id, reg_iq=reg_iq, reg_speed=reg_speed)

    print(f"Simulation complete. {len(machine_times)} data points.")

    # Extract data
    u_AG = watch_data[0]
    u_BG = watch_data[1]
    u_CG = watch_data[2]
    u_nG = watch_data[3]
    u_An = watch_data[4]
    u_Bn = watch_data[5]
    u_Cn = watch_data[6]
    Ta   = watch_data[7]
    Tb   = watch_data[8]
    Tc   = watch_data[9]
    line_AB = watch_data[14]
    line_BC = watch_data[15]
    rpm  = watch_data[19]

    t_ms = machine_times * 1000  # convert to ms

    # --- Plot ---
    fig, axes = plt.subplots(5, 1, figsize=(16, 18), sharex=True)
    fig.suptitle('SVPWM 电压波形分析 (Vdc = 150V)', fontsize=16, fontweight='bold')

    # 1. Terminal voltages (端电势)
    ax = axes[0]
    ax.plot(t_ms, u_AG, label='$u_{AG}$', alpha=0.7, linewidth=0.3, color='#e74c3c')
    ax.plot(t_ms, u_BG, label='$u_{BG}$', alpha=0.7, linewidth=0.3, color='#2ecc71')
    ax.plot(t_ms, u_CG, label='$u_{CG}$', alpha=0.7, linewidth=0.3, color='#3498db')
    ax.set_ylabel('端电势 [V]')
    ax.set_title('① 端电势 $u_{xG}$ (x=A,B,C)')
    ax.legend(loc='upper right')
    ax.set_ylim([-20, 170])
    ax.grid(True, alpha=0.3)

    # 2. u_nG (共模电压) — THE STAR
    ax = axes[1]
    ax.plot(t_ms, u_nG, color='#e67e22', linewidth=0.5, label='$u_{nG}$')
    ax.axhline(y=75, color='gray', linestyle='--', alpha=0.5, label='$V_{dc}/2$')
    ax.set_ylabel('$u_{nG}$ [V]')
    ax.set_title('② 共模电压 $u_{nG} = (u_{AG}+u_{BG}+u_{CG})/3$  ← 马鞍形波形')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # 3. Phase voltages (相电压 = 端电势 - 共模)
    ax = axes[2]
    ax.plot(t_ms, u_An, label='$u_{An}=u_{AG}-u_{nG}$', alpha=0.7, linewidth=0.3, color='#e74c3c')
    ax.plot(t_ms, u_Bn, label='$u_{Bn}=u_{BG}-u_{nG}$', alpha=0.7, linewidth=0.3, color='#2ecc71')
    ax.plot(t_ms, u_Cn, label='$u_{Cn}=u_{CG}-u_{nG}$', alpha=0.7, linewidth=0.3, color='#3498db')
    ax.set_ylabel('相电压 [V]')
    ax.set_title('③ 相电压 $u_{xn}$ (绕组上的实际电压)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # 4. Line-to-line voltages (线电压)
    ax = axes[3]
    ax.plot(t_ms, line_AB, label='$u_{AB}$', alpha=0.7, linewidth=0.3, color='#9b59b6')
    ax.plot(t_ms, line_BC, label='$u_{BC}$', alpha=0.7, linewidth=0.3, color='#1abc9c')
    ax.set_ylabel('线电压 [V]')
    ax.set_title('④ 线电压 (共模已被消除)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # 5. Duty cycles
    ax = axes[4]
    ax.plot(t_ms, Ta, label='$T_a$ (A相占空比)', alpha=0.8, linewidth=0.5, color='#e74c3c')
    ax.plot(t_ms, Tb, label='$T_b$ (B相占空比)', alpha=0.8, linewidth=0.5, color='#2ecc71')
    ax.plot(t_ms, Tc, label='$T_c$ (C相占空比)', alpha=0.8, linewidth=0.5, color='#3498db')
    ax.set_ylabel('占空比')
    ax.set_xlabel('时间 [ms]')
    ax.set_title('⑤ SVPWM 占空比 Ta, Tb, Tc')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([-0.05, 1.05])

    plt.tight_layout()
    output_path = os.path.join(os.path.dirname(__file__), 'svpwm_unG_waveforms.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to: {output_path}")

    # --- Zoomed plot for u_nG detail ---
    fig2, axes2 = plt.subplots(3, 1, figsize=(16, 10), sharex=True)
    fig2.suptitle('SVPWM u_nG 局部放大 (观察马鞍形特征)', fontsize=14, fontweight='bold')

    # Zoom to ~2 electrical cycles (pick a window where motor is running)
    # Electrical frequency at 200rpm with 21 pole pairs:
    # f_e = 200/60 * 21 = 70 Hz, T_e = 14.3 ms
    # Show ~ 30ms
    t_start_ms = 50
    t_end_ms = 90
    mask = (t_ms >= t_start_ms) & (t_ms <= t_end_ms)

    ax = axes2[0]
    ax.plot(t_ms[mask], u_AG[mask], label='$u_{AG}$', alpha=0.7, linewidth=0.4, color='#e74c3c')
    ax.plot(t_ms[mask], u_BG[mask], label='$u_{BG}$', alpha=0.7, linewidth=0.4, color='#2ecc71')
    ax.plot(t_ms[mask], u_CG[mask], label='$u_{CG}$', alpha=0.7, linewidth=0.4, color='#3498db')
    ax.set_ylabel('端电势 [V]')
    ax.set_title('端电势 (放大)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    ax = axes2[1]
    ax.plot(t_ms[mask], u_nG[mask], color='#e67e22', linewidth=0.6, label='$u_{nG}$')
    ax.axhline(y=75, color='gray', linestyle='--', alpha=0.5, label='$V_{dc}/2 = 75V$')
    ax.axhline(y=150/6, color='lightblue', linestyle=':', alpha=0.5, label='$V_{dc}/6 = 25V$')
    ax.axhline(y=150*5/6, color='lightcoral', linestyle=':', alpha=0.5, label='$5V_{dc}/6 = 125V$')
    ax.set_ylabel('$u_{nG}$ [V]')
    ax.set_title('共模电压 u_nG (放大) — 观察马鞍形包络')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    ax = axes2[2]
    ax.plot(t_ms[mask], u_An[mask], label='$u_{An}$', alpha=0.7, linewidth=0.4, color='#e74c3c')
    ax.plot(t_ms[mask], u_Bn[mask], label='$u_{Bn}$', alpha=0.7, linewidth=0.4, color='#2ecc71')
    ax.plot(t_ms[mask], u_Cn[mask], label='$u_{Cn}$', alpha=0.7, linewidth=0.4, color='#3498db')
    ax.set_ylabel('相电压 [V]')
    ax.set_xlabel('时间 [ms]')
    ax.set_title('相电压 (放大) — 减去 u_nG 后的准正弦波')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_path2 = os.path.join(os.path.dirname(__file__), 'svpwm_unG_zoomed.png')
    plt.savefig(output_path2, dpi=150, bbox_inches='tight')
    print(f"Zoomed figure saved to: {output_path2}")

    # --- Exactly 2 electrical cycles ---
    # f_e = 200/60 * 21 = 70 Hz, T_e = 1/70 ≈ 14.286 ms
    # 2 cycles ≈ 28.57 ms
    T_e_ms = 1000.0 / (200.0 / 60.0 * 21)  # ~14.286 ms
    t_start_2c = 60.0  # start at 60ms (motor already running steadily)
    t_end_2c = t_start_2c + 2 * T_e_ms

    mask2 = (t_ms >= t_start_2c) & (t_ms <= t_end_2c)

    fig3, axes3 = plt.subplots(5, 1, figsize=(18, 16), sharex=True,
                                gridspec_kw={'height_ratios': [1, 1, 1, 0.6, 0.6]})
    fig3.suptitle(f'SVPWM 精确 2 个电周期 ({t_start_2c:.1f}ms ~ {t_end_2c:.1f}ms, $T_e$={T_e_ms:.2f}ms)',
                  fontsize=14, fontweight='bold')

    # (a) Terminal voltages
    ax = axes3[0]
    ax.plot(t_ms[mask2], u_AG[mask2], label='$u_{AG}$', linewidth=0.5, color='#e74c3c')
    ax.plot(t_ms[mask2], u_BG[mask2], label='$u_{BG}$', linewidth=0.5, color='#2ecc71')
    ax.plot(t_ms[mask2], u_CG[mask2], label='$u_{CG}$', linewidth=0.5, color='#3498db')
    ax.set_ylabel('端电势 [V]')
    ax.set_title('(a) 端电势 $u_{xG}$：开关管导通时为 $V_{dc}$，关断时为 0')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim([-10, 160])
    ax.grid(True, alpha=0.3)

    # (b) u_nG - THE STAR
    ax = axes3[1]
    ax.plot(t_ms[mask2], u_nG[mask2], color='#e67e22', linewidth=0.6, label='$u_{nG}$')
    ax.axhline(y=75, color='gray', linestyle='--', alpha=0.5, linewidth=1, label='$V_{dc}/2 = 75V$')
    ax.axhline(y=150/6, color='#3498db', linestyle=':', alpha=0.6, linewidth=1, label='$V_{dc}/6 = 25V$')
    ax.axhline(y=150*5/6, color='#e74c3c', linestyle=':', alpha=0.6, linewidth=1, label='$5V_{dc}/6 = 125V$')
    ax.set_ylabel('$u_{nG}$ [V]')
    ax.set_title('(b) 共模电压 $u_{nG} = (u_{AG}+u_{BG}+u_{CG})/3$：马鞍形包络，含 3 次谐波')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim([-10, 160])
    ax.grid(True, alpha=0.3)

    # (c) Phase voltages
    ax = axes3[2]
    ax.plot(t_ms[mask2], u_An[mask2], label='$u_{An}=u_{AG}-u_{nG}$', linewidth=0.5, color='#e74c3c')
    ax.plot(t_ms[mask2], u_Bn[mask2], label='$u_{Bn}=u_{BG}-u_{nG}$', linewidth=0.5, color='#2ecc71')
    ax.plot(t_ms[mask2], u_Cn[mask2], label='$u_{Cn}=u_{CG}-u_{nG}$', linewidth=0.5, color='#3498db')
    ax.set_ylabel('相电压 [V]')
    ax.set_title('(c) 相电压：减去共模后，PWM 等效正弦波形')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    # (d) Duty cycles
    ax = axes3[3]
    ax.plot(t_ms[mask2], Ta[mask2], label='$T_a$', linewidth=0.8, color='#e74c3c')
    ax.plot(t_ms[mask2], Tb[mask2], label='$T_b$', linewidth=0.8, color='#2ecc71')
    ax.plot(t_ms[mask2], Tc[mask2], label='$T_c$', linewidth=0.8, color='#3498db')
    ax.set_ylabel('占空比')
    ax.set_title('(d) SVPWM 占空比：注意马鞍形调制波（三次谐波注入的体现）')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim([-0.05, 1.05])
    ax.grid(True, alpha=0.3)

    # (e) Line-to-line voltages
    ax = axes3[4]
    ax.plot(t_ms[mask2], line_AB[mask2], label='$u_{AB}$', linewidth=0.5, color='#9b59b6')
    ax.plot(t_ms[mask2], line_BC[mask2], label='$u_{BC}$', linewidth=0.5, color='#1abc9c')
    ax.set_ylabel('线电压 [V]')
    ax.set_xlabel('时间 [ms]')
    ax.set_title('(e) 线电压：共模 $u_{nG}$ 完全消除，仅含差模分量')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_path3 = os.path.join(os.path.dirname(__file__), 'svpwm_unG_2cycles.png')
    plt.savefig(output_path3, dpi=150, bbox_inches='tight')
    print(f"2-cycle figure saved to: {output_path3}")

    print("\nDone!")


if __name__ == '__main__':
    main()
