from tutorials_ep9_flux_estimator import *

# User input:
d = d_user_input_motor_dict = {
    # Timing
    'CL_TS': 1e-4,
    'TIME_SLICE': 0.2,
    'NUMBER_OF_SLICES': 2,
    'VL_EXE_PER_CL_EXE': 5,
    'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
    'CTRL.bool_apply_speed_closed_loop_control': True,
    'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
    'CTRL.bool_apply_sweeping_frequency_excitation': False,
    'CTRL.bool_overwrite_speed_commands': True,
    'CTRL.bool_zero_id_control': True,
    'FOC_delta': 10, # 25, # 6.5
    'FOC_desired_VLBW_HZ': 120, # 60
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

图 = 1 # 空载加速、加载、反转
# 小电感电机
d['init_npp'] = 22
d['init_IN'] = 1.3*6/1.414
d['init_R'] = 0.035
d['init_Ld'] = 1*0.036*1e-3
d['init_Lq'] = 1*0.036*1e-3
d['init_KE'] = 0.0125
d['init_Rreq'] = 0.0
d['init_Js'] = 0.44*1e-4
d['DC_BUS_VOLTAGE'] = 10
d['user_system_input_code'] = '''if ii < 1: CTRL.cmd_idq[0] = 0.0; CTRL.cmd_rpm = 150 \nelif ii <5: ACM.TLoad = 0.2 \nelif ii <100: CTRL.cmd_rpm = -150'''
def 图1画图代码():
    plt.style.use('bmh') # https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
    mpl.rc('font', family='Times New Roman', size=10.0)
    mpl.rc('legend', fontsize=10)
    mpl.rcParams['lines.linewidth'] = 0.75 # mpl.rc('lines', linewidth=4, linestyle='-.')
    mpl.rcParams['mathtext.fontset'] = 'stix'

    fig, axes = plt.subplots(nrows=7, ncols=1, dpi=150, facecolor='w', figsize=(8,12), sharex=True)

    ax = axes[0]
    ax.plot(global_machine_times, gdd['CTRL.cmd_rpm'], label=r'$\omega_r^*$')
    ax.plot(global_machine_times, gdd['CTRL.omega_r_mech'], label=r'$\omega_r$')
    ax.set_ylabel(r'Speed [r/min]', multialignment='center') #) #, fontdict=font)
    # ax.legend(loc=2, prop={'size': 6})
    ax.legend(loc=1, fontsize=6)

    ax = axes[1]
    ax.plot(global_machine_times, gdd['CTRL.cmd_idq[0]'], label=r'$i_d^*$')
    ax.plot(global_machine_times, gdd['CTRL.idq[0]'], label=r'$i_d$')
    # ax.plot(global_machine_times, gdd['ACM.iD'])
    ax.set_ylabel(r'$i_d$ [A]', multialignment='center') #, fontdict=font)

    ax = axes[2]
    ax.plot(global_machine_times, gdd['CTRL.cmd_idq[1]'], label=r'$i_q^*$')
    ax.plot(global_machine_times, gdd['CTRL.idq[1]'], label=r'$i_q$')
    ax.set_ylabel(r'$i_q$ [A]', multialignment='center') #, fontdict=font)

    ax = axes[3]
    ax.plot(global_machine_times, gdd['ACM.Tem'], label=r'ACM.$T_{\rm em}$')
    ax.plot(global_machine_times, gdd['CTRL.Tem'], label=r'CTRL.$T_{\rm em}$')
    ax.set_ylabel(r'$T_{\rm em}$ [Nm]', multialignment='center') #, fontdict=font)

    ax = axes[4]
    ax.plot(global_machine_times, (gdd['ACM.udq[0]']), label=r'$u_d$') # lpf1_inverter
    ax.plot(global_machine_times, gdd['CTRL.cmd_udq[0]'], label=r'$u_d^*$')
    ax.set_ylabel(r'$u_d$ [V]', multialignment='center') #, fontdict=font)

    ax = axes[5]
    ax.plot(global_machine_times, (gdd['ACM.udq[1]']), label=r'$u_q$') # lpf1_inverter
    ax.plot(global_machine_times, gdd['CTRL.cmd_udq[1]'], label=r'$u_q^*$')
    ax.set_ylabel(r'$u_q$ [V]', multialignment='center') #, fontdict=font)

    ax = axes[6]
    # ax.plot(global_machine_times, gdd['CTRL.cmd_uab[0]'], label=r'$u_\alpha$')
    # ax.plot(global_machine_times, gdd['CTRL.cmd_uab[1]'], label=r'$u_\beta$')
    # ax.set_ylabel(r'$u_{\alpha,\beta}$ [V]', multialignment='center') #, fontdict=font)
    ax.plot(global_machine_times, gdd['fe_htz.psi_2[0]'], label=r'$\psi_\alpha$')
    ax.plot(global_machine_times, gdd['fe_htz.psi_2[1]'], label=r'$\psi_\beta$')
    ax.plot(global_machine_times, gdd['fe_htz.u_offset[0]'], label=r'$u_{{\rm offset},\alpha}$')
    ax.plot(global_machine_times, gdd['fe_htz.u_offset[1]'], label=r'$u_{{\rm offset},\beta}$')
    # ax.set_ylabel(r'$\psi_2$ [Wb]', multialignment='center') #, fontdict=font)

    for ax in axes:
        ax.grid(True)
        ax.legend(loc=1)
        # for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        #     tick.label.set_font(font)
    axes[-1].set_xlabel('Time [s]') #, fontdict=font)
    return fig
sim1 = Simulation_Benchmark(d); gdd, global_machine_times = sim1.gdd, sim1.global_machine_times; fig = 图1画图代码(); # fig.savefig(f'SliceFSPM-fig-{图}.pdf', dpi=400, bbox_inches='tight', pad_inches=0)

plt.show()
quit()

图 = 2
d['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'] = 500
sim1 = Simulation_Benchmark(d); gdd, global_machine_times = sim1.gdd, sim1.global_machine_times; fig = 图1画图代码(); # fig.savefig(f'SliceFSPM-fig-{图}.pdf', dpi=400, bbox_inches='tight', pad_inches=0)
d['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'] = 1

plt.show()
quit()

图 = 3 # 空载加速、加载、反转
# 普通伺服电机
d['init_npp'] = 4
d['init_IN'] = 4
d['init_R'] = 1.1
d['init_Ld'] = 5e-3
d['init_Lq'] = 6e-3
d['init_KE'] = 0.1
d['init_Rreq'] = 0.0
d['init_Js'] = 0.008
d['DC_BUS_VOLTAGE'] = 80
d['user_system_input_code'] = '''if ii < 1: CTRL.cmd_idq[0] = 0.0; CTRL.cmd_rpm = 150 \nelif ii <5: ACM.TLoad = 1.27 \nelif ii <100: CTRL.cmd_rpm = -150'''
d['CL_SERIES_KP'] = None # activate auto-tuner
sim1 = Simulation_Benchmark(d); gdd, global_machine_times = sim1.gdd, sim1.global_machine_times; fig = 图1画图代码()

图 = 4
d['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'] = 1
sim1 = Simulation_Benchmark(d); gdd, global_machine_times = sim1.gdd, sim1.global_machine_times; fig = 图1画图代码(); # fig.savefig(f'SliceFSPM-fig-{图}.pdf', dpi=400, bbox_inches='tight', pad_inches=0)
d['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'] = 500

图 = 5 # 空载加速、加载、反转
# 大电感电机
d['init_npp'] = 24
d['init_IN'] = 4.93
d['init_R'] =  1.97
d['init_Ld'] = 0.1035
d['init_Lq'] = 0.1063
d['init_KE'] = 0.0745
d['init_Rreq'] = 0.0
d['init_Js'] = 1.5*0.051
d['DC_BUS_VOLTAGE'] = 800
d['user_system_input_code'] = '''if ii < 1: CTRL.cmd_idq[0] = 0.0; CTRL.cmd_rpm = 150 \nelif ii <5: ACM.TLoad = 20 \nelif ii <100: CTRL.cmd_rpm = -150'''
d['CL_SERIES_KP'] = None # activate auto-tuner
sim1 = Simulation_Benchmark(d); gdd, global_machine_times = sim1.gdd, sim1.global_machine_times; fig = 图1画图代码()

图 = 6
d['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'] = 1
sim1 = Simulation_Benchmark(d); gdd, global_machine_times = sim1.gdd, sim1.global_machine_times; fig = 图1画图代码(); # fig.savefig(f'SliceFSPM-fig-{图}.pdf', dpi=400, bbox_inches='tight', pad_inches=0)
d['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'] = 500



plt.show()
