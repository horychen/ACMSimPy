# %%
from rich import print

DISABLE_STUPID_WARNING = True


# data frame example
class CustomDataFrame:
    def __init__(self) -> None:
        self.plot_details = []

    def load(self, path1, path2):
        with open(path1, 'r', encoding='utf-8') as f:
            user_figs = f.read()
        with open(path2, 'r', encoding='utf-8') as f:
            signals_library = f.read()
        user_figs = user_figs.split('\n')
        signals_library = signals_library.split('\n')
        try:
            for i in range(len(user_figs)):
                user_fig = user_figs[i].split(',')
                for j in range(len(user_fig)):
                    user_fig[j] = user_fig[j].strip()
                if user_fig[0] == '':
                    continue
                self.plot_details.append({'data_axis': user_fig[0],
                                          'data_signal': user_fig[1:],
                                          'data_signal_label': user_fig[1:],
                                          'data_signal_num': len(user_fig[1:])})
            for i in range(len(signals_library)):
                signal = signals_library[i].split(',')
                for j in range(len(signal)):
                    signal[j] = signal[j].strip()
                if signal[0] == '':
                    continue
                for k in range(len(self.plot_details)):
                    for l in range(len(self.plot_details[k]['data_signal'])):
                        if signal[0] == self.plot_details[k]['data_signal'][l]:
                            self.plot_details[k]['data_signal_label'][l] = signal[3]
        except:
            raise Exception('user_cjh.txt or signals_library is not in the correct format.')

    def generate_function(self):
        with open(os.path.dirname(__file__) + '/collect_data.py', 'w') as f:
            f.write(
                f'''import numpy as np\ndef collect_data(watch_data, watch_index, CTRL, ACM, reg_id, reg_iq, reg_speed, fe_htz):\n''')
            index = 0
            for i in range(len(self.plot_details)):
                for j in range(len(self.plot_details[i]['data_signal'])):
                    f.write(f'\twatch_data[{index}][watch_index] = {self.plot_details[i]["data_signal"][j]}\n')
                    index += 1
            f.write(f'''\twatch_index += 1\n\treturn watch_index''')

    def plot(self, machine_times, watch_data):
        plt.style.use('bmh')  # https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
        mpl.rc('font', family='Times New Roman', size=10.0)
        mpl.rc('legend', fontsize=16)
        mpl.rcParams['lines.linewidth'] = 0.75  # mpl.rc('lines', linewidth=4, linestyle='-.')
        mpl.rcParams['mathtext.fontset'] = 'stix'

        total = len(self.plot_details)
        index = 0
        figure_index = 0

        if total < 6:
            fig, axes = plt.subplots(nrows=total, ncols=1, dpi=150, facecolor='w', figsize=(8, 12), sharex=True)
        else:
            fig, axes_ = plt.subplots(nrows=(total + 1) // 2, ncols=2, dpi=150, facecolor='w', figsize=(8, 12),
                                      sharex=True)
            axes = np.ravel(axes_)

        for plot_detail in self.plot_details:
            ax = axes[figure_index]
            figure_index += 1
            for i in range(len(plot_detail['data_signal'])):
                ax.plot(machine_times, watch_data[index], label=plot_detail['data_signal_label'][i])
                index += 1
            ax.set_ylabel(plot_detail['data_axis'], multialignment='center')
            ax.legend(loc=1, fontsize=12)
            ax.grid(True)
        axes[-1].set_xlabel('Time [s]')
        # plt.show()
        # return
        return {self.plot_details[i]['data_signal'][j]: watch_data[i * len(self.plot_details[i]['data_signal']) + j] for
                i in range(len(self.plot_details)) for j in range(len(self.plot_details[i]['data_signal']))}


    def lissajou(self, watch_data_as_dict, period, path):
        with open(path, 'r', encoding='utf-8') as f:
            user_figs = f.read()
        user_figs = user_figs.split('\n')
        user_fig_config = []
        try:
            for i in range(len(user_figs)):
                user_fig = user_figs[i].split(',')
                for j in range(len(user_fig)):
                    user_fig[j] = user_fig[j].strip()
                if user_fig[0]:
                    user_fig_config.append(user_fig)
        except:
            raise Exception('user_yzz.txt is not in the correct format.')
        plt.style.use('bmh')  # https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
        mpl.rc('font', family='Times New Roman', size=10.0)
        mpl.rc('legend', fontsize=16)
        mpl.rcParams['lines.linewidth'] = 0.75  # mpl.rc('lines', linewidth=4, linestyle='-.')
        mpl.rcParams['mathtext.fontset'] = 'stix'


        total = len(user_fig_config)
        figure_index = 0

        if total < 6:
            fig, axes = plt.subplots(nrows=total, ncols=1, dpi=150, facecolor='w', figsize=(8, 12))
        else:
            fig, axes_ = plt.subplots(nrows=(total + 1) // 2, ncols=2, dpi=150, facecolor='w', figsize=(8, 12))
            axes = np.ravel(axes_)

        for i in range(len(user_fig_config)):
            ax = axes[figure_index]
            figure_index += 1
            x_lim_low = min(watch_data_as_dict[user_fig_config[i][0]][int(int(user_fig_config[i][4])/period): int(int(user_fig_config[i][5])/period)])
            x_lim_high = max(watch_data_as_dict[user_fig_config[i][0]][int(int(user_fig_config[i][4])/period): int(int(user_fig_config[i][5])/period)])
            x_lim_shift = (x_lim_high - x_lim_low) / 10
            ax.plot(watch_data_as_dict[user_fig_config[i][0]][int(int(user_fig_config[i][4])/period): int(int(user_fig_config[i][5])/period)],
                    watch_data_as_dict[user_fig_config[i][1]][int(int(user_fig_config[i][4])/period): int(int(user_fig_config[i][5])/period)])
            ax.set_xlabel(user_fig_config[i][2], multialignment='center')
            ax.set_ylabel(user_fig_config[i][3], multialignment='center')
            ax.set_xlim(x_lim_low - x_lim_shift, x_lim_high + x_lim_shift)
            ax.set_ylim(x_lim_low - x_lim_shift, x_lim_high + x_lim_shift)
            ax.grid(True)
        plt.show()
        return None



custom = CustomDataFrame()
import os
custom.load(os.path.dirname(__file__) + '/user_cjh.txt', 
            os.path.dirname(__file__) + '/signals_library.txt')
custom.generate_function()

d = d_user_input_motor_dict = {
    # Timing
    'CL_TS': 1e-4,
    'TIME_SLICE': 10,
    'NUMBER_OF_SLICES': 1,
    'VL_EXE_PER_CL_EXE': 5,
    'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
    'CTRL.bool_apply_speed_closed_loop_control': True,
    'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
    'CTRL.bool_apply_sweeping_frequency_excitation': False,
    'CTRL.bool_overwrite_speed_commands': True,
    'CTRL.bool_zero_id_control': True,
    'FOC_delta': 10,  # 25, # 6.5
    'FOC_desired_VLBW_HZ': 120,  # 60
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
# 小电感电机
d['init_npp'] = 4
d['init_IN'] = 3
d['init_R'] = 1.10
d['init_Ld'] = 0.00496
d['init_Lq'] = 0.00496
d['init_KE'] = 0.1
d['init_KA'] = 0.1
d['init_Rreq'] = 0.0
d['init_Js'] = 0.000617
d['DC_BUS_VOLTAGE'] = 110

from tutorials_ep9_flux_estimator import *

print('Simulation_Benchmark')
# Auto-tuning PI
if d['CL_SERIES_KP'] is None:
    # sys.path.append(os.path.join(os.path.dirname(__file__), "tuner"))
    import tuner

    tuner.tunner_wrapper(d)
    print('\tAuto tuning...')
    print(f'\t{d=}\n')
else:
    print('\tSkip tuning.')

CTRL = The_Motor_Controller(CL_TS=d['CL_TS'],
                            VL_TS=d['VL_EXE_PER_CL_EXE'] * d['CL_TS'],
                            init_npp=d['init_npp'],
                            init_IN=d['init_IN'],
                            init_R=d['init_R'],
                            init_Ld=d['init_Ld'],
                            init_Lq=d['init_Lq'],
                            init_KE=d['init_KE'],
                            init_Rreq=d['init_Rreq'],
                            init_Js=d['init_Js'],
                            DC_BUS_VOLTAGE=d['DC_BUS_VOLTAGE'])
CTRL.bool_apply_decoupling_voltages_to_current_regulation = d[
    'CTRL.bool_apply_decoupling_voltages_to_current_regulation']
CTRL.bool_apply_sweeping_frequency_excitation = d['CTRL.bool_apply_sweeping_frequency_excitation']
CTRL.bool_yanzhengzhang = False  ####
# CTRL.bool_overwrite_speed_commands = d['CTRL.bool_overwrite_speed_commands']
CTRL.bool_zero_id_control = d['CTRL.bool_zero_id_control']
ACM = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=d['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])

fe_htz = Variables_FluxEstimator_Holtz03(CTRL.R)


reg_dispX = The_PID_Regulator(d['disp.Kp'], d['disp.Ki'], d['disp.Kd'], d['disp.tau'], d['disp.OutLimit'],
                              d['disp.IntLimit'], d['CL_TS'])
reg_dispY = The_PID_Regulator(d['disp.Kp'], d['disp.Ki'], d['disp.Kd'], d['disp.tau'], d['disp.OutLimit'],
                              d['disp.IntLimit'], d['CL_TS'])

if False:
    # Use incremental_pi codes
    reg_id = The_PI_Regulator(d['CL_SERIES_KP'], d['CL_SERIES_KP'] * d['CL_SERIES_KI'] * CTRL.CL_TS,
                              d['DC_BUS_VOLTAGE'] / 1.732)  # 我们假设调制方式是SVPWM，所以母线电压就是输出电压的线电压最大值，而我们用的是恒相幅值变换，所以限幅是相电压。
    reg_iq = The_PI_Regulator(d['CL_SERIES_KP'], d['CL_SERIES_KP'] * d['CL_SERIES_KI'] * CTRL.CL_TS,
                              d['DC_BUS_VOLTAGE'] / 1.732)  # 我们假设调制方式是SVPWM，所以母线电压就是输出电压的线电压最大值，而我们用的是恒相幅值变换，所以限幅是相电压。
    reg_speed = The_PI_Regulator(d['VL_SERIES_KP'], d['VL_SERIES_KP'] * d['VL_SERIES_KI'] * CTRL.VL_TS,
                                 d['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * d['init_IN'])  # IN 是线电流有效值，我们这边限幅是用的电流幅值。
else:
    # Use tustin_pi codes
    local_Kp = d['CL_SERIES_KP']
    if d['CTRL.bool_apply_decoupling_voltages_to_current_regulation'] == False:
        local_Ki = d['CL_SERIES_KP'] * d['CL_SERIES_KI'] * d[
            'FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False']
        print(
            '\tNote bool_apply_decoupling_voltages_to_current_regulation is False, to improve the current regulator performance a factor of %g has been multiplied to CL KI.' % (
                d['FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False']))
    else:
        local_Ki = d['CL_SERIES_KP'] * d['CL_SERIES_KI']
    local_Kd = 0.0
    local_tau = 0.0
    local_OutLimit = d['DC_BUS_VOLTAGE'] / 1.732
    local_IntLimit = 1.0 * d[
        'DC_BUS_VOLTAGE'] / 1.732  # Integrator having a lower output limit makes no sense. For example, the q-axis current regulator needs to cancel back emf using the integrator output for almost full dc bus voltage at maximum speed.
    reg_id = The_PID_Regulator(local_Kp, local_Ki, local_Kd, local_tau, local_OutLimit, local_IntLimit, d['CL_TS'])
    reg_iq = The_PID_Regulator(local_Kp, local_Ki, local_Kd, local_tau, local_OutLimit, local_IntLimit, d['CL_TS'])
    print(f'\t{reg_id.OutLimit=} V')

    local_Kp = d['VL_SERIES_KP']
    local_Ki = d['VL_SERIES_KP'] * d['VL_SERIES_KI']
    local_Kd = 0.0
    local_tau = 0.0
    local_OutLimit = d['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * d['init_IN']
    local_IntLimit = 1.0 * d['VL_LIMIT_OVERLOAD_FACTOR'] * 1.414 * d['init_IN']
    reg_speed = The_PID_Regulator(local_Kp, local_Ki, local_Kd, local_tau, local_OutLimit, local_IntLimit, CTRL.VL_TS)
    print(f'\t{reg_speed.OutLimit=} A')

CTRL, ACM, reg_id, reg_iq, reg_speed, reg_dispX, reg_dispY, fe_htz

# simulate to generate NUMBER_OF_SLICES*TIME_SLICE sec of data
for ii in range(d['NUMBER_OF_SLICES']):
    # perform animation step
    machine_times, watch_data = ACMSimPyIncremental(t0=ii * d['TIME_SLICE'], TIME=d['TIME_SLICE'],
                                                    ACM=ACM,
                                                    CTRL=CTRL,
                                                    reg_id=reg_id,
                                                    reg_iq=reg_iq,
                                                    reg_speed=reg_speed,
                                                    fe_htz=fe_htz)

# TODO:  程序员大哥，给我个好字典，谢谢您了！
watch_data_as_dict = custom.plot(machine_times, watch_data)
custom.lissajou(watch_data_as_dict, d['CL_TS'], os.path.dirname(__file__) + '/user_yzz.txt')
# Lissajour plot
# plt.plot(watch_data_as_dict['fe_htz.psi_2[0]'], watch_data_as_dict['fe_htz.psi_2[1]'])
# plt.show()












quit()
# which algorithm for torque and speed control
# list_of_control_strategies = {
#     "classic": ["ifoc", "dfoc"],
#     "advanced": ["sensorless ifoc",]
# }#1f INCREMENTAL PID
# print(list_of_control_strategies)

# which algorithm for flux estimation
# list_of_flux_estimator_algorithms = {
#     "voltage model": [
#         {
#             "saturation function based saturation time difference correction": 'fe01_satime',
#         }
#     ],
#     "full order model":[],
# }
# print(list_of_flux_estimator_algorithms)

# which algorithm for speed estimation

# which signals to look at
# parse : 

# fe_htz.psi_2[1], '[Wb]', r'$\psi_\beta$', dataframe, 

# votlage, '[V]', r'$\psi_\beta$', dataframe, 

# eval
# exec

# specify motor parameters and tune contol coefficients if any

# specify working conditions if any

图 = 1  # 空载加速、加载、反转
# 小电感电机
d['init_npp'] = 22
d['init_IN'] = 1.3 * 6 / 1.414
d['init_R'] = 0.035
d['init_Ld'] = 1 * 0.036 * 1e-3
d['init_Lq'] = 1 * 0.036 * 1e-3
d['init_KE'] = 0.0125
d['init_KA'] = 0.0125
d['init_Rreq'] = 0.0
d['init_Js'] = 0.44 * 1e-4
d['DC_BUS_VOLTAGE'] = 10
d[
    'user_system_input_code'] = '''if ii < 1: CTRL.cmd_idq[0] = 0.0; CTRL.cmd_rpm = 150 \nelif ii <5: ACM.TLoad = 0.2 \nelif ii <100: CTRL.cmd_rpm = -150'''


def 图1画图代码():
    plt.style.use('bmh')  # https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
    mpl.rc('font', family='Times New Roman', size=10.0)
    mpl.rc('legend', fontsize=10)
    mpl.rcParams['lines.linewidth'] = 0.75  # mpl.rc('lines', linewidth=4, linestyle='-.')
    mpl.rcParams['mathtext.fontset'] = 'stix'

    fig, axes = plt.subplots(nrows=7, ncols=1, dpi=150, facecolor='w', figsize=(8, 12), sharex=True)

    ax = axes[0]
    ax.plot(global_machine_times, gdd['CTRL.cmd_rpm'], label=r'$\omega_r^*$')
    ax.plot(global_machine_times, gdd['CTRL.omega_r_mech'], label=r'$\omega_r$')
    ax.set_ylabel(r'Speed [r/min]', multialignment='center')  # ) #, fontdict=font)
    # ax.legend(loc=2, prop={'size': 6})
    ax.legend(loc=1, fontsize=6)

    ax = axes[1]
    ax.plot(global_machine_times, gdd['CTRL.cmd_idq[0]'], label=r'$i_d^*$')
    ax.plot(global_machine_times, gdd['CTRL.idq[0]'], label=r'$i_d$')
    # ax.plot(global_machine_times, gdd['ACM.iD'])
    ax.set_ylabel(r'$i_d$ [A]', multialignment='center')  # , fontdict=font)

    ax = axes[2]
    ax.plot(global_machine_times, gdd['CTRL.cmd_idq[1]'], label=r'$i_q^*$')
    ax.plot(global_machine_times, gdd['CTRL.idq[1]'], label=r'$i_q$')
    ax.set_ylabel(r'$i_q$ [A]', multialignment='center')  # , fontdict=font)

    ax = axes[3]
    ax.plot(global_machine_times, gdd['ACM.Tem'], label=r'ACM.$T_{\rm em}$')
    ax.plot(global_machine_times, gdd['CTRL.Tem'], label=r'CTRL.$T_{\rm em}$')
    ax.set_ylabel(r'$T_{\rm em}$ [Nm]', multialignment='center')  # , fontdict=font)

    ax = axes[4]
    ax.plot(global_machine_times, (gdd['ACM.udq[0]']), label=r'$u_d$')  # lpf1_inverter
    ax.plot(global_machine_times, gdd['CTRL.cmd_udq[0]'], label=r'$u_d^*$')
    ax.set_ylabel(r'$u_d$ [V]', multialignment='center')  # , fontdict=font)

    ax = axes[5]
    # ax.plot(global_machine_times, (gdd['ACM.udq[1]']), label=r'$u_q$') # lpf1_inverter
    # ax.plot(global_machine_times, gdd['CTRL.cmd_udq[1]'], label=r'$u_q^*$')
    ax.plot(global_machine_times, gdd['fe_htz.u_offset[0]'], label=r'$u_{{\rm offset},\alpha}$')
    ax.plot(global_machine_times, gdd['fe_htz.u_offset[1]'], label=r'$u_{{\rm offset},\beta}$')
    ax.set_ylabel(r'u offset [V]', multialignment='center')  # , fontdict=font)

    ax = axes[6]
    # ax.plot(global_machine_times, gdd['CTRL.cmd_uab[0]'], label=r'$u_\alpha$')
    # ax.plot(global_machine_times, gdd['CTRL.cmd_uab[1]'], label=r'$u_\beta$')
    # ax.set_ylabel(r'$u_{\alpha,\beta}$ [V]', multialignment='center') #, fontdict=font)
    ax.plot(global_machine_times, gdd['fe_htz.psi_2[0]'], label=r'$\psi_\alpha$')
    ax.plot(global_machine_times, gdd['fe_htz.psi_2[1]'], label=r'$\psi_\beta$')
    ax.set_ylabel(r'$\psi_2$ [Wb]', multialignment='center')  # , fontdict=font)

    for ax in axes:
        ax.grid(True)
        ax.legend(loc=1)
        # for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        #     tick.label.set_font(font)
    axes[-1].set_xlabel('Time [s]')  # , fontdict=font)
    return fig
