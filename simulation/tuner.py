# -*- coding: utf-8 -*-
from pylab import np, plt, mpl
import control

def get_coeffs_dc_motor_current_regulator(R, L, Bandwidth_Hz):
    Kp = Bandwidth_Hz * 2 * np.pi * L
    Ki = R / L
    # print('DEBUG', Ki, R, L)
    return Kp, Ki

def get_coeffs_dc_motor_SPEED_regulator(J_s, n_pp, KA, delta, currentBandwidth_radPerSec):
    speedKi = currentBandwidth_radPerSec / delta**2
    speedKp = J_s/n_pp / (1.5*n_pp*KA) * delta * speedKi
    return speedKp, speedKi

def 逆上位机速度PI系数转换CODE(iSMC_speedKp, iSMC_speedKiCode, VL_TS, J_s):
    上位机速度KP = iSMC_speedKp / (2*np.pi * J_s)
    iSMC_speedKi = iSMC_speedKiCode/iSMC_speedKp/VL_TS
    上位机速度KI = iSMC_speedKi / (上位机速度KP * np.pi/2) * 1000
    return 上位机速度KP, 上位机速度KI

def 上位机速度PI系数转换CODE(上位机速度KP, 上位机速度KI, VL_TS, J_s):
    # 源代码
        # Kvp_in = 上位机速度KP
        # Kvi_in = 上位机速度KI/1000 * kvp_in * np.pi/2
        # iSMC_speedKp = kvp_in * 2*np.pi * J_s
        # iSMC_speedKi = iSMC_speedKp * kvi_in * VL_TS
    # 简洁点：
    iSMC_speedKp = 上位机速度KP * 2*np.pi * J_s
    iSMC_speedKi = 上位机速度KI/1000 * 上位机速度KP * np.pi/2
    iSMC_speedKiCode = iSMC_speedKi * iSMC_speedKp  * VL_TS
    return iSMC_speedKp, iSMC_speedKi, iSMC_speedKiCode

# current reference to current measurement
def c2c_design(R, L, CLBW_Hz=1000, CL_TS=1/20e3):
    currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R, L, CLBW_Hz)
    currentKiCode = currentKi * currentKp * CL_TS
    if True:
        # 这里打印的用于实验中CCS的debug窗口检查电流环PI系数
        上位机电流KP = CLBW_Hz
        上位机电流KI = 1000
        iSMC_currentKp = 上位机电流KP * L * 2*np.pi
        iSMC_currentKi = 上位机电流KI/1000 * R/L
        iSMC_currentKiCode = iSMC_currentKi * CL_TS * iSMC_currentKp
        print(f'\tiSMC_currentKp={iSMC_currentKp:g}, \
                  iSMC_currentKi={iSMC_currentKi:g}, \
                  iSMC_currentKiCode={iSMC_currentKiCode:g}')
        print(f'\tSimC_currentKp={currentKp:g}, \
                  SimC_currentKi={currentKi:g}, \
                  SimC_currentKiCode={currentKiCode:g}')
        print(f'\t上位机电流KP={上位机电流KP:g}, \
                  上位机电流KI={上位机电流KI:g}')
    Gi_closed = control.tf([1], [L/currentKp, 1]) # current loop zero-pole cancelled already
    currentBandwidth_radPerSec = currentKp/L

    # 注意，我们研究的开环传递函数是以电流给定为输入的，而不是以转速控制误差为输入，这样仿真和实验容易实现一点。
    c2c_tf = Gi_closed

    # fig5 = plt.figure(fignum)
    # plt.title('Designed Current Ref. to Velocity Meas. Transfer Function')
    mag, phase, omega = control.bode_plot(c2c_tf, 2*np.pi*np.logspace(0,4,500), dB=1, Hz=1, deg=1, lw='0.5', label=f'{CLBW_Hz:g} Hz')
    open_cutoff_frequency_HZ = omega[(np.abs(mag-0.0)).argmin()]/2/np.pi
    # print('\tCut-off frequency (without speed PI regulator):', open_cutoff_frequency_HZ, 'Hz')
    return  (currentKp, currentKi), \
            (上位机电流KP, 上位机电流KI), \
            (mag, phase, omega)

# current reference to velocity measaurement (this is not velocity open loop, because speed PI is not considered)
def c2v_design(R, L, n_pp, J_s, KA, B=0, CLBW_Hz=1000, CL_TS=1/20e3, fignum=5):

    currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R, L, CLBW_Hz)
    currentKiCode = currentKi * currentKp * CL_TS
    if True:
        # 这里打印的用于实验中CCS的debug窗口检查电流环PI系数
        上位机电流KP = CLBW_Hz
        上位机电流KI = 1000
        iSMC_currentKp = 上位机电流KP * L * 2*np.pi
        iSMC_currentKi = 上位机电流KI/1000 * R/L
        iSMC_currentKiCode = iSMC_currentKi * CL_TS * iSMC_currentKp
        print(f'\tiSMC_currentKp={iSMC_currentKp:g}, \
                  iSMC_currentKi={iSMC_currentKi:g}, \
                  iSMC_currentKiCode={iSMC_currentKiCode:g}')
        print(f'\tSimC_currentKp={currentKp:g}, \
                  SimC_currentKi={currentKi:g}, \
                  SimC_currentKiCode={currentKiCode:g}')
        print(f'\t上位机电流KP={上位机电流KP:g}, \
                  上位机电流KI={上位机电流KI:g}')
    Gi_closed = control.tf([1], [L/currentKp, 1]) # current loop zero-pole cancelled already
    currentBandwidth_radPerSec = currentKp/L

    KT = 1.5*n_pp*KA
    dc_motor_motion = control.tf([KT], [J_s/n_pp, B]) # [Apk] to [elec.rad/s]
    print(dc_motor_motion)
    # quit()

    # 注意，我们研究的开环传递函数是以电流给定为输入的，而不是以转速控制误差为输入，这样仿真和实验容易实现一点。
    # Gw_open = dc_motor_motion * Gi_closed * speedPI
    c2v_tf = dc_motor_motion * Gi_closed

    # fig5 = plt.figure(fignum)
    # plt.title('Designed Current Ref. to Velocity Meas. Transfer Function')
    mag, phase, omega = control.bode_plot(c2v_tf, 2*np.pi*np.logspace(0,4,500), dB=1, Hz=1, deg=1, lw='0.5', label=f'{CLBW_Hz:g} Hz')
    open_cutoff_frequency_HZ = omega[(np.abs(mag-0.0)).argmin()]/2/np.pi
    # print('\tCut-off frequency (without speed PI regulator):', open_cutoff_frequency_HZ, 'Hz')
    return  (currentKp, currentKi), \
            (上位机电流KP, 上位机电流KI), \
            (mag, phase, omega)

# velocity reference to velocity measaurement
def iterate_for_desired_bandwidth( delta, desired_VLBW_Hz, motor_dict, CLBW_Hz_initial=1000, CLBW_Hz_stepSize=100):

    R          = motor_dict['Rs']
    L          = motor_dict['Ls']
    J_s        = motor_dict['J_s']
    JLoadRatio = motor_dict['JLoadRatio']
    n_pp       = motor_dict['n_pp']
    KE         = motor_dict['KE']
    if KE != 0:
        KA = KE
    else:
        KA = motor_dict['KA']
    CL_TS      = motor_dict['CL_TS']
    VL_TS      = motor_dict['VL_TS']
    J_total    = J_s*(1+JLoadRatio) 

    CLBW_Hz = CLBW_Hz_initial  #100 # Hz (initial)
    VLBW_Hz = -10  # Hz (initial)
    count = 0

    while True:
        count += 1
        if count>20:
            msg = f'Loop count 20 is reached. Step size is {CLBW_Hz_stepSize} Hz.'
            print(msg)
            # raise Exception()
            break

        # Current loop (Tune its bandwidth to support required speed response)
        if abs(VLBW_Hz - desired_VLBW_Hz)<=10: # Hz
            break
        else:
            if VLBW_Hz > desired_VLBW_Hz:
                CLBW_Hz -= CLBW_Hz_stepSize # Hz
                if CLBW_Hz<=0:
                    raise Exception(f'Negative CLBW_Hz. Maybe change the step size of "CLBW_Hz" ({CLBW_Hz_stepSize} Hz) and try again.')
                    break
            else:
                CLBW_Hz += CLBW_Hz_stepSize # Hz
        # print(f'CLBW_Hz = {CLBW_Hz}')

        currentKp, currentKi = get_coeffs_dc_motor_current_regulator(R, L, CLBW_Hz)
        currentKiCode = currentKi * currentKp * CL_TS
        if True:
            # 这里打印的用于实验中CCS的debug窗口检查电流环PI系数
            上位机电流KP = CLBW_Hz
            上位机电流KI = 1000
            iSMC_currentKp = 上位机电流KP * L * 2*np.pi
            iSMC_currentKi = 上位机电流KI/1000 * R/L
            iSMC_currentKiCode = iSMC_currentKi * CL_TS * iSMC_currentKp

        Gi_closed = control.tf([1], [L/currentKp, 1]) # current loop zero-pole cancelled already
        currentBandwidth_radPerSec = currentKp/L

        # Speed loop
        KT = 1.5*n_pp*KA
        dc_motor_motion = control.tf([KT*n_pp/J_total], [1, 0])
        speedKp, speedKi = get_coeffs_dc_motor_SPEED_regulator(J_total, n_pp, KA, delta, currentBandwidth_radPerSec)
        speedKiCode = speedKi * speedKp * VL_TS
        if True:
            # 这里打印的用于实验中TI的debug窗口检查系数

            上位机速度KP, 上位机速度KI = 逆上位机速度PI系数转换CODE(speedKp, speedKiCode, VL_TS, J_total)

            iSMC_speedKp, iSMC_speedKi, iSMC_speedKiCode = 上位机速度PI系数转换CODE(上位机速度KP, 上位机速度KI, VL_TS, J_total)

        speedPI = control.tf([speedKp, speedKp*speedKi], [1, 0])
        Gw_open = dc_motor_motion * Gi_closed * speedPI

        # C2C
        # fig = plt.figure(101)
        c2c_tf = Gi_closed
        mag, phase, omega = control.bode_plot(c2c_tf, 2*np.pi*np.logspace(0,4,500), dB=1, Hz=1, deg=1, lw='0.5', label=f'{CLBW_Hz:g} Hz', plot=False)
        CLBW_Hz = omega[(np.abs(mag-0.707)).argmin()]/2/np.pi
        C2C_designedMagPhaseOmega = mag, phase, omega
        # plt.close(fig)

        # C2V
        # fig = plt.figure(102)
        c2v_tf = dc_motor_motion * Gi_closed
        mag, phase, omega = control.bode_plot(c2v_tf, 2*np.pi*np.logspace(0,4,500), dB=1, Hz=1, deg=1, lw='0.5', label=f'{CLBW_Hz:g} Hz', plot=False)
        open_cutoff_frequency_HZ = omega[(np.abs(mag-0.0)).argmin()]/2/np.pi
        C2V_designedMagPhaseOmega = mag, phase, omega
        # plt.close(fig)

        # V2V
        # fig = plt.figure(103)
        Gw_closed = Gw_open / (1+Gw_open)
        mag, phase, omega = control.bode_plot(Gw_closed, 2*np.pi*np.logspace(0,4,500), dB=1, Hz=1, deg=1, lw='0.5', label=f'{delta:g}', plot=False)
        VLBW_Hz = omega[(np.abs(mag-0.707)).argmin()]/2/np.pi
        V2V_designedMagPhaseOmega = mag, phase, omega
        # plt.close(fig)

        # print(Gw_closed)
        # print('\tSpeed loop bandwidth:', VLBW_Hz, 'Hz')

    return  (currentKp, currentKi), \
            (speedKp, speedKi), \
            (上位机电流KP, 上位机电流KI), \
            (上位机速度KP, 上位机速度KI), \
            (C2C_designedMagPhaseOmega, C2V_designedMagPhaseOmega, V2V_designedMagPhaseOmega), \
            (CLBW_Hz, VLBW_Hz, open_cutoff_frequency_HZ)


def tunner_wrapper(d):
    motor_dict = dict()
    motor_dict['Rs'] =         d['init_R']
    motor_dict['Ls'] =         d['init_Lq']
    motor_dict['J_s'] =        d['init_Js']
    motor_dict['JLoadRatio'] = 0.0
    motor_dict['n_pp'] =       d['init_npp']
    motor_dict['KE'] =         d['init_KE']
    motor_dict['KA'] =         d['init_KA']
    motor_dict['CL_TS'] =      d['CL_TS']
    motor_dict['VL_TS'] =      d['CL_TS'] * d['VL_EXE_PER_CL_EXE']
    delta = d['FOC_delta']
    desired_VLBW_HZ = d['FOC_desired_VLBW_HZ']

    currentPI, speedPI, 上位机电流PI, 上位机速度PI, tuple_designedMagPhaseOmega, BW_in_Hz = iterate_for_desired_bandwidth(delta, desired_VLBW_HZ, motor_dict)
    # print(currentPI, speedPI, 上位机电流PI, 上位机速度PI, BW_in_Hz)

    d['CL_SERIES_KP'] = currentPI[0]
    d['CL_SERIES_KI'] = currentPI[1]
    d['VL_SERIES_KP'] = speedPI[0]
    d['VL_SERIES_KI'] = speedPI[1]
    d['CLBW_HZ'] = BW_in_Hz[0]
    d['VLBW_HZ'] = BW_in_Hz[1]
    return d

if __name__ == '__main__':

    d = d_user_input = {
            'CL_TS': 1e-4,
            'TIME_SLICE': 0.1,
            'NUMBER_OF_SLICES': 6,
            'VL_EXE_PER_CL_EXE': 5,
            'init_npp': 21,
            'init_IN': 72/1.414,
            'init_R': 0.1222,
            'init_Ld': 0.0011, # 0.007834, # 0.000502, # 4*0.000502, # 
            'init_Lq': 0.00125, # 0.0089, # 0.000571, # 4*0.000571, # 
            'init_KE': 0.127, # 150 / 1.732 / (450/60*6.28*21), # 285 / 1.5 / 72 / 21,
            # 'init_KE': 0.087559479,
            'init_Rreq': 0.0,
            'init_Js': 0.203,
            'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
            'CTRL.bool_apply_sweeping_frequency_excitation': False,
            'CTRL.bool_overwrite_speed_commands': True,
            'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,  # 500,
            'DC_BUS_VOLTAGE': 150*1.732,
            'FOC_delta': 6.5,
            'FOC_desired_VLBW_HZ': 40,
            'CL_SERIES_KP': None, # 1.61377, # 11.49, # [BW=207.316Hz] # 6.00116,  # [BW=106.6Hz]
            'CL_SERIES_KI': None, # 97.76,
            'VL_SERIES_KP': None, # 0.479932, # 0.445651, # [BW=38.6522Hz] # 0.250665 # [BW=22.2Hz]
            'VL_SERIES_KI': None, # 30.5565,
            'VL_LIMIT_OVERLOAD_FACTOR': 1,
        }
    d = tunner_wrapper(d)

    print(d)
