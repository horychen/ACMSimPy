# %%
############################################# PACKAGES
from numba.experimental import jitclass
from numba import njit, int32, float64
from pylab import np, plt, mpl
plt.style.use('ggplot')
import collect_data
import humans_give_commands

NS_GLOBAL = 6

############################################# CLASS DEFINITION 
class The_Motor_Controller:
    def __init__(self, 
        ELL_param = 0.019,
        **kwargs
    ):
        ''' MOTOR '''
        self.npp  = kwargs.get('init_npp', 26)
        self.IN   = kwargs.get('init_IN', 17)
        self.R    = kwargs.get('init_R', 0.12)
        self.Ld   = kwargs.get('init_Ld', 0.00046)
        self.Lq   = kwargs.get('init_Lq', 0.00056)
        self.KE   = kwargs.get('init_KE', 0.019)
        self.Rreq = kwargs.get('init_Rreq', 0)
        self.Js   = kwargs.get('init_Js', 0.000364)
        self.DC_BUS_VOLTAGE = kwargs.get('DC_BUS_VOLTAGE', 48)
        self.Js_inv = 1 / self.Js
        self.Lq_inv = 1 / self.Lq
        ''' CONTROL '''
        # constants
        self.CL_TS = kwargs.get('CL_TS', 1e-4)
        self.VL_TS = kwargs.get('VL_TS', 1e-4)
        self.velocity_loop_ceiling = self.VL_TS / self.CL_TS
        self.velocity_loop_counter = self.velocity_loop_ceiling - 1
        print('\tCTRL.velocity_loop_ceiling =', self.velocity_loop_ceiling)
        # feedback / input
        self.theta_d = kwargs.get('theta_d', 0.0)
        self. thetaerror = kwargs.get('thetaerror', 0.0)
        self.omega_r_elec = kwargs.get('omega_r_elec', 0.0)
        self.omega_syn = kwargs.get('omega_syn', 0.0)
        self.omega_slip = kwargs.get('omega_slip', 0.0)
        self.uab      = kwargs.get('uab', np.zeros(2, dtype=np.float64))
        # self.uab_prev = np.zeros(2, dtype=np.float64)
        # self.uab_curr = np.zeros(2, dtype=np.float64)
        self.iab      = kwargs.get('iab', np.zeros(2, dtype=np.float64))
        self.iab_prev = kwargs.get('iab', np.zeros(2, dtype=np.float64))
        self.iab_curr = kwargs.get('iab', np.zeros(2, dtype=np.float64))
        # states
        self.timebase = kwargs.get('timebase', 0.0)
        self.KA = kwargs.get('init_KE', 0.0)
        self.Tem = kwargs.get('Tem', 0.0)
        self.cosT = kwargs.get('cosT', 1.0)
        self.sinT = kwargs.get('sinT', 0.0)
        # commands
        self.cmd_idq = kwargs.get('cmd_idq', np.zeros(2, dtype=np.float64))
        self.cmd_udq = kwargs.get('cmd_udq', np.zeros(2, dtype=np.float64)) 
        self.cmd_uab = kwargs.get('cmd_uab', np.zeros(2, dtype=np.float64))
        self.cmd_rpm = kwargs.get('cmd_rpm', 0.0)
        if self.Rreq >0:
            self.cmd_psi = 0.9 # [Wb]
        else:
            self.cmd_psi = kwargs.get('init_KE', 0.1)# [Wb]
        self.index_voltage_model_flux_estimation = kwargs.get('CTRL.index_voltage_model_flux_estimation', 1)
        self.index_separate_speed_estimation = kwargs.get('CTRL.index_separate_speed_estimation', 0)
        self.use_disturbance_feedforward_rejection = kwargs.get('use_disturbance_feedforward_rejection', 0)
        self.bool_apply_decoupling_voltages_to_current_regulation = kwargs.get('CTRL.bool_apply_decoupling_voltages_to_current_regulation', True)
        self.bool_apply_speed_closed_loop_control = kwargs.get('CTRL.bool_apply_speed_closed_loop_control', True)
        self.bool_zero_id_control = kwargs.get('CTRL.bool_zero_id_control', True)
        self.bool_reverse_rotation = kwargs.get('CTRL.bool_reverse_rotation', True)
        self.flag_reverse_rotation = kwargs.get('flag_reverse_rotation', True)
        self.counter_rotation = kwargs.get('counter_rotation', 10)
        # sweep frequency
        self.bool_apply_sweeping_frequency_excitation = kwargs.get('CTRL.bool_apply_sweeping_frequency_excitation', False)
        self.bool_overwrite_speed_commands = kwargs.get('CTRL.bool_overwrite_speed_commands', True)
        self.CMD_CURRENT_SINE_AMPERE = kwargs.get('CMD_CURRENT_SINE_AMPERE', 1) # [A]
        self.CMD_SPEED_SINE_RPM = kwargs.get('CMD_SPEED_SINE_RPM', 100) # [rpm]
        self.CMD_SPEED_SINE_HZ = kwargs.get('CMD_SPEED_SINE_HZ', 0) # [Hz]
        self.CMD_SPEED_SINE_STEP_SIZE = kwargs.get('CMD_SPEED_SINE_STEP_SIZE', 1) # [Hz]
        self.CMD_SPEED_SINE_LAST_END_TIME = kwargs.get('CMD_SPEED_SINE_LAST_END_TIME', 0.0)
        self.CMD_SPEED_SINE_END_TIME = kwargs.get('CMD_SPEED_SINE_END_TIME', 0.0)
        self.CMD_SPEED_SINE_HZ_CEILING = kwargs.get('CMD_SPEED_SINE_HZ_CEILING', 0.0)
        ''' tools for error calculation '''
        #psi error calculate
        self.bool_counter = kwargs.get('bool_counter', False)
        self.counter_psi = kwargs.get('psi_max', 0)
        self.psi_max = kwargs.get('counter_psi', 0)
        self.psi_min = kwargs.get('psi_min', 0)
        self.psi_sum = kwargs.get('psi_sum', 0)
        self.psi_avg = kwargs.get('psi_avg', 0)
        self.psi_max_fin = kwargs.get('psi_max_fin', 0)
        self.psi_min_fin = kwargs.get('psi_min_fin', 0)
        self.psi_avg_fin = kwargs.get('psi_avg_fin', 0)
        #theta error calculate
        self.bool_counter_theta_error = kwargs.get('bool_counter_theta_error', False)
        self.counter_theta_error = kwargs.get('counter_theta_error', 0)
        self.thetaerror_max = kwargs.get('thetaerror_max', 0)
        self.thetaerror_min = kwargs.get('thetaerror_min', 0)
        self.thetaerror_sum = kwargs.get('thetaerror_sum', 0)
        self.thetaerror_avg = kwargs.get('thetaerror_avg', 0)
        self.thetaerror_max_fin = kwargs.get('thetaerror_max_fin', 0)
        self.thetaerror_min_fin = kwargs.get('thetaerror_min_fin', 0)
        self.thetaerror_avg_fin = kwargs.get('thetaerror_avg_fin', 0)
        ''' OBSERVER '''
        # feedback / input
        self.idq = kwargs.get('idq', np.zeros(2, dtype=np.float64))
        # state
        # self.NS_SPEED  = 6 # = max(NS_SPEED, NS_FLUX)
        self.xSpeed    = kwargs.get('xSpeed', np.zeros(NS_GLOBAL, dtype=np.float64))
        self.xTorque   = kwargs.get('xTorque', np.zeros(NS_GLOBAL, dtype=np.float64))
        # outputs
        self.speed_observer_output_error = kwargs.get('speed_observer_output_error', 0.0)
        self.vartheta_d = kwargs.get('vartheta_d', 0.0)
        self.total_disrubance_feedforward = kwargs.get('total_disrubance_feedforward', 0.0)
        # gains
        omega_ob = 5000 # [rad/s]
        self.ell1 = kwargs.get('ell1', 0.0)
        self.ell2 = kwargs.get('ell2', 0.0)
        self.ell3 = kwargs.get('ell3', 0.0)
        self.ell4 = kwargs.get('ell4', 0.0)
        if False: # 2nd-order speed observer (assuming speed feedback)
            self.ell2 = 2 * omega_ob
            self.ell3 =     omega_ob**2 * self.Js/self.npp
        elif False: # 2nd-order position observer
            self.ell1 = 2 * omega_ob
            self.ell2 =     omega_ob**2 * self.Js/self.npp
        elif True: # 3rd-order position observer
            self.ell1 = 3 * omega_ob
            self.ell2 = 3 * omega_ob**2
            self.ell3 =     omega_ob**3 * self.Js/self.npp
        else: # 4th-order position observer
            self.ell1 = 4 * omega_ob
            self.ell2 = 6 * omega_ob**2
            self.ell3 = 4 * omega_ob**3 * self.Js/self.npp
            self.ell4 =     omega_ob**4

        self.one_over_six = kwargs.get('one_over_six', 1 / 6)
        self.use_encoder_angle_no_matter_what = kwargs.get('CTRL.use_encoder_angle_no_matter_what', True)
        self.flux_estimate_amplitude = kwargs.get('init_KE', 0.019)
        # mixed holtz and boldea to estimate K_E
        self.ell = kwargs.get('init_KE', 0.0)
        self.ell_prev = kwargs.get('init_KE', 0.0)
        self.ell = ELL_param
        self.ell_prev = ELL_param
        self.psi_A_amplitude = kwargs.get('psi_A_amplitude', 0.0)
        self.Kp_KE_estimate = kwargs.get('Kp_KE_estimate', -150)
        self.K1_for_max_ell = kwargs.get('K1_for_max_ell', 1000)
        self.K2_for_min_ell = kwargs.get('K2_for_min_ell', -0.005)
        # boldea 2008
        self.rotor_flux_error = kwargs.get('rotor_flux_error', np.zeros(NS_GLOBAL, dtype=np.float64))
        self.OFFSET_VOLTAGE_ALPHA = kwargs.get('OFFSET_VOLTAGE_ALPHA', 0.0)
        self.OFFSET_VOLTAGE_BETA = kwargs.get('OFFSET_VOLTAGE_BETA', 0.0)
        self.VM_PROPOSED_PI_CORRECTION_GAIN_P = kwargs.get('VM_PROPOSED_PI_CORRECTION_GAIN_P', 0.1)
        self.VM_PROPOSED_PI_CORRECTION_GAIN_I = kwargs.get('VM_PROPOSED_PI_CORRECTION_GAIN_I', 0.01)
        self.correction_integral_term = kwargs.get('correction_integral_term', np.zeros(NS_GLOBAL, dtype=np.float64))
        self.emf_stator = kwargs.get('emf_stator', np.zeros(NS_GLOBAL, dtype=np.float64))
        self.cmd_psi_mu = kwargs.get('cmd_psi_mu', np.zeros(NS_GLOBAL, dtype=np.float64))
        # natural speed observer
        self.nsoaf_KP = kwargs.get('nsoaf_KP', 0.0)
        self.nsoaf_KI = kwargs.get('nsoaf_KI', 0.0)
        self.nsoaf_KD = kwargs.get('nsoaf_KD', 0.0)
        self.nsoaf_omega_ob = kwargs.get('nsoaf_omega_ob', 0.0)
        self.nsoaf_set_omega_ob = kwargs.get('nsoaf_set_omega_ob', 200)
        self.TUNING_IGNORE_UQ  = kwargs.get('TUNING_IGNORE_UQ', False)
        self.tuning_nsoaf = kwargs.get('tuning_nsoaf', False)
        self.nsoaf_output_error = kwargs.get('nsoaf_output_error', 0.0)
        self.udq = kwargs.get('udq', np.zeros(2, dtype=np.float64))
        self.NSOAF_SPMSM_OR_IPMSM = kwargs.get('NSOAF_SPMSM_OR_IPMSM', 1)
        self.nsoaf_uQ = kwargs.get('nsoaf_uQ', 0.0)
        self.active_power_real = kwargs.get('active_power_real', 0.0)
        self.active_power_est = kwargs.get('active_power_est', 0.0)
        self.active_power_error = kwargs.get('active_power_error', 0.0)
        self.nsoaf_xTem = kwargs.get('nsoaf_xTem', 0.0)
        self.CLARKE_TRANS_TORQUE_GAIN = kwargs.get('CLARKE_TRANS_TORQUE_GAIN', 1.5)
        self.nsoaf_xIq      = kwargs.get('nsoaf_xIq', 0.0)
        self.nsoaf_xOmg     = kwargs.get('nsoaf_xOmg', 0.0)
        self.nsoaf_xTL      = kwargs.get('nsoaf_xTL', 0.0)
        self.nsoaf_xSpeed = kwargs.get('nsoaf_xSpeed', np.zeros(NS_GLOBAL, dtype=np.float64))
        self.uQ_now_filtered = kwargs.get('uQ_now_filtered', 0.0)
        self.idq_c = kwargs.get('idq_c', np.zeros(2, dtype=np.float64))
class The_AC_Machine:
    def __init__(self, CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1, ACM_param=1.0):
        # name plate data
        self.npp = CTRL.npp
        self.npp_inv = 1.0/self.npp
        self.IN  = CTRL.IN
        # electrical parameters
        self.R   = CTRL.R
        self.Ld  = CTRL.Ld
        self.Lq  = CTRL.Lq * ACM_param
        print(f'ACM: {self.Lq=}')

        self.KE  = CTRL.KE
        self.Rreq  = CTRL.Rreq
        # mechanical parameters
        self.Js  = CTRL.Js # kg.m^2

        # states
        self.NS = 5
        self.x = np.zeros(self.NS, dtype=np.float64)
        self.x[2] = CTRL.KA
        # inputs
        self.uab = np.zeros(2, dtype=np.float64)
        self.udq = np.zeros(2, dtype=np.float64)
        self.TLoad = 0
        # output
        self.omega_slip = 0.0
        self.omega_r_elec = 0.0
        self.omega_r_mech = 0.0
        self.omega_syn = 0.0
        self.theta_d = 0.0
        self.theta_d_mech = 0.0
        self.KA = CTRL.KA
        self.iD = 0.0
        self.iQ = 0.0
        self.iAlfa = 0.0
        self.iBeta = 0.0
        self.ia = 0.0
        self.ib = 0.0
        self.ic = 0.0
        self.Tem = 0.0
        self.cosT = 1.0
        self.sinT = 0.0
        self.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD = MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
        self.bool_apply_load_model = False

class The_PI_Regulator:
    def __init__(self, KP_CODE, KI_CODE, OUTPUT_LIMIT):
        self.Kp = KP_CODE
        self.Ki = KI_CODE
        self.Err      = 0.0
        self.setpoint = 0.0
        self.measurement = 0.0
        self.Out      = 0.0
        self.OutLimit = OUTPUT_LIMIT
        self.ErrPrev  = 0.0
        self.OutPrev  = 0.0

class The_PID_Regulator:
    def __init__(self, Kp, Ki, Kd, tau, OutLimit, IntLimit, T):

        # Regulator gains */
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd

        # Derivative low-pass filter time constant */
        self.tau = tau

        # Output limits */
        self.OutLimit = OutLimit

        # Integrator limits */
        self.IntLimit = IntLimit

        # Sample time (in seconds) */
        self.T = T

        # Regulator "memory" */
        self.integrator = 0.0
        self.prevError = 0.0            # Required for integrator */
        self.differentiator = 0.0
        self.prevMeasurement = 0.0      # Required for differentiator */

        # Regulator output */
        self.Out = 0.0

        # Regulator input */
        self.setpoint = 0.0
        self.measurement = 0.0

class SVgen_Object:
    def __init__(self, CPU_TICK_PER_SAMPLING_PERIOD):
        self.Ualfa = 0.0
        self.Ubeta = 0.0
        self.Unot = 0.0
        self.Ta = 0.5
        self.Tb = 0.5
        self.Tc = 0.5
        self.SYSTEM_MAX_PWM_DUTY_LIMATATION = 0.96
        self.SYSTEM_MIN_PWM_DUTY_LIMATATION = 0.04

        # Those variables are only needed in simulation
        self.bool_interupt_event = False
        self.bool_counting_down = False
        self.bool_RisingEdgeDelay_is_active  = np.zeros(3, dtype=np.float64)
        self.bool_FallingEdgeDelay_is_active = np.zeros(3, dtype=np.float64)
        self.carrier_counter = 0
        self.deadtime_counter = np.zeros(3, dtype=np.float64)
        self.S1, self.S2, self.S3, self.S4, self.S5, self.S6 = 0,0,0,0,0,0
        self.EPwm1Regs_CMPA_bit_CMPA = 0.5 * (0.5*CPU_TICK_PER_SAMPLING_PERIOD) # half of up/down counting maximum
        self.EPwm2Regs_CMPA_bit_CMPA = 0.5 * (0.5*CPU_TICK_PER_SAMPLING_PERIOD) # half of up/down counting maximum
        self.EPwm3Regs_CMPA_bit_CMPA = 0.5 * (0.5*CPU_TICK_PER_SAMPLING_PERIOD) # half of up/down counting maximum
        self.phase_U_gate_signal = 0
        self.phase_V_gate_signal = 0
        self.phase_W_gate_signal = 0
        self.voltage_potential_at_terminal = np.zeros(3, dtype=np.float64)
        self.line_to_line_voltage_AC = 0.0
        self.line_to_line_voltage_BC = 0.0
        self.line_to_line_voltage_AB = 0.0

class Variables_FluxEstimator_Holtz03:
    def __init__(self, IM_STAOTR_RESISTANCE, init_KE):

        self.xFlux = np.zeros(NS_GLOBAL, dtype=np.float64)
        self.xFlux[0] = init_KE
        self.xFlux[1] = 0
        self.xFlux[4] = init_KE
        self.psi_1 = np.zeros(2, dtype=np.float64)
        self.psi_2 = np.zeros(2, dtype=np.float64)
        self.psi_A = np.zeros(2, dtype=np.float64)
        self.psi_s = np.zeros(2, dtype=np.float64)
        self.psi_2_prev= np.zeros(2, dtype=np.float64)
        self.psi_e = 0
        self.psi_1_nonSat= np.zeros(2, dtype=np.float64)
        self.psi_2_nonSat= np.zeros(2, dtype=np.float64)

        self.psi_1_min= np.zeros(2, dtype=np.float64)
        self.psi_1_max= np.zeros(2, dtype=np.float64)
        self.psi_2_min= np.zeros(2, dtype=np.float64)
        self.psi_2_max= np.zeros(2, dtype=np.float64)

        self.rs_est   = IM_STAOTR_RESISTANCE
        # self.rreq_est = IM_ROTOR_RESISTANCE
        self.theta_d = 0.0
        self.Delta_t = 1
        self.u_offset= np.zeros(2, dtype=np.float64)

        self.u_off_original_lpf_input= np.zeros(2, dtype=np.float64) # holtz03 original (but I uses int32egrator instead of LPF)
        self.u_off_saturation_time_correction= np.zeros(2, dtype=np.float64) # exact offset calculation for compensation
        self.u_off_calculated_increment= np.zeros(2, dtype=np.float64)    # saturation time based correction
        self.GAIN_OFFSET_INIT = 10.0
        self.gain_off = self.GAIN_OFFSET_INIT  # HOLTZ_2002_GAIN_OFFSET; # 5; -> slow but stable // 50.1 // 20 -> too large then speed will oscillate during reversal near zero
        self.GAIN_OFFSET_REALTIME = 0.0

        self.flag_pos2negLevelA= np.zeros(2, dtype=np.int32)
        self.flag_pos2negLevelB= np.zeros(2, dtype=np.int32)
        self.time_pos2neg= np.zeros(2, dtype=np.float64)
        self.time_pos2neg_prev= np.zeros(2, dtype=np.float64)

        self.flag_neg2posLevelA= np.zeros(2, dtype=np.int32)
        self.flag_neg2posLevelB= np.zeros(2, dtype=np.int32)
        self.time_neg2pos= np.zeros(2, dtype=np.float64)
        self.time_neg2pos_prev= np.zeros(2, dtype=np.float64)

        self.psi_aster_max = 0.9 #IM_FLUX_COMMAND_DC_PART + IM_FLUX_COMMAND_SINE_PART

        self.maximum_of_sat_min_time= np.zeros(2, dtype=np.float64)
        self.maximum_of_sat_max_time= np.zeros(2, dtype=np.float64)
        self.sat_min_time= np.zeros(2, dtype=np.float64)
        self.sat_max_time= np.zeros(2, dtype=np.float64)
        self.sat_min_time_reg= np.zeros(2, dtype=np.float64)
        self.sat_max_time_reg= np.zeros(2, dtype=np.float64)
        self.extra_limit = 0.0
        self.flag_limit_too_low = False

        self.negative_cycle_in_count = np.zeros(2, dtype=np.float64)
        self.positive_cycle_in_count = np.zeros(2, dtype=np.float64)
        self.count_positive_in_one_cycle = np.zeros(2, dtype=np.float64)
        self.count_negative_in_one_cycle = np.zeros(2, dtype=np.float64)
        self.count_positive_cycle = 0
        self.count_negative_cycle = 0
        self.u_off_direct_calculated = np.zeros(2, dtype=np.float64)
        self.sign__u_off_saturation_time_correction = np.zeros(2, dtype=np.float64)
        self.sat_time_offset = np.zeros(2, dtype=np.float64)
        self.ell_compensation_flag = 0
        self.ell_compensation_flag_alpha = 0 
        self.ell_compensation_flag_beta = 0
        
        # no saturation time based 
        self.u_offset_correction_factor = 100
        self.u_offset_filered = np.zeros(2, dtype=np.float64)
############################################# OBSERVERS SECTION
def DYNAMICS_SpeedObserver(x, CTRL, SO_param=1.0):
    fx = np.zeros(NS_GLOBAL)

    # [rad]
    # output_error = np.sin(CTRL.theta_d - x[0])
    output_error = angle_diff(CTRL.theta_d, x[0]) # OE version 2
        # CTRL.output_error = np.sin(CTRL.theta_d - CTRL.xSpeed[0]) # OE version 1 simple and silly
        # CTRL.output_error = angle_diff(CTRL.theta_d - CTRL.xSpeed[0]) # OE version 2
        # CTRL.output_error = q-axis component # OE version 3 Boldea
    CTRL.speed_observer_output_error = output_error

    # 机械子系统 (omega_r_elec, theta_d, theta_r_mech)
    fx[0] = CTRL.ell1*output_error + x[1]
    fx[1] = CTRL.ell2*output_error + (CTRL.Tem + x[2]) * CTRL.npp/CTRL.Js # elec. angular rotor speed
    fx[2] = CTRL.ell3*output_error + x[3]
    fx[3] = CTRL.ell4*output_error + 0.0
    return fx

def DYNAMICS_FluxEstimator(x, CTRL, FE_param=1.0):
    fx = np.zeros(NS_GLOBAL)
    fx[0] = CTRL.uab[0] - CTRL.R * FE_param * CTRL.iab[0] - x[2] + 0.1
    fx[1] = CTRL.uab[1] - CTRL.R * FE_param * CTRL.iab[1] - x[3] 

    uhf_alfa = CTRL.cosT * (CTRL.ell - CTRL.flux_estimate_amplitude) * CTRL.Kp_KE_estimate
    uhf_beta = CTRL.sinT * (CTRL.ell - CTRL.flux_estimate_amplitude) * CTRL.Kp_KE_estimate

    fx[4] = CTRL.uab[0] - CTRL.R * FE_param * CTRL.iab[0] - (x[2] + uhf_alfa)
    fx[5] = CTRL.uab[1] - CTRL.R * FE_param * CTRL.iab[1] - (x[3] + uhf_beta)
    
    
    return fx
# nature speed oberver init
def init_nsoaf(CTRL):
    NSOAF_TL_P = 1
    NSOAF_TL_I =20
    NSOAF_TL_D = 0
    NSOAF_OMEGA_OBSERVER = 200

    CTRL.nsoaf_KP = NSOAF_TL_P
    CTRL.nsoaf_KI = NSOAF_TL_I
    CTRL.nsoaf_KD = NSOAF_TL_D
    CTRL.nsoaf_set_omega_ob = NSOAF_OMEGA_OBSERVER

    CTRL.KA = CTRL.KE + (CTRL.Ld - CTRL.Lq) * CTRL.idq[0]
    CTRL.nsoaf_omega_ob = CTRL.nsoaf_set_omega_ob
    nso_one_parameter_tuning(CTRL)
# nature speed oberver tuning
def nso_one_parameter_tuning(CTRL):
    if CTRL.nsoaf_omega_ob < 170:
        CTRL.nsoaf_set_omega_ob = 170
    one_over__npp_divided_by_Js__times__Lq_id_plus_KActive = 1 / (CTRL.npp * CTRL.Js_inv * CTRL.Lq_inv * (CTRL.Lq * CTRL.cmd_idq[0] + CTRL.KA) )
    # uq_inv = 1.0 / CTRL.cmd_udq[1]
    uq_inv = 1.0
    if CTRL.TUNING_IGNORE_UQ :
        uq_inv = 1.0
    CTRL.nsoaf_KD = (3 * CTRL.nsoaf_omega_ob - CTRL.R * CTRL.Lq_inv) * one_over__npp_divided_by_Js__times__Lq_id_plus_KActive                       * uq_inv
    CTRL.nsoaf_KP = ( (3 * CTRL.nsoaf_omega_ob*CTRL.nsoaf_omega_ob)         * one_over__npp_divided_by_Js__times__Lq_id_plus_KActive - 1.5 * CTRL.npp * CTRL.KA ) * uq_inv
    CTRL.nsoaf_KI = CTRL.nsoaf_omega_ob * CTRL.nsoaf_omega_ob * CTRL.nsoaf_omega_ob      * one_over__npp_divided_by_Js__times__Lq_id_plus_KActive                    * uq_inv

# 低通滤波器
def _lpf(x, y, tau_inv,CTRL):
    return y + tau_inv * (x - y) * CTRL.CL_TS
# nature speed oberver dynamics
def NSO_Dynamics(x, CTRL, param = 1.0):
    fx = np.zeros(NS_GLOBAL)
    xIq  =  x[0]
    xOmg = x[1]
    xTL  =  x[2]
    CTRL.nsoaf_uQ = CTRL.udq[1]
    CTRL.nsoaf_output_error = CTRL.idq_c[1] - xIq
    if CTRL.Ld-CTRL.Lq != 0:
        CTRL.uQ_now_filtered = _lpf(CTRL.nsoaf_uQ, CTRL.uQ_now_filtered, 30, CTRL)

        # CTRL.active_power_real =  abs(CTRL.nsoaf_uQ ) * CTRL.idq_c[1]
        # CTRL.active_power_est = abs(CTRL.nsoaf_uQ ) *xIq
        # CTRL.active_power_error =  abs(CTRL.nsoaf_uQ ) * CTRL.nsoaf_output_error
        CTRL.active_power_real =  1 * CTRL.idq_c[1]
        CTRL.active_power_est = 1 *xIq
        CTRL.active_power_error =  1 * CTRL.nsoaf_output_error
    else:
        CTRL.active_power_real =  CTRL.idq_c[1]
        CTRL.active_power_est =  xIq
        CTRL.active_power_error =  CTRL.nsoaf_output_error

    fx[0] = CTRL.Lq_inv * (CTRL.udq[1] - CTRL.R * xIq - xOmg * (CTRL.KE + CTRL.Ld * CTRL.idq_c[0]) ) - CTRL.nsoaf_KD * CTRL.active_power_error

    CTRL.nsoaf_xTem = CTRL.CLARKE_TRANS_TORQUE_GAIN * CTRL.npp * CTRL.KA * xIq
    fx[1] = CTRL.Js_inv * CTRL.npp * (CTRL.nsoaf_xTem - xTL -CTRL.nsoaf_KP * CTRL.active_power_error)
    fx[2] = CTRL.nsoaf_KI * CTRL.active_power_error
    return fx

def RK4_ObserverSolver_CJH_Style(THE_DYNAMICS, x, hs, CTRL, param=1.0):

    k1, k2, k3, k4 = np.zeros(NS_GLOBAL), np.zeros(NS_GLOBAL), np.zeros(NS_GLOBAL), np.zeros(NS_GLOBAL) # incrementals at 4 stages
    xk, fx = np.zeros(NS_GLOBAL), np.zeros(NS_GLOBAL) # state x for stage 2/3/4, state derivative

    CTRL.uab[0] = CTRL.cmd_uab[0]
    CTRL.uab[1] = CTRL.cmd_uab[1]
    CTRL.iab[0] = CTRL.iab_prev[0]
    CTRL.iab[1] = CTRL.iab_prev[1]
    fx = THE_DYNAMICS(x, CTRL, param)
    for i in range(0, NS_GLOBAL):
        k1[i] = fx[i] * hs
        xk[i] = x[i] + k1[i]*0.5

    CTRL.iab[0] = 0.5*(CTRL.iab_prev[0]+CTRL.iab_curr[0])
    CTRL.iab[1] = 0.5*(CTRL.iab_prev[1]+CTRL.iab_curr[1])
    fx = THE_DYNAMICS(xk, CTRL, param)
    for i in range(0, NS_GLOBAL):
        k2[i] = fx[i] * hs
        xk[i] = x[i] + k2[i]*0.5

    fx = THE_DYNAMICS(xk, CTRL, param)
    for i in range(0, NS_GLOBAL):
        k3[i] = fx[i] * hs
        xk[i] = x[i] + k3[i]

    CTRL.iab[0] = CTRL.iab_curr[0]
    CTRL.iab[1] = CTRL.iab_curr[1]
    fx = THE_DYNAMICS(xk, CTRL, param)
    for i in range(0, NS_GLOBAL):
        k4[i] = fx[i] * hs
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i]) * CTRL.one_over_six

def rhf_CmdErrFdkCorInFrameRho_Dynamics(x, CTRL,FE_param=1):
    
    fx = np.zeros(NS_GLOBAL)
    
    CTRL.rotor_flux_error[0] = ( CTRL.cmd_psi_mu[0] - (x[0]-CTRL.Lq * CTRL.iab[0]) )
    CTRL.rotor_flux_error[1] = ( CTRL.cmd_psi_mu[1] - (x[1]-CTRL.Lq * CTRL.iab[1]) )

    CTRL.emf_stator[0] = CTRL.uab[0] - CTRL.R * FE_param * CTRL.iab[0] + CTRL.OFFSET_VOLTAGE_ALPHA + CTRL.VM_PROPOSED_PI_CORRECTION_GAIN_P * CTRL.rotor_flux_error[0] + x[2]
    CTRL.emf_stator[1] = CTRL.uab[1] - CTRL.R * FE_param * CTRL.iab[1] + CTRL.OFFSET_VOLTAGE_BETA  + CTRL.VM_PROPOSED_PI_CORRECTION_GAIN_P * CTRL.rotor_flux_error[1] + x[3]
    fx[0] = CTRL.emf_stator[0]
    fx[1] = CTRL.emf_stator[1]
    fx[2] = CTRL.VM_PROPOSED_PI_CORRECTION_GAIN_I * CTRL.rotor_flux_error[0]
    fx[3] = CTRL.VM_PROPOSED_PI_CORRECTION_GAIN_I * CTRL.rotor_flux_error[1]
    return fx

def angle_diff(a,b):
    # ''' a and b must be within [0, 2*np.pi]'''
    _, a = divmod(a, 2*np.pi)
    _, b = divmod(b, 2*np.pi)
    d1 = a-b
    if d1 > 0:
        d2 = a - (b + 2*np.pi) # d2 is negative
    else:
        d2 = (2*np.pi + a) - b # d2 is positive
    if np.abs(d1) < np.abs(d2):
        return d1
    else:
        return d2
# 计算角度最大最小值的峰峰值和平均值
def cal_theta_error(CTRL):
    if CTRL.bool_counter_theta_error == True:
        if CTRL.counter_theta_error < 3000: 
            if CTRL.thetaerror_max  < CTRL.thetaerror:
                CTRL.thetaerror_max =  CTRL.thetaerror   
                CTRL.thetaerror_max_fin = CTRL.thetaerror_max 
            if CTRL.thetaerror_min > CTRL.thetaerror:
                CTRL.thetaerror_min = CTRL.thetaerror
                CTRL.thetaerror_min_fin = CTRL.thetaerror_min
            CTRL.thetaerror_sum += CTRL.thetaerror
            CTRL.counter_theta_error += 1
        if CTRL.counter_theta_error == 3000:
            CTRL.thetaerror_avg = CTRL.thetaerror_sum / 3000
            CTRL.counter_theta_error = 0
            CTRL.thetaerror_sum = 0
            CTRL.thetaerror_min = 0
            CTRL.thetaerror_max = 0
            CTRL.bool_counter_theta_error = False
# 计算磁链最大最小值的峰峰值和平均值
def cal_psi_error(CTRL, fe_htz, ACM):
    fe_htz.psi_e = ACM.KA * np.cos(ACM.theta_d) - fe_htz.psi_A[0]
    if CTRL.bool_counter == True:
        if CTRL.counter_psi < 3000: 
            if CTRL.psi_max  < fe_htz.psi_e:
                CTRL.psi_max =  fe_htz.psi_e    
                CTRL.psi_max_fin = CTRL.psi_max 
            if CTRL.psi_min > fe_htz.psi_e:
                CTRL.psi_min = fe_htz.psi_e
                CTRL.psi_min_fin = CTRL.psi_min
            CTRL.psi_sum += fe_htz.psi_e
            CTRL.counter_psi += 1
        if CTRL.counter_psi == 3000:
            CTRL.psi_avg = CTRL.psi_sum / 3000
            CTRL.counter_psi = 0
            CTRL.psi_sum = 0
            CTRL.psi_min = 0
            CTRL.psi_max = 0
            CTRL.bool_counter = False
'''flux estimators'''
# According to CTRL.index_voltage_model_flux_estimation, chose your flux estimator 
#0. encoder 
#1. mixed holtz and boldea to estimate K_E
#2. B. Proposed Flux Command Error feedback PI Correction in Controller Frame
#3. chen saturation time based
#4. cancel saturation and estimate K_E and u_offset from \psi_max and \psi_min

def mixed_holtz_and_boldea_to_estimate_K_E(fe_htz, CTRL, ACM, FE_param):
    RK4_ObserverSolver_CJH_Style(DYNAMICS_FluxEstimator, fe_htz.xFlux, CTRL.CL_TS, CTRL, FE_param)
    fe_htz.psi_1[0] = fe_htz.xFlux[0]
    fe_htz.psi_1[1] = fe_htz.xFlux[1]
    fe_htz.psi_s[0] = fe_htz.xFlux[4]
    fe_htz.psi_s[1] = fe_htz.xFlux[5]

    # 没有限幅的
    fe_htz.psi_A[0] = fe_htz.psi_s[0] - CTRL.Lq*CTRL.iab[0]
    fe_htz.psi_A[1] = fe_htz.psi_s[1] - CTRL.Lq*CTRL.iab[1]
    
    fe_htz.psi_A_amplitude = np.sqrt(fe_htz.psi_A[0]**2 + fe_htz.psi_A[1]**2)
    # if fe_htz.psi_A[0] < fe_htz.psi_A_amplitude and fe_htz.psi_A[0] > -1 * fe_htz.psi_A_amplitude:
    #     fe_htz.ell_compensation_flag_alpha = 1
    # if fe_htz.psi_A[1] < fe_htz.psi_A_amplitude and fe_htz.psi_A[1] > -1 * fe_htz.psi_A_amplitude:
    #     fe_htz.ell_compensation_flag_beta = 1
    # if fe_htz.ell_compensation_flag_alpha == 1 and fe_htz.ell_compensation_flag_beta == 1:
    #     fe_htz.ell_compensation_flag = 1
    # 疯狂限幅的
    fe_htz.psi_2[0] = fe_htz.psi_1[0] - CTRL.Lq*CTRL.iab[0]
    fe_htz.psi_2[1] = fe_htz.psi_1[1] - CTRL.Lq*CTRL.iab[1]
    # fe_htz.psi_2_ampl = np.sqrt(fe_htz.psi_2[0]*fe_htz.psi_2[0]+fe_htz.psi_2[1]*fe_htz.psi_2[1])

    # 限幅前求角度还是应该限幅后？
    # fe_htz.theta_d = np.arctan2(fe_htz.psi_2[1], fe_htz.psi_2[0])
    # fe_htz.cosT = np.cos(fe_htz.theta_d)
    # fe_htz.sinT = np.sin(fe_htz.theta_d)

    # fe_htz.psi_1_nonSat[0] += CTRL.CL_TS*(fe_htz.emf_stator[0])
    # fe_htz.psi_1_nonSat[1] += CTRL.CL_TS*(fe_htz.emf_stator[1])
    # fe_htz.psi_2_nonSat[0] = fe_htz.psi_1_nonSat[0] - CTRL.Lq*CTRL.iab[0]
    # fe_htz.psi_2_nonSat[1] = fe_htz.psi_1_nonSat[1] - CTRL.Lq*CTRL.iab[1]

    fe_htz.psi_aster_max = CTRL.ell

    # 限幅是针对转子磁链限幅的
    for ind in range(0,2):
        if CTRL.cmd_rpm != 0.0:
            if fe_htz.psi_2[ind]    > fe_htz.psi_aster_max: # TODO BUG呀！这里怎么可以是>应该是大于等于啊！
                fe_htz.psi_2[ind]   = fe_htz.psi_aster_max
                fe_htz.sat_max_time[ind] += CTRL.CL_TS
            elif fe_htz.psi_2[ind] < -fe_htz.psi_aster_max:
                fe_htz.psi_2[ind]   = -fe_htz.psi_aster_max
                fe_htz.sat_min_time[ind] += CTRL.CL_TS
            else:
                # 这样可以及时清零饱和时间
                if fe_htz.sat_max_time[ind]>0: fe_htz.sat_max_time[ind] -= CTRL.CL_TS
                if fe_htz.sat_min_time[ind]>0: fe_htz.sat_min_time[ind] -= CTRL.CL_TS

        # 上限饱和减去下限饱和作为误差，主要为了消除实际磁链幅值大于给定的情况，实际上这种现象在常见工况下出现次数不多。
        fe_htz.u_off_saturation_time_correction[ind] = fe_htz.sat_max_time[ind] - fe_htz.sat_min_time[ind]
        # u_offset波动会导致sat_min_time和sat_max_time的波动，这个时候最有效的办法是减少gain_off。
        # 但是同时，观察饱和时间sat_min_time等的波形，可以发现它里面也会出现一个正弦波包络线。
        if fe_htz.sat_min_time[ind] > fe_htz.maximum_of_sat_min_time[ind]: fe_htz.maximum_of_sat_min_time[ind] = fe_htz.sat_min_time[ind]
        if fe_htz.sat_max_time[ind] > fe_htz.maximum_of_sat_max_time[ind]: fe_htz.maximum_of_sat_max_time[ind] = fe_htz.sat_max_time[ind]

    # 数数，算磁链周期
    if fe_htz.psi_2[0]    > 0.0:
        fe_htz.count_positive_in_one_cycle[0] += 1
        if fe_htz.count_negative_in_one_cycle[0]!=0: 
            fe_htz.negative_cycle_in_count[0] = fe_htz.count_negative_in_one_cycle[0]; fe_htz.count_negative_in_one_cycle[0] = 0
    elif fe_htz.psi_2[0] < -0.0:
        fe_htz.count_negative_in_one_cycle[0] += 1
        if fe_htz.count_positive_in_one_cycle[0]!=0: 
            fe_htz.positive_cycle_in_count[0] = fe_htz.count_positive_in_one_cycle[0]; fe_htz.count_positive_in_one_cycle[0] = 0

    if fe_htz.psi_2[1]    > 0.0:
        fe_htz.count_positive_in_one_cycle[1] += 1
        if fe_htz.count_negative_in_one_cycle[1]!=0: 
            fe_htz.negative_cycle_in_count[1] = fe_htz.count_negative_in_one_cycle[1]; fe_htz.count_negative_in_one_cycle[1] = 0
    elif fe_htz.psi_2[1] < -0.0:
        fe_htz.count_negative_in_one_cycle[1] += 1
        if fe_htz.count_positive_in_one_cycle[1]!=0: 
            fe_htz.positive_cycle_in_count[1] = fe_htz.count_positive_in_one_cycle[1]; fe_htz.count_positive_in_one_cycle[1] = 0

    # 限幅后的转子磁链，再求取限幅后的定子磁链
    fe_htz.psi_1[0] = fe_htz.psi_2[0] + CTRL.Lq*CTRL.iab[0]
    fe_htz.psi_1[1] = fe_htz.psi_2[1] + CTRL.Lq*CTRL.iab[1]
    fe_htz.xFlux[0] = fe_htz.psi_1[0]
    fe_htz.xFlux[1] = fe_htz.psi_1[1]

    # Speed Estimation
    # if True:
    #     temp = (fe_htz.psi_1[0]*fe_htz.psi_1[0]+fe_htz.psi_1[1]*fe_htz.psi_1[1])
    #     if(temp>0.001):
    #         fe_htz.field_speed_est = - (fe_htz.psi_1[0]*-fe_htz.emf_stator[1] + fe_htz.psi_1[1]*fe_htz.emf_stator[0]) / temp
    #     temp = (fe_htz.psi_2[0]*fe_htz.psi_2[0]+fe_htz.psi_2[1]*fe_htz.psi_2[1])
    #     if(temp>0.001):
    #         fe_htz.slip_est = CTRL.motor->Rreq*(CTRL.iab[0]*-fe_htz.psi_2[1]+CTRL.iab[1]*fe_htz.psi_2[0]) / temp
    #     fe_htz.omg_est = fe_htz.field_speed_est - fe_htz.slip_est

    # TODO My proposed saturation time based correction method NOTE VERY COOL

    # Loop for alpha & beta components # destroy integer outside this loop to avoid accidentally usage 
    for ind in range(0,2):

        # /* 必须先检查是否进入levelA */
        if fe_htz.flag_pos2negLevelA[ind] == True: 
            if fe_htz.psi_2_prev[ind]<0 and fe_htz.psi_2[ind]<0: # 二次检查，磁链已经是负的了  <- 可以改为施密特触发器
                if fe_htz.flag_pos2negLevelB[ind] == False:
                    fe_htz.count_negative_cycle+=1 # fe_htz.count_positive_cycle = 0

                    # 第一次进入寻找最小值的levelB，说明最大值已经检测到。
                    fe_htz.psi_1_max[ind] = fe_htz.psi_2_max[ind] # 不区别定转子磁链，区别：psi_2是连续更新的，而psi_1是离散更新的。
                    fe_htz.Delta_t_last = fe_htz.Delta_t
                    fe_htz.Delta_t = fe_htz.time_pos2neg[ind] - fe_htz.time_pos2neg_prev[ind]
                    fe_htz.time_pos2neg_prev[ind] = fe_htz.time_pos2neg[ind] # 备份作为下次耗时参考点
                    # 初始化
                    fe_htz.flag_neg2posLevelA[ind] = False
                    fe_htz.flag_neg2posLevelB[ind] = False

                    # 注意这里是正半周到负半周切换的时候才执行一次的哦！
                    # CALCULATE_OFFSET_VOLTAGE_COMPENSATION_TERMS
                    if True:
                        fe_htz.u_off_original_lpf_input[ind]         = 0.5*(fe_htz.psi_2_min[ind] + fe_htz.psi_2_max[ind]) /  (fe_htz.Delta_t+fe_htz.Delta_t_last) 
                        fe_htz.u_off_calculated_increment[ind]       = 0.5*(fe_htz.psi_2_min[ind] + fe_htz.psi_2_max[ind]) / ((fe_htz.Delta_t+fe_htz.Delta_t_last) - (fe_htz.sat_max_time[ind]+fe_htz.sat_min_time[ind])) 
                        fe_htz.u_off_saturation_time_correction[ind] = fe_htz.sat_max_time[ind] - fe_htz.sat_min_time[ind] 
                        fe_htz.u_off_direct_calculated[ind] += (fe_htz.count_negative_cycle+fe_htz.count_positive_cycle>4) * fe_htz.u_off_calculated_increment[ind] # if(BOOL_USE_METHOD_DIFFERENCE_INPUT) 
                        # 引入 count：刚起动时的几个磁链正负半周里，Delta_t_last 存在巨大的计算误差，所以要放弃更新哦。

                    # fe_htz.accumulated__u_off_saturation_time_correction[ind] += fe_htz.u_off_saturation_time_correction[ind]
                    fe_htz.sign__u_off_saturation_time_correction[ind] = -1.0
                    # 饱和时间的正弦包络线的正负半周的频率比磁链频率低多啦！需要再额外加一个低频u_offset校正
                    # 以下自适应律是在速度处于稳态情况下运行
                    fe_htz.sat_time_offset[ind] = fe_htz.maximum_of_sat_max_time[ind] - fe_htz.maximum_of_sat_min_time[ind]
                
                    if fe_htz.maximum_of_sat_max_time[ind] != 0 and fe_htz.maximum_of_sat_min_time[ind] != 0:
                        CTRL.ell = CTRL.ell + CTRL.K1_for_max_ell * CTRL.CL_TS * fe_htz.maximum_of_sat_max_time[ind] + 1000 * CTRL.CL_TS * fe_htz.maximum_of_sat_min_time[ind]
                    elif (fe_htz.maximum_of_sat_max_time[ind] == 0 and fe_htz.maximum_of_sat_min_time[ind] == 0) or fe_htz.ell_compensation_flag == 1:
                        CTRL.ell = CTRL.ell + CTRL.K2_for_min_ell * CTRL.CL_TS * (fe_htz.negative_cycle_in_count[0] + fe_htz.positive_cycle_in_count[0] + fe_htz.negative_cycle_in_count[1] + fe_htz.positive_cycle_in_count[1])
                        # CTRL.ell = CTRL.ell + CTRL.K2_for_min_ell * CTRL.CL_TS * fe_htz.positive_cycle_in_count[0]
                    # if ind==0
                        #     print(f'{CTRL.timebase}')
                        # if CTRL.timebase < 1:
                    #     if fe_htz.maximum_of_sat_max_time[ind] != 0 and fe_htz.maximum_of_sat_min_time[ind] != 0:
                    #         CTRL.ell = CTRL.ell + 1000 * CTRL.CL_TS * fe_htz.maximum_of_sat_max_time[ind] + 1000 * CTRL.CL_TS * fe_htz.maximum_of_sat_min_time[ind]
                    #     elif fe_htz.maximum_of_sat_max_time[ind] == 0 and fe_htz.maximum_of_sat_min_time[ind] == 0:
                    #         CTRL.ell = CTRL.ell - 10 * CTRL.CL_TS - 10 * CTRL.CL_TS     
                    # else:    
                    # if fe_htz.maximum_of_sat_max_time[ind] != 0 and fe_htz.maximum_of_sat_min_time[ind] != 0:
                    #     CTRL.ell = CTRL.ell + 0.003 * fe_htz.psi_A_amplitude 
                    # elif fe_htz.maximum_of_sat_max_time[ind] == 0 and fe_htz.maximum_of_sat_min_time[ind] == 0:
                    #     CTRL.ell = CTRL.ell_prev - 0.003 * np.abs(((CTRL.ell + fe_htz.psi_2_min[0]) + (CTRL.ell - fe_htz.psi_2_max[0])) * 0.5 )
                    #     CTRL.ell_prev = CTRL.ell

                    fe_htz.maximum_of_sat_max_time[ind] = 0.0
                    fe_htz.maximum_of_sat_min_time[ind] = 0.0

                    fe_htz.psi_1_min[ind] = 0.0
                    fe_htz.psi_2_min[ind] = 0.0
                    if False: # BOOL_TURN_ON_ADAPTIVE_EXTRA_LIMIT):
                        fe_htz.sat_min_time_reg[ind] = fe_htz.sat_min_time[ind]
                        if(fe_htz.sat_max_time_reg[ind]>CL_TS and fe_htz.sat_min_time_reg[ind]>CL_TS):
                            fe_htz.flag_limit_too_low = True
                            fe_htz.extra_limit += 1e-2 * (fe_htz.sat_max_time_reg[ind] + fe_htz.sat_min_time_reg[ind]) / fe_htz.Delta_t 
                        else:
                            fe_htz.flag_limit_too_low = False
                            fe_htz.extra_limit -= 2e-4 * fe_htz.Delta_t
                            if(bool_positive_extra_limit):
                                if(fe_htz.extra_limit<0.0):
                                    fe_htz.extra_limit = 0.0

                        fe_htz.sat_max_time_reg[ind] = 0.0

                    fe_htz.sat_min_time[ind] = 0.0

                fe_htz.flag_pos2negLevelB[ind] = True
                if fe_htz.flag_pos2negLevelB[ind] == True: # 寻找磁链最小值
                    if fe_htz.psi_2[ind] < fe_htz.psi_2_min[ind]:
                        fe_htz.psi_2_min[ind] = fe_htz.psi_2[ind]

            else: # 磁链还没有变负，说明是虚假过零，比如在震荡，fe_htz.psi_2[0]>0
                fe_htz.flag_pos2negLevelA[ind] = False # /* 震荡的话，另一方的检测就有可能被触动？ */

        if fe_htz.psi_2_prev[ind]>0 and fe_htz.psi_2[ind]<0: # 发现磁链由正变负的时刻
            fe_htz.flag_pos2negLevelA[ind] = True
            fe_htz.time_pos2neg[ind] = CTRL.timebase

        if fe_htz.flag_neg2posLevelA[ind] == True:
            if fe_htz.psi_2_prev[ind]>0 and fe_htz.psi_2[ind]>0: # 二次检查，磁链已经是正的了
                if fe_htz.flag_neg2posLevelB[ind] == False:
                    fe_htz.count_positive_cycle+=1 # fe_htz.count_negative_cycle = 0
                    # 第一次进入寻找最大值的levelB，说明最小值已经检测到。
                    fe_htz.psi_1_min[ind] = fe_htz.psi_2_min[ind] # 不区别定转子磁链，区别：psi_2是连续更新的，而psi_1是离散更新的。
                    fe_htz.Delta_t_last = fe_htz.Delta_t
                    fe_htz.Delta_t = fe_htz.time_neg2pos[ind] - fe_htz.time_neg2pos_prev[ind]
                    fe_htz.time_neg2pos_prev[ind] = fe_htz.time_neg2pos[ind] # 备份作为下次耗时参考点
                    # 初始化
                    fe_htz.flag_pos2negLevelA[ind] = False
                    fe_htz.flag_pos2negLevelB[ind] = False

                    if True:
                        fe_htz.u_off_original_lpf_input[ind]         = 0.5*(fe_htz.psi_2_min[ind] + fe_htz.psi_2_max[ind]) /  (fe_htz.Delta_t+fe_htz.Delta_t_last) 
                        fe_htz.u_off_calculated_increment[ind]       = 0.5*(fe_htz.psi_2_min[ind] + fe_htz.psi_2_max[ind]) / ((fe_htz.Delta_t+fe_htz.Delta_t_last) - (fe_htz.sat_max_time[ind]+fe_htz.sat_min_time[ind])) 
                        fe_htz.u_off_saturation_time_correction[ind] = fe_htz.sat_max_time[ind] - fe_htz.sat_min_time[ind] 
                        fe_htz.u_off_direct_calculated[ind] += (fe_htz.count_negative_cycle+fe_htz.count_positive_cycle>4) * fe_htz.u_off_calculated_increment[ind] # if(BOOL_USE_METHOD_DIFFERENCE_INPUT) 
                        # 引入 count：刚起动时的几个磁链正负半周里，Delta_t_last 存在巨大的计算误差，所以要放弃更新哦。

                    # fe_htz.accumulated__u_off_saturation_time_correction[ind] += fe_htz.u_off_saturation_time_correction[ind]
                    fe_htz.sign__u_off_saturation_time_correction[ind] = 1.0

                    fe_htz.psi_1_max[ind] = 0.0
                    fe_htz.psi_2_max[ind] = 0.0
                    if False: # BOOL_TURN_ON_ADAPTIVE_EXTRA_LIMIT):
                        fe_htz.sat_max_time_reg[ind] = fe_htz.sat_max_time[ind]
                        if(fe_htz.sat_min_time_reg[ind]>CL_TS and fe_htz.sat_max_time_reg[ind]>CL_TS):
                            fe_htz.flag_limit_too_low = True
                            fe_htz.extra_limit += 1e-2 * (fe_htz.sat_min_time_reg[ind] + fe_htz.sat_max_time_reg[ind]) / fe_htz.Delta_t 
                        else:
                            fe_htz.flag_limit_too_low = False
                            fe_htz.extra_limit -= 2e-4 * fe_htz.Delta_t
                            if(fe_htz.extra_limit<0.0):
                                fe_htz.extra_limit = 0.0

                        fe_htz.sat_min_time_reg[ind] = 0.0

                    fe_htz.sat_max_time[ind] = 0.0

                fe_htz.flag_neg2posLevelB[ind] = True 
                if fe_htz.flag_neg2posLevelB[ind] == True: # 寻找磁链最大值
                    if fe_htz.psi_2[ind] > fe_htz.psi_2_max[ind]:
                        fe_htz.psi_2_max[ind] = fe_htz.psi_2[ind]

            else: # 磁链还没有变正，说明是虚假过零，比如在震荡，fe_htz.psi_2[0]<0
                fe_htz.flag_neg2posLevelA[ind] = False

        if fe_htz.psi_2_prev[ind]<0 and fe_htz.psi_2[ind]>0: # 发现磁链由负变正的时刻
            fe_htz.flag_neg2posLevelA[ind] = True
            fe_htz.time_neg2pos[ind] = CTRL.timebase

    # /*这里一共有四种方案，积分两种，LPF两种：
    # 1. Holtz03原版是用u_off_original_lpf_input过LPF，
    # 2. 我发现u_off_original_lpf_input过积分器才能完全补偿偏置电压，
    # 3. 我还提出可以直接算出偏置电压补偿误差（可加LPF），
    # 4. 我还提出了用饱和时间去做校正的方法*/

    INTEGRAL_INPUT_ALPHA = fe_htz.u_off_saturation_time_correction[0] # exact offset calculation for compensation
    INTEGRAL_INPUT_BETA  = fe_htz.u_off_saturation_time_correction[1] # exact offset calculation for compensation

    if fe_htz.GAIN_OFFSET_REALTIME != 0.0:
        integer_local_sum = fe_htz.negative_cycle_in_count[0] + fe_htz.positive_cycle_in_count[0] + fe_htz.negative_cycle_in_count[1] + fe_htz.positive_cycle_in_count[1]
        if integer_local_sum>0:
            fe_htz.gain_off = fe_htz.GAIN_OFFSET_REALTIME * fe_htz.GAIN_OFFSET_INIT / (integer_local_sum*CTRL.CL_TS)
    else:
        fe_htz.gain_off = fe_htz.GAIN_OFFSET_INIT

    fe_htz.u_offset[0] += 0.9*fe_htz.gain_off * 100 * CTRL.CL_TS * INTEGRAL_INPUT_ALPHA
    fe_htz.u_offset[1] += 0.9*fe_htz.gain_off * 100 * CTRL.CL_TS * INTEGRAL_INPUT_BETA
    fe_htz.xFlux[2] = fe_htz.u_offset[0]
    fe_htz.xFlux[3] = fe_htz.u_offset[1]

    fe_htz.psi_2_prev[0] = fe_htz.psi_2[0]
    fe_htz.psi_2_prev[1] = fe_htz.psi_2[1]

    # psi_2_ampl 在限幅前已经算过了，还有必要限幅后在这里再算一次吗？
    # fe_htz.psi_A_ampl = np.sqrt(fe_htz.psi_A[0] **2 + fe_htz.psi_A[1] **2)
    # if fe_htz.psi_A_ampl == 0:
    #     fe_htz.psi_A_ampl = 1.0
    # amplitude_inverse = 1.0 / fe_htz.psi_A_ampl
    fe_htz.psi_A_ampl = np.sqrt(fe_htz.psi_A[0] **2 + fe_htz.psi_A[1] **2)
    CTRL.flux_estimate_amplitude = fe_htz.psi_A_ampl
    if fe_htz.psi_A_ampl == 0:
        fe_htz.psi_A_ampl = 1.0
    amplitude_inverse = 1.0 / fe_htz.psi_A_ampl

    CTRL.cosT = fe_htz.psi_A[0] * amplitude_inverse
    CTRL.sinT = fe_htz.psi_A[1] * amplitude_inverse
    fe_htz.theta_d = np.arctan2(fe_htz.psi_A[1], fe_htz.psi_A[0]) # Costly operation, but it is needed only once per control interrupt
    CTRL.theta_d = fe_htz.theta_d
    CTRL.cosT = np.cos(CTRL.theta_d)
    CTRL.sinT = np.sin(CTRL.theta_d)

    while ACM.theta_d> np.pi: ACM.theta_d -= 2*np.pi
    while ACM.theta_d<-np.pi: ACM.theta_d += 2*np.pi

    if CTRL.theta_d * ACM.theta_d > 0:
        CTRL.thetaerror = np.sin(ACM.theta_d - CTRL.theta_d)
    elif CTRL.theta_d * ACM.theta_d < 0:
        CTRL.thetaerror = np.sin(CTRL.theta_d + ACM.theta_d)
    
    cal_theta_error(CTRL)
    cal_psi_error(CTRL, fe_htz, ACM)
    if CTRL.use_encoder_angle_no_matter_what == True:
        CTRL.theta_d = ACM.theta_d
        CTRL.cosT = np.cos(CTRL.theta_d)
        CTRL.sinT = np.sin(CTRL.theta_d)

def B_Proposed_Flux_Command_Error_feedback_PI_Correction(fe_htz, CTRL, ACM, FE_param):

    CTRL.cmd_psi_mu[0] = CTRL.cmd_psi * CTRL.cosT
    CTRL.cmd_psi_mu[1] = CTRL.cmd_psi * CTRL.sinT 
    RK4_ObserverSolver_CJH_Style(rhf_CmdErrFdkCorInFrameRho_Dynamics, fe_htz.xFlux, CTRL.CL_TS, CTRL, FE_param)
    #// Unpack x
    fe_htz.psi_1[0]                         = fe_htz.xFlux[0]
    fe_htz.psi_1[1]                         = fe_htz.xFlux[1]
    CTRL.correction_integral_term[0]        = fe_htz.xFlux[2]
    CTRL.correction_integral_term[1]        = fe_htz.xFlux[3]
    fe_htz.u_offset[0] = CTRL.correction_integral_term[0]
    fe_htz.u_offset[1] = CTRL.correction_integral_term[1]

    #// rotor flux updates

    fe_htz.psi_2[0] = fe_htz.psi_1[0] - CTRL.Lq * CTRL.iab_curr[0]
    fe_htz.psi_2[1] = fe_htz.psi_1[1] - CTRL.Lq * CTRL.iab_curr[1]

    CTRL.theta_d = np.arctan2(fe_htz.psi_2[1], fe_htz.psi_2[0]) 

    while ACM.theta_d> np.pi: ACM.theta_d -= 2*np.pi
    while ACM.theta_d<-np.pi: ACM.theta_d += 2*np.pi
    fe_htz.theta_d = ACM.theta_d
    if CTRL.theta_d * fe_htz.theta_d > 0:
        CTRL.thetaerror = np.sin(CTRL.theta_d - fe_htz.theta_d)
    elif CTRL.theta_d * fe_htz.theta_d < 0:
        CTRL.thetaerror = np.sin(CTRL.theta_d + fe_htz.theta_d)
    
    cal_theta_error(CTRL)
    cal_psi_error(CTRL, fe_htz, ACM)

    if CTRL.use_encoder_angle_no_matter_what == True:
        CTRL.theta_d = ACM.theta_d
        CTRL.cosT = np.cos(CTRL.theta_d)
        CTRL.sinT = np.sin(CTRL.theta_d)

def chen_saturation_time_based(fe_htz, CTRL, ACM, FE_param, ELL_param):
    # sensorless
    RK4_ObserverSolver_CJH_Style(DYNAMICS_FluxEstimator, fe_htz.xFlux, CTRL.CL_TS, CTRL, FE_param)
    fe_htz.psi_1[0] = fe_htz.xFlux[0]
    fe_htz.psi_1[1] = fe_htz.xFlux[1]

    # 疯狂限幅的
    fe_htz.psi_2[0] = fe_htz.psi_1[0] - CTRL.Lq*CTRL.iab[0]
    fe_htz.psi_2[1] = fe_htz.psi_1[1] - CTRL.Lq*CTRL.iab[1]
    # fe_htz.psi_2_ampl = np.sqrt(fe_htz.psi_2[0]*fe_htz.psi_2[0]+fe_htz.psi_2[1]*fe_htz.psi_2[1])
    fe_htz.psi_A[0] = fe_htz.psi_2[0]
    fe_htz.psi_A[1] = fe_htz.psi_2[1]
    # 限幅前求角度还是应该限幅后？
    # fe_htz.theta_d = np.arctan2(fe_htz.psi_2[1], fe_htz.psi_2[0])
    # fe_htz.cosT = np.cos(fe_htz.theta_d)
    # fe_htz.sinT = np.sin(fe_htz.theta_d)

    # fe_htz.psi_1_nonSat[0] += CTRL.CL_TS*(fe_htz.emf_stator[0])
    # fe_htz.psi_1_nonSat[1] += CTRL.CL_TS*(fe_htz.emf_stator[1])
    # fe_htz.psi_2_nonSat[0] = fe_htz.psi_1_nonSat[0] - CTRL.Lq*CTRL.iab[0]
    # fe_htz.psi_2_nonSat[1] = fe_htz.psi_1_nonSat[1] - CTRL.Lq*CTRL.iab[1]

    fe_htz.psi_aster_max = ELL_param

    # 限幅是针对转子磁链限幅的
    for ind in range(0,2):
        if CTRL.cmd_rpm != 0.0:
            if fe_htz.psi_2[ind]    > fe_htz.psi_aster_max: # TODO BUG呀！这里怎么可以是>应该是大于等于啊！
                fe_htz.psi_2[ind]   = fe_htz.psi_aster_max
                fe_htz.sat_max_time[ind] += CTRL.CL_TS
            elif fe_htz.psi_2[ind] < -fe_htz.psi_aster_max:
                fe_htz.psi_2[ind]   = -fe_htz.psi_aster_max
                fe_htz.sat_min_time[ind] += CTRL.CL_TS
            else:
                # 这样可以及时清零饱和时间
                if fe_htz.sat_max_time[ind]>0: fe_htz.sat_max_time[ind] -= CTRL.CL_TS
                if fe_htz.sat_min_time[ind]>0: fe_htz.sat_min_time[ind] -= CTRL.CL_TS

        # 上限饱和减去下限饱和作为误差，主要为了消除实际磁链幅值大于给定的情况，实际上这种现象在常见工况下出现次数不多。
        fe_htz.u_off_saturation_time_correction[ind] = fe_htz.sat_max_time[ind] - fe_htz.sat_min_time[ind]
        # u_offset波动会导致sat_min_time和sat_max_time的波动，这个时候最有效的办法是减少gain_off。
        # 但是同时，观察饱和时间sat_min_time等的波形，可以发现它里面也会出现一个正弦波包络线。
        if fe_htz.sat_min_time[ind] > fe_htz.maximum_of_sat_min_time[ind]: fe_htz.maximum_of_sat_min_time[ind] = fe_htz.sat_min_time[ind]
        if fe_htz.sat_max_time[ind] > fe_htz.maximum_of_sat_max_time[ind]: fe_htz.maximum_of_sat_max_time[ind] = fe_htz.sat_max_time[ind]

    # 数数，算磁链周期
    if fe_htz.psi_2[0]    > 0.0:
        fe_htz.count_positive_in_one_cycle[0] += 1
        if fe_htz.count_negative_in_one_cycle[0]!=0: 
            fe_htz.negative_cycle_in_count[0] = fe_htz.count_negative_in_one_cycle[0]; fe_htz.count_negative_in_one_cycle[0] = 0
    elif fe_htz.psi_2[0] < -0.0:
        fe_htz.count_negative_in_one_cycle[0] += 1
        if fe_htz.count_positive_in_one_cycle[0]!=0: 
            fe_htz.positive_cycle_in_count[0] = fe_htz.count_positive_in_one_cycle[0]; fe_htz.count_positive_in_one_cycle[0] = 0

    if fe_htz.psi_2[1]    > 0.0:
        fe_htz.count_positive_in_one_cycle[1] += 1
        if fe_htz.count_negative_in_one_cycle[1]!=0: 
            fe_htz.negative_cycle_in_count[1] = fe_htz.count_negative_in_one_cycle[1]; fe_htz.count_negative_in_one_cycle[1] = 0
    elif fe_htz.psi_2[1] < -0.0:
        fe_htz.count_negative_in_one_cycle[1] += 1
        if fe_htz.count_positive_in_one_cycle[1]!=0: 
            fe_htz.positive_cycle_in_count[1] = fe_htz.count_positive_in_one_cycle[1]; fe_htz.count_positive_in_one_cycle[1] = 0

    # 限幅后的转子磁链，再求取限幅后的定子磁链
    fe_htz.psi_1[0] = fe_htz.psi_2[0] + CTRL.Lq*CTRL.iab[0]
    fe_htz.psi_1[1] = fe_htz.psi_2[1] + CTRL.Lq*CTRL.iab[1]
    fe_htz.xFlux[0] = fe_htz.psi_1[0]
    fe_htz.xFlux[1] = fe_htz.psi_1[1]

    # Speed Estimation
    # if True:
    #     temp = (fe_htz.psi_1[0]*fe_htz.psi_1[0]+fe_htz.psi_1[1]*fe_htz.psi_1[1])
    #     if(temp>0.001):
    #         fe_htz.field_speed_est = - (fe_htz.psi_1[0]*-fe_htz.emf_stator[1] + fe_htz.psi_1[1]*fe_htz.emf_stator[0]) / temp
    #     temp = (fe_htz.psi_2[0]*fe_htz.psi_2[0]+fe_htz.psi_2[1]*fe_htz.psi_2[1])
    #     if(temp>0.001):
    #         fe_htz.slip_est = CTRL.motor->Rreq*(CTRL.iab[0]*-fe_htz.psi_2[1]+CTRL.iab[1]*fe_htz.psi_2[0]) / temp
    #     fe_htz.omg_est = fe_htz.field_speed_est - fe_htz.slip_est

    # TODO My proposed saturation time based correction method NOTE VERY COOL

    # Loop for alpha & beta components # destroy integer outside this loop to avoid accidentally usage 
    for ind in range(0,2):

        # /* 必须先检查是否进入levelA */
        if fe_htz.flag_pos2negLevelA[ind] == True: 
            if fe_htz.psi_2_prev[ind]<0 and fe_htz.psi_2[ind]<0: # 二次检查，磁链已经是负的了  <- 可以改为施密特触发器
                if fe_htz.flag_pos2negLevelB[ind] == False:
                    fe_htz.count_negative_cycle+=1 # fe_htz.count_positive_cycle = 0

                    # 第一次进入寻找最小值的levelB，说明最大值已经检测到。
                    fe_htz.psi_1_max[ind] = fe_htz.psi_2_max[ind] # 不区别定转子磁链，区别：psi_2是连续更新的，而psi_1是离散更新的。
                    fe_htz.Delta_t_last = fe_htz.Delta_t
                    fe_htz.Delta_t = fe_htz.time_pos2neg[ind] - fe_htz.time_pos2neg_prev[ind]
                    fe_htz.time_pos2neg_prev[ind] = fe_htz.time_pos2neg[ind] # 备份作为下次耗时参考点
                    # 初始化
                    fe_htz.flag_neg2posLevelA[ind] = False
                    fe_htz.flag_neg2posLevelB[ind] = False

                    # 注意这里是正半周到负半周切换的时候才执行一次的哦！
                    # CALCULATE_OFFSET_VOLTAGE_COMPENSATION_TERMS
                    if True:
                        fe_htz.u_off_original_lpf_input[ind]         = 0.5*(fe_htz.psi_2_min[ind] + fe_htz.psi_2_max[ind]) /  (fe_htz.Delta_t+fe_htz.Delta_t_last) 
                        fe_htz.u_off_calculated_increment[ind]       = 0.5*(fe_htz.psi_2_min[ind] + fe_htz.psi_2_max[ind]) / ((fe_htz.Delta_t+fe_htz.Delta_t_last) - (fe_htz.sat_max_time[ind]+fe_htz.sat_min_time[ind])) 
                        fe_htz.u_off_saturation_time_correction[ind] = fe_htz.sat_max_time[ind] - fe_htz.sat_min_time[ind] 
                        fe_htz.u_off_direct_calculated[ind] += (fe_htz.count_negative_cycle+fe_htz.count_positive_cycle>4) * fe_htz.u_off_calculated_increment[ind] # if(BOOL_USE_METHOD_DIFFERENCE_INPUT) 
                        # 引入 count：刚起动时的几个磁链正负半周里，Delta_t_last 存在巨大的计算误差，所以要放弃更新哦。

                    # fe_htz.accumulated__u_off_saturation_time_correction[ind] += fe_htz.u_off_saturation_time_correction[ind]
                    fe_htz.sign__u_off_saturation_time_correction[ind] = -1.0
                    # 饱和时间的正弦包络线的正负半周的频率比磁链频率低多啦！需要再额外加一个低频u_offset校正
                    fe_htz.sat_time_offset[ind] = fe_htz.maximum_of_sat_max_time[ind] - fe_htz.maximum_of_sat_min_time[ind]
                    fe_htz.maximum_of_sat_max_time[ind] = 0.0
                    fe_htz.maximum_of_sat_min_time[ind] = 0.0

                    fe_htz.psi_1_min[ind] = 0.0
                    fe_htz.psi_2_min[ind] = 0.0
                    if False: # BOOL_TURN_ON_ADAPTIVE_EXTRA_LIMIT):
                        fe_htz.sat_min_time_reg[ind] = fe_htz.sat_min_time[ind]
                        if(fe_htz.sat_max_time_reg[ind]>CL_TS and fe_htz.sat_min_time_reg[ind]>CL_TS):
                            fe_htz.flag_limit_too_low = True
                            fe_htz.extra_limit += 1e-2 * (fe_htz.sat_max_time_reg[ind] + fe_htz.sat_min_time_reg[ind]) / fe_htz.Delta_t 
                        else:
                            fe_htz.flag_limit_too_low = False
                            fe_htz.extra_limit -= 2e-4 * fe_htz.Delta_t
                            if(bool_positive_extra_limit):
                                if(fe_htz.extra_limit<0.0):
                                    fe_htz.extra_limit = 0.0

                        fe_htz.sat_max_time_reg[ind] = 0.0

                    fe_htz.sat_min_time[ind] = 0.0

                fe_htz.flag_pos2negLevelB[ind] = True
                if fe_htz.flag_pos2negLevelB[ind] == True: # 寻找磁链最小值
                    if fe_htz.psi_2[ind] < fe_htz.psi_2_min[ind]:
                        fe_htz.psi_2_min[ind] = fe_htz.psi_2[ind]

            else: # 磁链还没有变负，说明是虚假过零，比如在震荡，fe_htz.psi_2[0]>0
                fe_htz.flag_pos2negLevelA[ind] = False # /* 震荡的话，另一方的检测就有可能被触动？ */

        if fe_htz.psi_2_prev[ind]>0 and fe_htz.psi_2[ind]<0: # 发现磁链由正变负的时刻
            fe_htz.flag_pos2negLevelA[ind] = True
            fe_htz.time_pos2neg[ind] = CTRL.timebase

        if fe_htz.flag_neg2posLevelA[ind] == True:
            if fe_htz.psi_2_prev[ind]>0 and fe_htz.psi_2[ind]>0: # 二次检查，磁链已经是正的了
                if fe_htz.flag_neg2posLevelB[ind] == False:
                    fe_htz.count_positive_cycle+=1 # fe_htz.count_negative_cycle = 0
                    # 第一次进入寻找最大值的levelB，说明最小值已经检测到。
                    fe_htz.psi_1_min[ind] = fe_htz.psi_2_min[ind] # 不区别定转子磁链，区别：psi_2是连续更新的，而psi_1是离散更新的。
                    fe_htz.Delta_t_last = fe_htz.Delta_t
                    fe_htz.Delta_t = fe_htz.time_neg2pos[ind] - fe_htz.time_neg2pos_prev[ind]
                    fe_htz.time_neg2pos_prev[ind] = fe_htz.time_neg2pos[ind] # 备份作为下次耗时参考点
                    # 初始化
                    fe_htz.flag_pos2negLevelA[ind] = False
                    fe_htz.flag_pos2negLevelB[ind] = False

                    if True:
                        fe_htz.u_off_original_lpf_input[ind]         = 0.5*(fe_htz.psi_2_min[ind] + fe_htz.psi_2_max[ind]) /  (fe_htz.Delta_t+fe_htz.Delta_t_last) 
                        fe_htz.u_off_calculated_increment[ind]       = 0.5*(fe_htz.psi_2_min[ind] + fe_htz.psi_2_max[ind]) / ((fe_htz.Delta_t+fe_htz.Delta_t_last) - (fe_htz.sat_max_time[ind]+fe_htz.sat_min_time[ind])) 
                        fe_htz.u_off_saturation_time_correction[ind] = fe_htz.sat_max_time[ind] - fe_htz.sat_min_time[ind] 
                        fe_htz.u_off_direct_calculated[ind] += (fe_htz.count_negative_cycle+fe_htz.count_positive_cycle>4) * fe_htz.u_off_calculated_increment[ind] # if(BOOL_USE_METHOD_DIFFERENCE_INPUT) 
                        # 引入 count：刚起动时的几个磁链正负半周里，Delta_t_last 存在巨大的计算误差，所以要放弃更新哦。

                    # fe_htz.accumulated__u_off_saturation_time_correction[ind] += fe_htz.u_off_saturation_time_correction[ind]
                    fe_htz.sign__u_off_saturation_time_correction[ind] = 1.0

                    fe_htz.psi_1_max[ind] = 0.0
                    fe_htz.psi_2_max[ind] = 0.0
                    if False: # BOOL_TURN_ON_ADAPTIVE_EXTRA_LIMIT):
                        fe_htz.sat_max_time_reg[ind] = fe_htz.sat_max_time[ind]
                        if(fe_htz.sat_min_time_reg[ind]>CL_TS and fe_htz.sat_max_time_reg[ind]>CL_TS):
                            fe_htz.flag_limit_too_low = True
                            fe_htz.extra_limit += 1e-2 * (fe_htz.sat_min_time_reg[ind] + fe_htz.sat_max_time_reg[ind]) / fe_htz.Delta_t 
                        else:
                            fe_htz.flag_limit_too_low = False
                            fe_htz.extra_limit -= 2e-4 * fe_htz.Delta_t
                            if(fe_htz.extra_limit<0.0):
                                fe_htz.extra_limit = 0.0

                        fe_htz.sat_min_time_reg[ind] = 0.0

                    fe_htz.sat_max_time[ind] = 0.0

                fe_htz.flag_neg2posLevelB[ind] = True 
                if fe_htz.flag_neg2posLevelB[ind] == True: # 寻找磁链最大值
                    if fe_htz.psi_2[ind] > fe_htz.psi_2_max[ind]:
                        fe_htz.psi_2_max[ind] = fe_htz.psi_2[ind]

            else: # 磁链还没有变正，说明是虚假过零，比如在震荡，fe_htz.psi_2[0]<0
                fe_htz.flag_neg2posLevelA[ind] = False

        if fe_htz.psi_2_prev[ind]<0 and fe_htz.psi_2[ind]>0: # 发现磁链由负变正的时刻
            fe_htz.flag_neg2posLevelA[ind] = True
            fe_htz.time_neg2pos[ind] = CTRL.timebase

    # /*这里一共有四种方案，积分两种，LPF两种：
    # 1. Holtz03原版是用u_off_original_lpf_input过LPF，
    # 2. 我发现u_off_original_lpf_input过积分器才能完全补偿偏置电压，
    # 3. 我还提出可以直接算出偏置电压补偿误差（可加LPF），
    # 4. 我还提出了用饱和时间去做校正的方法*/

    INTEGRAL_INPUT_ALPHA = fe_htz.u_off_saturation_time_correction[0] # exact offset calculation for compensation
    INTEGRAL_INPUT_BETA  = fe_htz.u_off_saturation_time_correction[1] # exact offset calculation for compensation

    if fe_htz.GAIN_OFFSET_REALTIME != 0.0:
        integer_local_sum = fe_htz.negative_cycle_in_count[0] + fe_htz.positive_cycle_in_count[0] + fe_htz.negative_cycle_in_count[1] + fe_htz.positive_cycle_in_count[1]
        if integer_local_sum>0:
            fe_htz.gain_off = fe_htz.GAIN_OFFSET_REALTIME * fe_htz.GAIN_OFFSET_INIT / (integer_local_sum*CTRL.CL_TS)
    else:
        fe_htz.gain_off = fe_htz.GAIN_OFFSET_INIT

    fe_htz.u_offset[0] += 0.9*fe_htz.gain_off * 10 * CTRL.CL_TS * INTEGRAL_INPUT_ALPHA
    fe_htz.u_offset[1] += 0.9*fe_htz.gain_off * 10 * CTRL.CL_TS * INTEGRAL_INPUT_BETA
    fe_htz.xFlux[2] = fe_htz.u_offset[0]
    fe_htz.xFlux[3] = fe_htz.u_offset[1]

    fe_htz.psi_2_prev[0] = fe_htz.psi_2[0]
    fe_htz.psi_2_prev[1] = fe_htz.psi_2[1]

    # psi_2_ampl 在限幅前已经算过了，还有必要限幅后在这里再算一次吗？
    # fe_htz.psi_A_ampl = np.sqrt(fe_htz.psi_A[0] **2 + fe_htz.psi_A[1] **2)
    # if fe_htz.psi_A_ampl == 0:
    #     fe_htz.psi_A_ampl = 1.0
    # amplitude_inverse = 1.0 / fe_htz.psi_A_ampl
    fe_htz.psi_2_ampl = np.sqrt(fe_htz.psi_2[0] **2 + fe_htz.psi_2[1] **2)
    CTRL.flux_estimate_amplitude = fe_htz.psi_2_ampl
    if fe_htz.psi_2_ampl == 0:
        fe_htz.psi_2_ampl = 1.0
    amplitude_inverse = 1.0 / fe_htz.psi_2_ampl

    CTRL.cosT = fe_htz.psi_2[0] * amplitude_inverse
    CTRL.sinT = fe_htz.psi_2[1] * amplitude_inverse
    fe_htz.theta_d = np.arctan2(fe_htz.psi_2[1], fe_htz.psi_2[0]) # Costly operation, but it is needed only once per control interrupt
    CTRL.theta_d = fe_htz.theta_d
    CTRL.cosT = np.cos(CTRL.theta_d)
    CTRL.sinT = np.sin(CTRL.theta_d)

    while ACM.theta_d> np.pi: ACM.theta_d -= 2*np.pi
    while ACM.theta_d<-np.pi: ACM.theta_d += 2*np.pi
    
    if CTRL.theta_d * ACM.theta_d > 0:
        CTRL.thetaerror = np.sin(ACM.theta_d - CTRL.theta_d)
    elif CTRL.theta_d * ACM.theta_d < 0:
        CTRL.thetaerror = np.sin(CTRL.theta_d + ACM.theta_d)
    
    if CTRL.use_encoder_angle_no_matter_what == True:
        CTRL.theta_d = ACM.theta_d
        CTRL.cosT = np.cos(CTRL.theta_d)
        CTRL.sinT = np.sin(CTRL.theta_d)

def no_saturation_to_cal_KE_and_uoffset(fe_htz, CTRL, ACM, FE_param):
    RK4_ObserverSolver_CJH_Style(DYNAMICS_FluxEstimator, fe_htz.xFlux, CTRL.CL_TS, CTRL, FE_param)
    fe_htz.psi_1[0] = fe_htz.xFlux[0]
    fe_htz.psi_1[1] = fe_htz.xFlux[1]

    # 疯狂限幅的
    fe_htz.psi_2[0] = fe_htz.psi_1[0] - CTRL.Lq*CTRL.iab[0]
    fe_htz.psi_2[1] = fe_htz.psi_1[1] - CTRL.Lq*CTRL.iab[1]

    for ind in range(0,2):
        # /* 必须先检查是否进入levelA */
        if fe_htz.flag_pos2negLevelA[ind] == True: 
            if fe_htz.psi_2_prev[ind]<0 and fe_htz.psi_2[ind]<0: # 二次检查，磁链已经是负的了  <- 可以改为施密特触发器
                if fe_htz.flag_pos2negLevelB[ind] == False:
                    fe_htz.count_negative_cycle+=1 # fe_htz.count_positive_cycle = 0

                    # 第一次进入寻找最小值的levelB，说明最大值已经检测到。
                    fe_htz.psi_1_max[ind] = fe_htz.psi_2_max[ind] # 不区别定转子磁链，区别：psi_2是连续更新的，而psi_1是离散更新的。
                    fe_htz.Delta_t_last = fe_htz.Delta_t
                    fe_htz.Delta_t = fe_htz.time_pos2neg[ind] - fe_htz.time_pos2neg_prev[ind]
                    fe_htz.time_pos2neg_prev[ind] = fe_htz.time_pos2neg[ind] # 备份作为下次耗时参考点
                    # 初始化
                    fe_htz.flag_neg2posLevelA[ind] = False
                    fe_htz.flag_neg2posLevelB[ind] = False
                    CTRL.ell = (fe_htz.psi_2_max[ind] - fe_htz.psi_2_min[ind]) * 0.5
                    fe_htz.u_offset[0] = (fe_htz.psi_2_max[0] + fe_htz.psi_2_min[0]) * 0.5 * fe_htz.u_offset_correction_factor
                    fe_htz.u_offset[1] = (fe_htz.psi_2_max[1] + fe_htz.psi_2_min[1]) * 0.5 * fe_htz.u_offset_correction_factor

                    fe_htz.psi_1_min[ind] = 0.0
                    fe_htz.psi_2_min[ind] = 0.0

                fe_htz.flag_pos2negLevelB[ind] = True
                if fe_htz.flag_pos2negLevelB[ind] == True: # 寻找磁链最小值
                    if fe_htz.psi_2[ind] < fe_htz.psi_2_min[ind]:
                        fe_htz.psi_2_min[ind] = fe_htz.psi_2[ind]

            else: # 磁链还没有变负，说明是虚假过零，比如在震荡，fe_htz.psi_2[0]>0
                fe_htz.flag_pos2negLevelA[ind] = False # /* 震荡的话，另一方的检测就有可能被触动？ */

        if fe_htz.psi_2_prev[ind]>0 and fe_htz.psi_2[ind]<0: # 发现磁链由正变负的时刻
            fe_htz.flag_pos2negLevelA[ind] = True
            fe_htz.time_pos2neg[ind] = CTRL.timebase

        if fe_htz.flag_neg2posLevelA[ind] == True:
            if fe_htz.psi_2_prev[ind]>0 and fe_htz.psi_2[ind]>0: # 二次检查，磁链已经是正的了
                if fe_htz.flag_neg2posLevelB[ind] == False:
                    fe_htz.count_positive_cycle+=1 # fe_htz.count_negative_cycle = 0
                    # 第一次进入寻找最大值的levelB，说明最小值已经检测到。
                    fe_htz.psi_1_min[ind] = fe_htz.psi_2_min[ind] # 不区别定转子磁链，区别：psi_2是连续更新的，而psi_1是离散更新的。
                    fe_htz.Delta_t_last = fe_htz.Delta_t
                    fe_htz.Delta_t = fe_htz.time_neg2pos[ind] - fe_htz.time_neg2pos_prev[ind]
                    fe_htz.time_neg2pos_prev[ind] = fe_htz.time_neg2pos[ind] # 备份作为下次耗时参考点
                    # 初始化
                    fe_htz.flag_pos2negLevelA[ind] = False
                    fe_htz.flag_pos2negLevelB[ind] = False

                    # fe_htz.accumulated__u_off_saturation_time_correction[ind] += fe_htz.u_off_saturation_time_correction[ind]
                    CTRL.ell = (fe_htz.psi_2_max[ind] - fe_htz.psi_2_min[ind]) * 0.5
                    fe_htz.u_offset[0] = (fe_htz.psi_2_max[0] + fe_htz.psi_2_min[0]) * 0.5 * fe_htz.u_offset_correction_factor
                    fe_htz.u_offset[1] = (fe_htz.psi_2_max[1] + fe_htz.psi_2_min[1]) * 0.5 * fe_htz.u_offset_correction_factor
                    fe_htz.psi_1_max[ind] = 0.0
                    fe_htz.psi_2_max[ind] = 0.0

                fe_htz.flag_neg2posLevelB[ind] = True 
                if fe_htz.flag_neg2posLevelB[ind] == True: # 寻找磁链最大值
                    if fe_htz.psi_2[ind] > fe_htz.psi_2_max[ind]:
                        fe_htz.psi_2_max[ind] = fe_htz.psi_2[ind]

            else: # 磁链还没有变正，说明是虚假过零，比如在震荡，fe_htz.psi_2[0]<0
                fe_htz.flag_neg2posLevelA[ind] = False

        if fe_htz.psi_2_prev[ind]<0 and fe_htz.psi_2[ind]>0: # 发现磁链由负变正的时刻
            fe_htz.flag_neg2posLevelA[ind] = True
            fe_htz.time_neg2pos[ind] = CTRL.timebase
    #由于原本的uoffset所得的结果使得谐波有抖动，于是想滤波但是没滤成
    
    # fe_htz.u_offset_filered[0] =_lpf(fe_htz.u_offset[0], fe_htz.u_offset_filered[0], 1, CTRL)
    # fe_htz.u_offset_filered[1] =_lpf(fe_htz.u_offset[1], fe_htz.u_offset_filered[1], 1, CTRL)
    # fe_htz.xFlux[2] = fe_htz.u_offset_filered[0]
    # fe_htz.xFlux[3] = fe_htz.u_offset_filered[1]
    # fe_htz.u_offset[0] = fe_htz.u_offset_filered[0]
    # fe_htz.u_offset[1] = fe_htz.u_offset_filered[1]
    
    fe_htz.xFlux[2] = fe_htz.u_offset[0]
    fe_htz.xFlux[3] = fe_htz.u_offset[1]

    fe_htz.psi_2_prev[0] = fe_htz.psi_2[0]
    fe_htz.psi_2_prev[1] = fe_htz.psi_2[1]

    # psi_2_ampl 在限幅前已经算过了，还有必要限幅后在这里再算一次吗？
    # fe_htz.psi_A_ampl = np.sqrt(fe_htz.psi_A[0] **2 + fe_htz.psi_A[1] **2)
    # if fe_htz.psi_A_ampl == 0:
    #     fe_htz.psi_A_ampl = 1.0
    # amplitude_inverse = 1.0 / fe_htz.psi_A_ampl
    fe_htz.psi_2_ampl = np.sqrt(fe_htz.psi_2[0] **2 + fe_htz.psi_2[1] **2)
    CTRL.flux_estimate_amplitude = fe_htz.psi_2_ampl
    if fe_htz.psi_2_ampl == 0:
        fe_htz.psi_2_ampl = 1.0
    amplitude_inverse = 1.0 / fe_htz.psi_2_ampl

    CTRL.cosT = fe_htz.psi_2[0] * amplitude_inverse
    CTRL.sinT = fe_htz.psi_2[1] * amplitude_inverse
    while ACM.theta_d> np.pi: ACM.theta_d -= 2*np.pi
    while ACM.theta_d<-np.pi: ACM.theta_d += 2*np.pi
    #使用编码器进入仅观测状态
    if CTRL.use_encoder_angle_no_matter_what == True:
        CTRL.theta_d = ACM.theta_d
        CTRL.cosT = np.cos(CTRL.theta_d)
        CTRL.sinT = np.sin(CTRL.theta_d)
'''speed observers'''
# According to CTRL.index_separate_speed_estimation, chose your speed observer
#0. encoder for speed
#1. SEPARATE_SPEED_OBSERVER
#2. NATRUE_SPEED_OBSERVER

def SEPARATE_SPEED_OBSERVER(CTRL, FE_param):
    RK4_ObserverSolver_CJH_Style(DYNAMICS_SpeedObserver, CTRL.xSpeed, CTRL.CL_TS, CTRL, FE_param)
    while CTRL.xSpeed[0]> np.pi: CTRL.xSpeed[0] -= 2*np.pi
    while CTRL.xSpeed[0]<-np.pi: CTRL.xSpeed[0] += 2*np.pi
    # CTRL.uab_prev[0] = CTRL.uab_curr[0] # This is needed only if voltage is measured, e.g., by eCAP. Remember to update the code below marked by [$].
    # CTRL.uab_prev[1] = CTRL.uab_curr[1] # This is needed only if voltage is measured, e.g., by eCAP. Remember to update the code below marked by [$].
    """ Speed Observer Outputs """
    CTRL.vartheta_d = CTRL.xSpeed[0]
    CTRL.omega_r_elec = CTRL.xSpeed[1]
    if CTRL.use_disturbance_feedforward_rejection == 0:
        CTRL.total_disrubance_feedforward = 0.0
    if CTRL.use_disturbance_feedforward_rejection == 1:
        CTRL.total_disrubance_feedforward = CTRL.xSpeed[2]
    elif CTRL.use_disturbance_feedforward_rejection == 2:
        CTRL.total_disrubance_feedforward = CTRL.xSpeed[2] + CTRL.ell2*CTRL.speed_observer_output_error
    
def NATRUE_SPEED_OBSERVER(CTRL, FE_param):
    CTRL.idq_c[0] = CTRL.iab_curr[0] * CTRL.cosT + CTRL.iab_curr[1] * CTRL.sinT
    CTRL.idq_c[1] = CTRL.iab_curr[0] *-CTRL.sinT + CTRL.iab_curr[1] * CTRL.cosT
    if CTRL.nsoaf_set_omega_ob != CTRL.nsoaf_omega_ob:
        CTRL.nsoaf_omega_ob = CTRL.nsoaf_set_omega_ob
        nso_one_parameter_tuning(CTRL)
    # init_nsoaf(CTRL)
    CTRL.nsoaf_xIq      = CTRL.nsoaf_xSpeed[0]
    CTRL.nsoaf_xOmg     = CTRL.nsoaf_xSpeed[1]
    CTRL.nsoaf_xTL      = CTRL.nsoaf_xSpeed[2]
    CTRL.nsoaf_uQ = CTRL.udq[1]
    # CTRL.uQ_now_filtered = _lpf(CTRL.nsoaf_uQ, CTRL.uQ_now_filtered, 30, CTRL)
    CTRL.nsoaf_output_error = CTRL.idq_c[1] - CTRL.nsoaf_xIq
    
    CTRL.active_power_real =  1* CTRL.idq_c[1]
    CTRL.active_power_est =  1* CTRL.nsoaf_xIq
    CTRL.active_power_error =  1* CTRL.nsoaf_output_error
    CTRL.KA = CTRL.KE + (CTRL.Ld - CTRL.Lq ) * CTRL.idq_c[0]
    CTRL.nsoaf_xTem = CTRL.CLARKE_TRANS_TORQUE_GAIN * CTRL.npp * CTRL.KA * CTRL.nsoaf_xIq
    RK4_ObserverSolver_CJH_Style(NSO_Dynamics,CTRL.nsoaf_xSpeed,CTRL.CL_TS, CTRL, FE_param)
    
    """ Speed Observer Outputs """
    # CTRL.omega_r_elec = ACM.omega_r_elec 
    CTRL.omega_r_elec = CTRL.nsoaf_xOmg 
    
############################################# MACHINE SIMULATION SECTION
def DYNAMICS_MACHINE(t, x, ACM, CLARKE_TRANS_TORQUE_GAIN=1.5):
    fx = np.zeros(ACM.NS) # s x = f(x)

    # theta_d_mech = x[0]
    # omega_r_mech = x[1]
    KA    = x[2]
    iD    = x[3]
    iQ    = x[4]
    # ACM.theta_d = x[0]*ACM.npp
    # ACM.omega_r = x[1]*ACM.npp
    if KA==0.0:
        ACM.omega_slip = 0.0
    else:
        ACM.omega_slip = ACM.Rreq * iQ / KA
    ACM.omega_syn  = x[1]*ACM.npp + ACM.omega_slip

    # 电磁子系统 (KA, iD, iQ as x[2], x[3], x[4])
    if ACM.Rreq > 0:
        # s KA
        fx[2] = ACM.Rreq*iD - ACM.Rreq / (ACM.Ld - ACM.Lq) * KA # [Apply Park Transorm to (31b)]
        # s iD
        fx[3] = (ACM.udq[0] - ACM.R*iD + ACM.omega_syn*ACM.Lq*iQ - fx[2]) / ACM.Lq # (6a)
    elif ACM.Rreq < 0:
        raise Exception('ACM.Rreq is used to calculate slip so it must be zero for PMSM.')
    else:
            # note fx[3] * ACM.Lq = ACM.udq[0] - ACM.R*iD + omega*ACM.Lq*iQ - fx[2]
            #  =>  fx[3] * ACM.Lq = ACM.udq[0] - ACM.R*iD + omega*ACM.Lq*iQ - (ACM.Ld - ACM.Lq) * fx[3] - 0.0
            #  =>  fx[3] * ACM.Ld = ACM.udq[0] - ACM.R*iD + omega*ACM.Lq*iQ
            #  =>  s iD
        # s iD
        fx[3] = (ACM.udq[0] - ACM.R*iD + ACM.omega_syn*ACM.Lq*iQ) / ACM.Ld
        # s KA
        fx[2] = (ACM.Ld - ACM.Lq) * fx[3] + 0.0
    # s iQ
    fx[4] = (ACM.udq[1] - ACM.R*iQ - ACM.omega_syn*ACM.Lq*iD - ACM.omega_syn*ACM.KA) / ACM.Lq

    # 机械子系统 (theta_d_mech, omega_mech as x[0], x[1])
    ACM.Tem = CLARKE_TRANS_TORQUE_GAIN * ACM.npp * KA * iQ # 电磁转矩计算
    fx[0] = x[1] + ACM.omega_slip / ACM.npp # mech. angular rotor position (accumulated)
    fx[1] = (ACM.Tem - ACM.TLoad) / ACM.Js  # mech. angular rotor speed

    return fx

def RK4_MACHINE(t, ACM, hs): # 四阶龙格库塔法
    k1, k2, k3, k4 = np.zeros(ACM.NS), np.zeros(ACM.NS), np.zeros(ACM.NS), np.zeros(ACM.NS) # incrementals at 4 stages
    xk, fx = np.zeros(ACM.NS), np.zeros(ACM.NS) # state x for stage 2/3/4, state derivative

    if False:
        """ this is about twice slower than loop through the element one by one """ 
        fx = DYNAMICS_MACHINE(t, ACM.x, ACM) # @t
        k1 = fx * hs
        xk = ACM.x + k1*0.5

        fx = DYNAMICS_MACHINE(t, xk, ACM)  # @t+hs/2
        k2 = fx * hs
        xk = ACM.x + k2*0.5

        fx = DYNAMICS_MACHINE(t, xk, ACM)  # @t+hs/2
        k3 = fx * hs
        xk = ACM.x + k3

        fx = DYNAMICS_MACHINE(t, xk, ACM)  # @t+hs
        k4 = fx * hs
        ACM.x = ACM.x + (k1 + 2*(k2 + k3) + k4)/6.0
    else:
        fx = DYNAMICS_MACHINE(t, ACM.x, ACM) # @t
        for i in range(ACM.NS):
            k1[i] = fx[i] * hs
            xk[i] = ACM.x[i] + k1[i]*0.5

        fx = DYNAMICS_MACHINE(t, xk, ACM)  # @t+hs/2
        for i in range(ACM.NS):
            k2[i] = fx[i] * hs
            xk[i] = ACM.x[i] + k2[i]*0.5

        fx = DYNAMICS_MACHINE(t, xk, ACM)  # @t+hs/2
        for i in range(ACM.NS):
            k3[i] = fx[i] * hs
            xk[i] = ACM.x[i] + k3[i]

        fx = DYNAMICS_MACHINE(t, xk, ACM)  # @t+hs
        for i in range(ACM.NS):
            k4[i] = fx[i] * hs
            # ACM.x_dot[i] = (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0 / hs # derivatives
            ACM.x[i] = ACM.x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0

############################################# BASIC FOC SECTION
def incremental_pi(reg):
    reg.Err = reg.setpoint - reg.measurement
    reg.Out = reg.OutPrev + \
        reg.Kp * (reg.Err - reg.ErrPrev) + \
        reg.Ki * reg.Err
    if reg.Out >    reg.OutLimit:
        reg.Out =   reg.OutLimit
    elif reg.Out < -reg.OutLimit:
        reg.Out =  -reg.OutLimit
    reg.ErrPrev = reg.Err
    reg.OutPrev = reg.Out

def tustin_pid(reg):

    # Error signal
    error = reg.setpoint - reg.measurement

    # Proportional
    proportional = reg.Kp * error

    # Integral
    reg.integrator = reg.integrator + 0.5 * reg.Ki * reg.T * (error + reg.prevError) # Tustin
    # reg.integrator = reg.integrator + reg.Ki * reg.T * (error) # Euler

    # Anti-wind-up via integrator clamping */
    if reg.integrator  >  reg.IntLimit:
        reg.integrator =  reg.IntLimit
    elif reg.integrator< -reg.IntLimit:
        reg.integrator = -reg.IntLimit

    # Derivative (band-limited differentiator) # Note: derivative on measurement, therefore minus sign in front of equation! */
    reg.differentiator = -(2.0 * reg.Kd * (reg.measurement - reg.prevMeasurement) \
                        + (2.0 * reg.tau - reg.T) * reg.differentiator) \
                        / (2.0 * reg.tau + reg.T)

    # Compute output and apply limits
    reg.Out = proportional + reg.integrator + reg.differentiator

    if reg.Out  >  reg.OutLimit:
        reg.Out =  reg.OutLimit
    elif reg.Out< -reg.OutLimit:
        reg.Out = -reg.OutLimit

    # Store error and measurement for later use */
    reg.prevError       = error
    reg.prevMeasurement = reg.measurement

    # Implement dynamic clamping
    reg.IntLimit = reg.OutLimit - proportional 

    # Return controller output */
    return reg.Out

def FOC(CTRL, reg_speed, reg_id, reg_iq):
    # speed loop
    reg_speed.setpoint = CTRL.cmd_rpm / 60 * 2*np.pi * CTRL.npp # [elec.rad/s]
    reg_speed.measurement = CTRL.omega_r_elec # [elec.rad/s]
    CTRL.velocity_loop_counter += 1
    if CTRL.velocity_loop_counter >= CTRL.velocity_loop_ceiling:
        CTRL.velocity_loop_counter = 0
        # incremental_pi(reg_speed)
        tustin_pid(reg_speed)

    # dq-frame current commands
    # Q-axis current command
    if CTRL.bool_apply_speed_closed_loop_control == True:
        CTRL.cmd_idq[1] = reg_speed.Out
        # CTRL.cmd_idq[0] = 0.0 # for user specifying
    # D-axis current command
    if CTRL.Rreq>0: # IM
        # psi command divided by magntizing inductance
        CTRL.cmd_idq[0] = CTRL.cmd_psi / (CTRL.Ld - CTRL.Lq) # [Wb] / [H]
        # slip angular speed
        CTRL.omega_slip = CTRL.Rreq * CTRL.cmd_idq[1] / CTRL.KA # Use commands for calculation (base off Harnefors recommendations)
    else: # PMSM
        CTRL.omega_slip = 0.0

        # if CTRL.bool_zero_id_control == True:
        #     CTRL.cmd_idq[0] = 0
        # else:
        #     # Field weakening control (simple)
        #         # 电气转折速度 = CTRL.DC_BUS_VOLTAGE/1.732 / CTRL.KA * 0.7
        #         # print('Turning point', 电气转折速度 * 60 / (2*np.pi * CTRL.npp), 'r/min')
        #         # if CTRL.omega_r_elec > 电气转折速度:
        #         #     # 计算弱磁电流
        #         #     下一弱磁速度增量 = 10 / 60 * 2*np.pi * CTRL.npp
        #         #     CTRL.cmd_idq[0] = - (CTRL.KE - 电气转折速度 / (CTRL.omega_r_elec + 下一弱磁速度增量)) / (CTRL.Ld - CTRL.Lq)
        #         #     # 修改速度环输出限幅
        #         #     reg_speed.OutLimit = np.sqrt((CTRL.IN*1.414)**2 - CTRL.cmd_idq[0]**2)
        #         #     print(f'{reg_speed.OutLimit=}')
        #     当前速度 = CTRL.omega_r_elec*60/(2*np.pi*CTRL.npp)
        #     MAX_DEMAG_CURRENT = 60
        #     if 当前速度 < 450:
        #         CTRL.cmd_idq[0] = 0
        #     elif 当前速度 < 1000:
        #         CTRL.cmd_idq[0] = (当前速度 - 450) / (1000 - 450) * -MAX_DEMAG_CURRENT
        #     else:
        #         CTRL.cmd_idq[0] = -MAX_DEMAG_CURRENT
        #     if CTRL.IN*1.414 > CTRL.cmd_idq[0]:
        #         reg_speed.OutLimit = np.sqrt((CTRL.IN*1.414)**2 - CTRL.cmd_idq[0]**2)

        if CTRL.bool_zero_id_control == True:
            CTRL.cmd_idq[0] = 0
        else:
            # Field weakening control (simple)
                # 电气转折速度 = CTRL.DC_BUS_VOLTAGE/1.732 / CTRL.KA * 0.7
                # print('Turning point', 电气转折速度 * 60 / (2*np.pi * CTRL.npp), 'r/min')
                # if CTRL.omega_r_elec > 电气转折速度:
                #     # 计算弱磁电流
                #     下一弱磁速度增量 = 10 / 60 * 2*np.pi * CTRL.npp
                #     CTRL.cmd_idq[0] = - (CTRL.KE - 电气转折速度 / (CTRL.omega_r_elec + 下一弱磁速度增量)) / (CTRL.Ld - CTRL.Lq)
                #     # 修改速度环输出限幅
                #     reg_speed.OutLimit = np.sqrt((CTRL.IN*1.414)**2 - CTRL.cmd_idq[0]**2)
                #     print(f'{reg_speed.OutLimit=}')
            当前速度 = CTRL.omega_r_elec*60/(2*np.pi*CTRL.npp)
            MAX_DEMAG_CURRENT = 60
            if 当前速度 < 450:
                CTRL.cmd_idq[0] = 0
            elif 当前速度 < 1000:
                CTRL.cmd_idq[0] = (当前速度 - 450) / (1000 - 450) * -MAX_DEMAG_CURRENT
            else:
                CTRL.cmd_idq[0] = -MAX_DEMAG_CURRENT
            if CTRL.IN*1.414 > CTRL.cmd_idq[0]:
                reg_speed.OutLimit = np.sqrt((CTRL.IN*1.414)**2 - CTRL.cmd_idq[0]**2)

    CTRL.omega_syn = CTRL.omega_r_elec + CTRL.omega_slip

    # d-axis
    reg_id.setpoint = CTRL.cmd_idq[0]
    reg_id.measurement = CTRL.idq[0]
    # incremental_pi(reg_id)
    tustin_pid(reg_id)
    CTRL.cmd_udq[0] = reg_id.Out

        # if HUMAN.use_disturbance_feedforward_rejection == 0:
        #     CTRL.cmd_idq[1] = reg_speed.Out
        # else:
        #     CTRL.cmd_idq[1] = HUMAN.KP*(reg_speed.setpoint-reg_speed.measurement) + OB.total_disrubance_feedforward

    # q-axis
    reg_iq.setpoint = CTRL.cmd_idq[1]
    reg_iq.measurement = CTRL.idq[1]
    # incremental_pi(reg_iq)
    tustin_pid(reg_iq)
    CTRL.cmd_udq[1] = reg_iq.Out

    # Decoupling between two axes of current loop controller
    if CTRL.bool_apply_decoupling_voltages_to_current_regulation:
        decoupled_M_axis_voltage = -CTRL.omega_syn *             CTRL.Lq * CTRL.cmd_idq[1]
        decoupled_T_axis_voltage =  CTRL.omega_syn * ( CTRL.KA + CTRL.Lq * CTRL.cmd_idq[0])
        CTRL.cmd_udq[0] += decoupled_M_axis_voltage
        CTRL.cmd_udq[1] += decoupled_T_axis_voltage
        # BUG: 这里的电压是不受限幅影响的啊哈哈哈哈哈，如果没有仿真SVPWM和逆变器，那么这边可以产生任意大的输出电压
        # BUG: 这里的电压是不受限幅影响的啊哈哈哈哈哈，如果没有仿真SVPWM和逆变器，那么这边可以产生任意大的输出电压
        # BUG: 这里的电压是不受限幅影响的啊哈哈哈哈哈，如果没有仿真SVPWM和逆变器，那么这边可以产生任意大的输出电压
        if CTRL.cmd_udq[0]   >  reg_iq.OutLimit:
            CTRL.cmd_udq[0]  =  reg_iq.OutLimit
        elif CTRL.cmd_udq[0] < -reg_iq.OutLimit:
            CTRL.cmd_udq[0]  = -reg_iq.OutLimit
        if CTRL.cmd_udq[1]   >  reg_iq.OutLimit:
            CTRL.cmd_udq[1]  =  reg_iq.OutLimit
        elif CTRL.cmd_udq[1] < -reg_iq.OutLimit:
            CTRL.cmd_udq[1]  = -reg_iq.OutLimit

def SFOC_Dynamic(CTRL, reg_speed, reg_id, reg_iq):
    pass

    CTRL.cmd_udq[1] = reg_iq.Out

############################################# DSP SECTION
def DSP(ACM, CTRL, reg_speed, reg_id, reg_iq, fe_htz, FE_param=1.0,ELL_param = 0.019):
    CTRL.timebase += CTRL.CL_TS

    """ Current Measurement """
    CTRL.iab[0] = CTRL.iab_curr[0] = ACM.iAlfa
    CTRL.iab[1] = CTRL.iab_curr[1] = ACM.iBeta

    """ Angular Position Measurement """
    if CTRL.index_voltage_model_flux_estimation == 0:
        CTRL.theta_d = ACM.theta_d
        # do this once per control interrupt
        CTRL.cosT = np.cos(CTRL.theta_d)
        CTRL.sinT = np.sin(CTRL.theta_d)

    # mixed holtz and boldea to estimate K_E
    elif CTRL.index_voltage_model_flux_estimation == 1:
        mixed_holtz_and_boldea_to_estimate_K_E(fe_htz, CTRL, ACM, FE_param)

    # B. Proposed Flux Command Error feedback PI Correction in Controller Frame (xRho) */
    elif CTRL.index_voltage_model_flux_estimation == 2:
        B_Proposed_Flux_Command_Error_feedback_PI_Correction(fe_htz, CTRL, ACM, FE_param)

    # chen saturarion time based
    elif CTRL.index_voltage_model_flux_estimation == 3:
        chen_saturation_time_based(fe_htz, CTRL, ACM, FE_param, ELL_param)

    # cancel saturation and estimate K_E and u_offset from \psi_max and \psi_min
    elif CTRL.index_voltage_model_flux_estimation == 4:
        no_saturation_to_cal_KE_and_uoffset(fe_htz, CTRL, ACM, FE_param)
    """ Park Transformation Essentials """
    # Park transformation
    CTRL.idq[0] = CTRL.iab[0] * CTRL.cosT + CTRL.iab[1] * CTRL.sinT
    CTRL.idq[1] = CTRL.iab[0] *-CTRL.sinT + CTRL.iab[1] * CTRL.cosT
    CTRL.udq[0] = CTRL.cmd_uab[0] * CTRL.cosT + CTRL.cmd_uab[1] * CTRL.sinT
    CTRL.udq[1] = CTRL.cmd_uab[0] *-CTRL.sinT + CTRL.cmd_uab[1] * CTRL.cosT
    # now we are ready to calculate torque using dq-currents
    CTRL.KA = (CTRL.Ld - CTRL.Lq) * CTRL.idq[0] + CTRL.KE # 有功磁链计算
    CTRL.Tem =     1.5 * CTRL.npp * CTRL.idq[1] * CTRL.KA # 电磁转矩计算

    """ Speed Estimation """
    if CTRL.index_separate_speed_estimation == 0:
        #TODO simulate the encoder
        CTRL.omega_r_elec = ACM.omega_r_elec
    elif CTRL.index_separate_speed_estimation == 1:
        SEPARATE_SPEED_OBSERVER(CTRL, FE_param)

        """ Nature Speed Observer """
    elif CTRL.index_separate_speed_estimation == 2:
        NATRUE_SPEED_OBSERVER(CTRL, FE_param)

    """ (Optional) Do Park transformation again using the position estimate from the speed observer """
    pass

    # update previous current measurement for soeed observation and flux estimation 
    CTRL.iab_prev[0] = CTRL.iab_curr[0]
    CTRL.iab_prev[1] = CTRL.iab_curr[1]

    """ Speed and Current Controller (two cascaded closed loops) """
    FOC(CTRL, reg_speed, reg_id, reg_iq)
    reverse_rotation(CTRL, ACM)
    # [$] Inverse Park transformation: get voltage commands in alpha-beta frame as SVPWM input
    CTRL.cmd_uab[0] = CTRL.cmd_udq[0] * CTRL.cosT + CTRL.cmd_udq[1] *-CTRL.sinT
    CTRL.cmd_uab[1] = CTRL.cmd_udq[0] * CTRL.sinT + CTRL.cmd_udq[1] * CTRL.cosT

############################################# Inverter and PWM
def SVGEN_DQ(v, one_over_Vdc):

    # Normalization (which converts [Volt] into [s])
    Talfa = v.Ualfa * one_over_Vdc # v.Ualfa is in sense of amplitude invariant Clarke transformation
    Tbeta = v.Ubeta * one_over_Vdc # v.Ubeta is in sense of amplitude invariant Clarke transformation
    Tz    = v.Unot  * one_over_Vdc # duration of the added zero sequence voltage

    # Inverse clarke transformation??
    A = Tbeta # 0 degree line pointing at 0 degree
    C =  1.7320508*Talfa - Tbeta # C =  sin( 60/180*np.pi)*Talfa - sin(30/180*np.pi)*Tbeta
    B = -1.7320508*Talfa - Tbeta # B = -sin( 60/180*np.pi)*Talfa - sin(30/180*np.pi)*Tbeta

    # 60 degree Sector determination
    Sector = 0 
    if (A > 0): Sector = 1
    if (C > 0): Sector = Sector+2
    if (B > 0): Sector = Sector+4

    # X,Y,Z calculations (Note an additional factor of 1.7320508 is introduced to be equivalent to normalizing Ualfa and Ubeta to a base value of Vdc/sqrt(3))
    XXX =              Tbeta*1.7320508
    YYY =  1.5*Talfa + Tbeta*0.8660254
    ZZZ = -1.5*Talfa + Tbeta*0.8660254

    if Sector == 0: # Sector 0: this is special case for (Ualfa,Ubeta) = (0,0)*/
        v.Ta = 0.5
        v.Tb = 0.5
        v.Tc = 0.5
    if Sector == 1: #Sector 1: t1=Z and t2=Y (abc ---> Tb,Ta,Tc)*/
        t1 = ZZZ
        t2 = YYY
        v.Tb=(1-t1-t2)*0.5 + Tz*0.5
        v.Ta = v.Tb+t1              # taon = tbon+t1        */
        v.Tc = v.Ta+t2              # tcon = taon+t2        */
    elif Sector == 2:     # Sector 2: t1=Y and t2=-X (abc ---> Ta,Tc,Tb)*/
        t1 = YYY
        t2 = -XXX
        v.Ta=(1-t1-t2)*0.5 + Tz*0.5
        v.Tc = v.Ta+t1              #  tcon = taon+t1       */
        v.Tb = v.Tc+t2              #  tbon = tcon+t2       */
    elif Sector == 3:     # Sector 3: t1=-Z and t2=X (abc ---> Ta,Tb,Tc)*/
        t1 = -ZZZ
        t2 = XXX
        v.Ta=(1-t1-t2)*0.5 + Tz*0.5
        v.Tb = v.Ta+t1              #   tbon = taon+t1      */
        v.Tc = v.Tb+t2              #   tcon = tbon+t2      */
    elif Sector == 4:     # Sector 4: t1=-X and t2=Z (abc ---> Tc,Tb,Ta)*/
        t1 = -XXX
        t2 = ZZZ
        v.Tc=(1-t1-t2)*0.5 + Tz*0.5
        v.Tb = v.Tc+t1              #   tbon = tcon+t1      */
        v.Ta = v.Tb+t2              #   taon = tbon+t2      */
    elif Sector ==  5:    # Sector 5: t1=X and t2=-Y (abc ---> Tb,Tc,Ta)*/
        t1 = XXX
        t2 = -YYY                   #   tbon = (1-t1-t2)*0.5    */
        v.Tb=(1-t1-t2)*0.5 + Tz*0.5
        v.Tc = v.Tb+t1              #   taon = tcon+t2      */
        v.Ta = v.Tc+t2
    elif Sector == 6:     # Sector 6: t1=-Y and t2=-Z (abc ---> Tc,Ta,Tb)*/
        t1 = -YYY
        t2 = -ZZZ
        v.Tc=(1-t1-t2)*0.5 + Tz*0.5
        v.Ta = v.Tc+t1              #   taon = tcon+t1      */
        v.Tb = v.Ta+t2              #   tbon = taon+t2      */

    # 高低有效逻辑翻转
    v.Ta = 1-v.Ta
    v.Tb = 1-v.Tb
    v.Tc = 1-v.Tc

    # 考虑到输出功率时母线电压会跌落，不要用满占空比。
    if (v.Ta>v.SYSTEM_MAX_PWM_DUTY_LIMATATION): v.Ta=v.SYSTEM_MAX_PWM_DUTY_LIMATATION
    if (v.Tb>v.SYSTEM_MAX_PWM_DUTY_LIMATATION): v.Tb=v.SYSTEM_MAX_PWM_DUTY_LIMATATION
    if (v.Tc>v.SYSTEM_MAX_PWM_DUTY_LIMATATION): v.Tc=v.SYSTEM_MAX_PWM_DUTY_LIMATATION
    if (v.Ta<v.SYSTEM_MIN_PWM_DUTY_LIMATATION): v.Ta=v.SYSTEM_MIN_PWM_DUTY_LIMATATION
    if (v.Tb<v.SYSTEM_MIN_PWM_DUTY_LIMATATION): v.Tb=v.SYSTEM_MIN_PWM_DUTY_LIMATATION
    if (v.Tc<v.SYSTEM_MIN_PWM_DUTY_LIMATATION): v.Tc=v.SYSTEM_MIN_PWM_DUTY_LIMATATION

    return v

def gate_signal_generator(ii, v, CPU_TICK_PER_SAMPLING_PERIOD, DEAD_TIME_AS_COUNT):
    # 波谷中断 # if ii % CPU_TICK_PER_SAMPLING_PERIOD == 0:
    if v.bool_interupt_event:
        v.bool_interupt_event = False # this clause is one-time-execution code
        v.bool_counting_down = False # counting up first
        v.carrier_counter = 0 # reset main counter

        # dead time
        v.deadtime_counter[0] = 0
        v.deadtime_counter[1] = 0
        v.deadtime_counter[2] = 0
        v.bool_RisingEdgeDelay_is_active[0] = False
        v.bool_RisingEdgeDelay_is_active[1] = False
        v.bool_RisingEdgeDelay_is_active[2] = False
        v.bool_FallingEdgeDelay_is_active[0] = False
        v.bool_FallingEdgeDelay_is_active[1] = False
        v.bool_FallingEdgeDelay_is_active[2] = False

    # 波峰中断 # if ii % CPU_TICK_PER_SAMPLING_PERIOD == CPU_TICK_PER_SAMPLING_PERIOD * 0.5:
    if v.carrier_counter == CPU_TICK_PER_SAMPLING_PERIOD * 0.5:
        v.bool_counting_down = True

        # dead time
        v.deadtime_counter[0] = 0
        v.deadtime_counter[1] = 0
        v.deadtime_counter[2] = 0
        v.bool_RisingEdgeDelay_is_active[0] = False
        v.bool_RisingEdgeDelay_is_active[1] = False
        v.bool_RisingEdgeDelay_is_active[2] = False
        v.bool_FallingEdgeDelay_is_active[0] = False
        v.bool_FallingEdgeDelay_is_active[1] = False
        v.bool_FallingEdgeDelay_is_active[2] = False

    # 计数
    if v.bool_counting_down:
        v.carrier_counter -= 1
    else:
        v.carrier_counter += 1

    # 理想门极信号
    v.S1 = v.phase_U_gate_signal = True if v.carrier_counter >= v.EPwm1Regs_CMPA_bit_CMPA else False
    v.S2 = v.phase_V_gate_signal = True if v.carrier_counter >= v.EPwm2Regs_CMPA_bit_CMPA else False
    v.S3 = v.phase_W_gate_signal = True if v.carrier_counter >= v.EPwm3Regs_CMPA_bit_CMPA else False

    v.S4, v.S5, v.S6 = not v.S1, not v.S2, not v.S3

    # 应用死区时间，获得实际门极信号
    # Insert dead time based on Active Hgih Complementary (AHC)
    if v.bool_counting_down == False:

        if v.carrier_counter >= v.EPwm1Regs_CMPA_bit_CMPA:
            v.deadtime_counter[0] += 1
            if v.deadtime_counter[0] <= DEAD_TIME_AS_COUNT:
                v.bool_RisingEdgeDelay_is_active[0] = True # this boolean variable is not used
                v.S1 = False
            else:
                pass # False
        if v.carrier_counter >= v.EPwm2Regs_CMPA_bit_CMPA:
            v.deadtime_counter[1] += 1
            if v.deadtime_counter[1] <= DEAD_TIME_AS_COUNT:
                v.bool_RisingEdgeDelay_is_active[1] = True # this boolean variable is not used
                v.S2 = False
            else:
                pass # False
        if v.carrier_counter >= v.EPwm3Regs_CMPA_bit_CMPA:
            v.deadtime_counter[2] += 1
            if v.deadtime_counter[2] <= DEAD_TIME_AS_COUNT:
                v.bool_RisingEdgeDelay_is_active[2] = True # this boolean variable is not used
                v.S3 = False
            else:
                pass # False
    elif v.bool_counting_down == True:

        if v.carrier_counter < v.EPwm1Regs_CMPA_bit_CMPA:
            v.deadtime_counter[0] += 1
            if v.deadtime_counter[0] < DEAD_TIME_AS_COUNT:
                v.bool_FallingEdgeDelay_is_active[0] = True # this boolean variable is not used
                v.S4 = False
            else:
                pass # False
        if v.carrier_counter < v.EPwm2Regs_CMPA_bit_CMPA:
            v.deadtime_counter[1] += 1
            if v.deadtime_counter[1] < DEAD_TIME_AS_COUNT:
                v.bool_FallingEdgeDelay_is_active[1] = True # this boolean variable is not used
                v.S5 = False
            else:
                pass # False
        if v.carrier_counter < v.EPwm3Regs_CMPA_bit_CMPA:
            v.deadtime_counter[2] += 1
            if v.deadtime_counter[2] < DEAD_TIME_AS_COUNT:
                v.bool_FallingEdgeDelay_is_active[2] = True # this boolean variable is not used
                v.S6 = False
            else:
                pass # False

# 让电机正反转
def reverse_rotation(CTRL, ACM):
    if CTRL.bool_reverse_rotation == True:
            if ACM.theta_d > 2.5 and CTRL.flag_reverse_rotation == True:
                CTRL.counter_rotation = CTRL.counter_rotation + 1
                CTRL.flag_reverse_rotation = False
            elif ACM.theta_d < - 2.5 and CTRL.flag_reverse_rotation == False:
                CTRL.flag_reverse_rotation = True
            if CTRL.counter_rotation == 2 * ACM.npp:
                if CTRL.bool_apply_speed_closed_loop_control == True:
                    CTRL.counter_rotation = 0
                    CTRL.cmd_rpm  = -1 * CTRL.cmd_rpm
                else:
                    CTRL.counter_rotation = 0
                    CTRL.cmd_idq[1] = -1 * CTRL.cmd_idq[1]
############################################# Wrapper level 1 (Main simulation | Incremental Edition)
""" MAIN for  ('-time simulation """
def vehicel_load_model(t, ACM):
    EVM=1500    #####(车身质量)
    EVA=2.5     ####(迎风面积)
    EVCD=0.37   #####(风阻系数)
    EVF=0.015   #####(摩擦系数)
    EVR=0.297   #####(车轮转动半径)
    grav=9.8    #####(重力加速度g)
    VEV=ACM.omega_r_mech*60*EVR*60*1e-3   ####车速
    # VEV=60               ##### 暂时车速
    FW=EVCD*EVA*VEV*VEV/21.15    ##### 风阻
    FF=EVM*grav*EVF              ##### 滚阻
    FLoad=(FW+FF)*0.5            ##### 单侧阻力负载
    ACM.TLoad=FLoad*EVR          ##### 单侧转矩负载
    ACM.Js = EVJ = EVM*EVR*EVR*0.25  ##### 单轮等效转动惯量

def ACMSimPyIncremental(t0, TIME, ACM=None, CTRL=None, reg_id=None, reg_iq=None, reg_speed=None, fe_htz=None, FE_param = 1.0, ELL_param = 0.019):

    # RK4 simulation and controller execution relative freuqencies
    MACHINE_TS = CTRL.CL_TS / ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    controller_down_sampling_ceiling = int(CTRL.CL_TS / MACHINE_TS)

    # SVPWM
    CPU_TICK_PER_SAMPLING_PERIOD = ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    DEAD_TIME_AS_COUNT = int(200*0.5e-4*CPU_TICK_PER_SAMPLING_PERIOD) # 200 count for 0--5000--0 counting sequence
    # print(t0, 's', 'DEAD_TIME_AS_COUNT =', DEAD_TIME_AS_COUNT, )
    Vdc = CTRL.DC_BUS_VOLTAGE # Vdc is assumed measured and known
    one_over_Vdc = 1/Vdc
    svgen1 = SVgen_Object(CPU_TICK_PER_SAMPLING_PERIOD)
    # print('Vdc, CPU_TICK_PER_SAMPLING_PERIOD, controller_down_sampling_ceiling', Vdc, CPU_TICK_PER_SAMPLING_PERIOD, controller_down_sampling_ceiling)

    # watch variabels
    machine_times = np.arange(t0, t0+TIME, MACHINE_TS)
    watch_data    = np.zeros( (60, len(machine_times)) ) # new
    # control_times = np.arange(t0, t0+TIME, CTRL.CL_TS)
    # watch_data = np.zeros( (40, len(control_times)) ) # old

    # Main loop
    jj = controller_down_sampling_ceiling # run controller at step 1
    watch_index = 0
    for ii in range(len(machine_times)):

        t = machine_times[ii]
        # print(t)

        """ Machine Simulation @ MACHINE_TS """
        # Numerical Integration (ode4) with 5 states
        if ACM.bool_apply_load_model: vehicel_load_model(t, ACM)
        RK4_MACHINE(t, ACM, hs=MACHINE_TS)

        """ Machine Simulation Output @ MACHINE_TS """
        # Generate output variables for easy access
        ACM.theta_d_mech = ACM.x[0] # ACM.x[0] = ACM.x[0] - ACM.x[0]//(2*np.pi)*(2*np.pi)
        ACM.omega_r_mech = ACM.x[1]
        ACM.KA           = ACM.x[2]
        ACM.iD           = ACM.x[3]
        ACM.iQ           = ACM.x[4]
        ACM.theta_d      = ACM.theta_d_mech * ACM.npp
        ACM.omega_r_elec = ACM.omega_r_mech * ACM.npp
        ACM.omega_syn    = ACM.omega_r_elec + ACM.omega_slip

        # Inverse Park transformation
        ACM.cosT = np.cos(ACM.theta_d)
        ACM.sinT = np.sin(ACM.theta_d)
        ACM.iAlfa = ACM.iD * ACM.cosT + ACM.iQ *-ACM.sinT # as motor controller input
        ACM.iBeta = ACM.iD * ACM.sinT + ACM.iQ * ACM.cosT # as motor controller input

        jj += 1
        if jj >= controller_down_sampling_ceiling:
            jj = 0

            humans_give_commands.humans_give_commands(CTRL,ACM,t)

            """ DSP @ CL_TS """
            # print(ii+1)
            DSP(ACM=ACM,
                CTRL=CTRL,
                reg_speed=reg_speed,
                reg_id=reg_id,
                reg_iq=reg_iq,
                fe_htz=fe_htz, 
                FE_param=FE_param,
                ELL_param=ELL_param)

            # DEBUG
            # CTRL.cmd_uab[0] = 10*np.cos(5*2*np.pi*CTRL.timebase)
            # CTRL.cmd_uab[1] = 10*np.sin(5*2*np.pi*CTRL.timebase)

            # SVPWM for voltage source inverter
            svgen1.Ualfa = CTRL.cmd_uab[0]
            svgen1.Ubeta = CTRL.cmd_uab[1]
            SVGEN_DQ(svgen1, one_over_Vdc)
            # 高低有效逻辑翻转（实物用的光耦HCP2060带了一个非门，导致PWM输出是低有效，SVGEN_DQ是和实物一致的，但是在仿真里我们得马上反回来，否则PWM输出就反相了）
            svgen1.Ta = 1-svgen1.Ta
            svgen1.Tb = 1-svgen1.Tb
            svgen1.Tc = 1-svgen1.Tc
            svgen1.EPwm1Regs_CMPA_bit_CMPA = (int)(svgen1.Ta*CPU_TICK_PER_SAMPLING_PERIOD*0.5) # 0.5 for up and down counting # 50000000*CTRL.CL_TS)
            svgen1.EPwm2Regs_CMPA_bit_CMPA = (int)(svgen1.Tb*CPU_TICK_PER_SAMPLING_PERIOD*0.5) # 0.5 for up and down counting # 50000000*CTRL.CL_TS)
            svgen1.EPwm3Regs_CMPA_bit_CMPA = (int)(svgen1.Tc*CPU_TICK_PER_SAMPLING_PERIOD*0.5) # 0.5 for up and down counting # 50000000*CTRL.CL_TS)

            svgen1.bool_interupt_event = True

        """ Voltage Source Inverter (in alpha-beta frame) """
        if CPU_TICK_PER_SAMPLING_PERIOD >= 20: # implementing SVPWM

            # Amplitude invariant Clarke transformation
            ACM.ia = ACM.iAlfa
            ACM.ib = ACM.iAlfa*-0.5 + ACM.iBeta*0.8660254
            ACM.ic = ACM.iAlfa*-0.5 + ACM.iBeta*-0.8660254

            # Get S1 -- S6
            gate_signal_generator(ii, svgen1, CPU_TICK_PER_SAMPLING_PERIOD=CPU_TICK_PER_SAMPLING_PERIOD, DEAD_TIME_AS_COUNT=DEAD_TIME_AS_COUNT)

            # 端电势
            # inverter connects motor terminals to dc bus capacitor depending on gate signals and phase current (during dead zone)
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

            # 线电压 AC 和 BC
            svgen1.line_to_line_voltage_AC = svgen1.voltage_potential_at_terminal[0] - svgen1.voltage_potential_at_terminal[2]
            svgen1.line_to_line_voltage_BC = svgen1.voltage_potential_at_terminal[1] - svgen1.voltage_potential_at_terminal[2]
            svgen1.line_to_line_voltage_AB = svgen1.voltage_potential_at_terminal[0] - svgen1.voltage_potential_at_terminal[1]

            # 线电压 做 Amplitude invariant Clarke transformation 获得 alpha-beta 电压
            ACM.uab[0] = svgen1.line_to_line_voltage_AC*0.6666667 - (svgen1.line_to_line_voltage_BC + 0)*0.3333333
            ACM.uab[1] = 0.577350269 * (svgen1.line_to_line_voltage_BC - 0)

        else:
            # (no SVPWM, the discrepancy between CTRL.cosT and ACM.cosT will be simulated, i.e., the zero-hold feature of the inverter)
            ACM.uab[0] = CTRL.cmd_uab[0]
            ACM.uab[1] = CTRL.cmd_uab[1]

        # Park transformation
        ACM.udq[0] = ACM.uab[0] *  ACM.cosT + ACM.uab[1] * ACM.sinT
        ACM.udq[1] = ACM.uab[0] * -ACM.sinT + ACM.uab[1] * ACM.cosT
        import collect_data
        
        watch_index = collect_data.collect_data(watch_data, watch_index, CTRL, ACM, reg_id, reg_iq, reg_speed, fe_htz)

    # return machine_times, watch_data # old
    return machine_times, watch_data # new

def lpf1_inverter(array):
    y_tminus1 = 0.0
    new_array = []
    for x in array:
        new_x = y_tminus1 + 5* 0.00020828993959591752 * (x - y_tminus1)
        y_tminus1 = new_x
        new_array.append(y_tminus1)
    return new_array


# %%

