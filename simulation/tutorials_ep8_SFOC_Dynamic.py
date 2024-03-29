# %%
############################################# PACKAGES
from numba.experimental import jitclass
from numba import njit, int32, float64
from pylab import np, plt, mpl
plt.style.use('ggplot')

############################################# CLASS DEFINITION 
@jitclass(
    spec=[
        # CONTROL
            # constants
            ('CL_TS', float64),
            ('VL_TS', float64),
            ('velocity_loop_counter', float64),
            ('velocity_loop_ceiling', float64),
            # feedback / input
            ('theta_d', float64),
            ('theta_m', float64),
            ('omega_r_elec', float64),
            ('omega_syn', float64),
            ('omega_slip', float64),
            ('uab', float64[:]),
            ('psi_stator_MT_fb', float64[:]),
            ('psi_stator_ab_fb', float64[:]),
            ('kPFL', float64),
            ('kPCL', float64),
            ('Intergral_of_iT_error', float64),
            # ('uab_prev', float64[:]),
            # ('uab_curr', float64[:]),
            ('iab', float64[:]),
            ('iab_prev', float64[:]),
            ('iab_curr', float64[:]),
            # states
            ('timebase', float64),
            ('Tem', float64),
            ('cosT', float64),
            ('sinT', float64),
            ('cosT_MT', float64),
            ('sinT_MT', float64),
            ('lastcomplex', float64[:]),
            # commands
            ('cmd_iMT',float64[:]),
            ('cmd_uMT', float64[:]),
            ('cmd_idq', float64[:]),
            ('cmd_udq', float64[:]),
            ('cmd_uab', float64[:]),
            ('cmd_rpm', float64),
            ('cmd_KA',  float64),
            ('cmd_psi_Ms', float64),
            ('index_separate_speed_estimation', int32),
            ('use_disturbance_feedforward_rejection', int32),
            ('bool_apply_decoupling_voltages_to_current_regulation', int32),
            ('bool_apply_speed_closed_loop_control', int32),
            ('bool_zero_id_control', int32),
            ('bool_use_FOC_or_SFOC', int32),
            # commands (sweep freuqency)
            ('bool_apply_sweeping_frequency_excitation', int32),
            ('bool_overwrite_speed_commands', int32),
            ('CMD_CURRENT_SINE_AMPERE', float64),
            ('CMD_SPEED_SINE_RPM', float64),
            ('CMD_SPEED_SINE_HZ', float64),
            ('CMD_SPEED_SINE_STEP_SIZE', float64),
            ('CMD_SPEED_SINE_LAST_END_TIME', float64),
            ('CMD_SPEED_SINE_END_TIME', float64),
            ('CMD_SPEED_SINE_HZ_CEILING', float64),
        # MOTOR
            # name plate data
            ('npp',   int32),
            ('IN',  float64),
            ('DC_BUS_VOLTAGE', float64),
            # electrical parameters
            ('R',   float64),
            ('Ld',  float64),
            ('Lq',  float64),
            ('KE',  float64),
            ('KA',  float64),
            ('Rreq',float64),
            ('Lsigma', float64),
            # mechanical parameters
            ('Js',  float64),
        # OBSERVER
            # feedback / inputs
            ('idq', float64[:]),
            ('iMT', float64[:]),
            ('emf_stator', float64[:]),
            # states
            ('NS', int32),
            ('xS', float64[:]),
            ('xT', float64[:]),
            # outputs
            ('speed_observer_output_error', float64),
            ('vartheta_d', float64),
            ('total_disrubance_feedforward', float64),
            # gains
            ('ell1', float64),
            ('ell2', float64),
            ('ell3', float64),
            ('ell4', float64),
            #
            ('one_over_six', float64),
    ])
class The_Motor_Controller:
    def __init__(self, CL_TS, VL_TS,
        init_npp = 4,
        init_IN = 3,
        init_R = 1.1,
        init_Ld = 5e-3,
        init_Lq = 6e-3,
        init_KE = 0.095,
        init_KA = 0.095,
        init_Rreq = -1, # note division by 0 is equal to infinity
        init_Js = 0.0006168,
        DC_BUS_VOLTAGE = 300,
    ):
        ''' CONTROL '''
        # constants
        self.CL_TS = CL_TS
        self.VL_TS = VL_TS
        self.velocity_loop_ceiling = VL_TS / CL_TS
        self.velocity_loop_counter = self.velocity_loop_ceiling - 1
        print('\tCTRL.velocity_loop_ceiling =', self.velocity_loop_ceiling)
        # feedback / input
        self.theta_d = 0.0
        self.theta_m = 0.0
        self.omega_r_elec = 0.0
        self.omega_syn = 0.0
        self.omega_slip = 0.0
        self.uab      = np.zeros(2, dtype=np.float64)
        self.psi_stator_MT_fb = np.zeros(2, dtype=np.float64)
        self.psi_stator_ab_fb = np.zeros(2, dtype=np.float64)
        self.kPFL = 2 # 10 
        self.kPCL = 20 # 200
        self.Intergral_of_iT_error = 0.0
        # self.uab_prev = np.zeros(2, dtype=np.float64)
        # self.uab_curr = np.zeros(2, dtype=np.float64)
        self.iab      = np.zeros(2, dtype=np.float64)
        self.iab_prev = np.zeros(2, dtype=np.float64)
        self.iab_curr = np.zeros(2, dtype=np.float64)
        # states
        self.timebase = 0.0
        self.Tem = 0.0
        self.cosT = 1.0
        self.sinT = 0.0
        self.lastcomplex = np.zeros(2, dtype=np.float64)
        # commands
        self.cmd_iMT = np.zeros(2, dtype=np.float64)
        self.cmd_uMT = np.zeros(2, dtype=np.float64)
        self.cmd_idq = np.zeros(2, dtype=np.float64)
        self.cmd_udq = np.zeros(2, dtype=np.float64)
        self.cmd_uab = np.zeros(2, dtype=np.float64)
        self.cmd_rpm = 0.0
        self.cmd_KA  = init_KA
        self.cmd_psi_Ms = init_KA
        self.index_separate_speed_estimation = 0
        self.use_disturbance_feedforward_rejection = 0
        self.bool_apply_decoupling_voltages_to_current_regulation = False
        self.bool_apply_speed_closed_loop_control = True
        self.bool_zero_id_control = True
        self.bool_use_FOC_or_SFOC = True
        # sweep frequency
        self.bool_apply_sweeping_frequency_excitation = False
        self.bool_overwrite_speed_commands = False
        self.CMD_CURRENT_SINE_AMPERE = 1 # [A]
        self.CMD_SPEED_SINE_RPM = 100 # [r/min]
        self.CMD_SPEED_SINE_HZ = 0 # [Hz]
        self.CMD_SPEED_SINE_STEP_SIZE = 1 # [Hz]
        self.CMD_SPEED_SINE_LAST_END_TIME = 0.0
        self.CMD_SPEED_SINE_END_TIME = 0.0
        self.CMD_SPEED_SINE_HZ_CEILING = 100
        ''' MOTOR '''
        self.npp  = init_npp
        self.IN   = init_IN
        self.R    = init_R
        self.Ld   = init_Ld
        self.Lq   = init_Lq
        self.KE   = init_KE
        self.KA   = init_KA
        self.Rreq = init_Rreq
        self.Js   = init_Js
        self.DC_BUS_VOLTAGE = DC_BUS_VOLTAGE
        self.Lsigma = 22*1e-3

        ''' OBSERVER '''
        # feedback / input
        self.idq = np.zeros(2, dtype=np.float64)
        self.iMT = np.zeros(2, dtype=np.float64)
        self.emf_stator = np.zeros(2, dtype=np.float64)
        # state
        self.NS   = 6 # = max(NS_SPEED, NS_FLUX)
        self.xS   = np.zeros(self.NS, dtype=np.float64) # the internal states of speed estimator
        self.xT   = np.zeros(self.NS, dtype=np.float64) # the internal states of torque estimator
        # outputs
        self.speed_observer_output_error = 0.0
        self.vartheta_d = 0.0
        self.total_disrubance_feedforward = 0.0

        # gains
        omega_ob = 100 # [rad/s]
        self.ell1 = 0.0
        self.ell2 = 0.0
        self.ell3 = 0.0
        self.ell4 = 0.0
        if False: # 2nd-order speed observer (assuming speed feedback)
            self.ell2 = 2 * omega_ob
            self.ell3 =     omega_ob**2 * init_Js/init_npp
        elif False: # 2nd-order position observer
            self.ell1 = 2 * omega_ob
            self.ell2 =     omega_ob**2 * init_Js/init_npp
        elif True: # 3rd-order position observer
            self.ell1 = 3 * omega_ob
            self.ell2 = 3 * omega_ob**2
            self.ell3 =     omega_ob**3 * init_Js/init_npp
        else: # 4th-order position observer
            self.ell1 = 4 * omega_ob
            self.ell2 = 6 * omega_ob**2
            self.ell3 = 4 * omega_ob**3 * init_Js/init_npp
            self.ell4 =     omega_ob**4

        self.one_over_six = 1.0 / 6.0

@jitclass(
    spec=[
        # name plate data
        ('npp',   int32),
        ('npp_inv', float64),
        ('IN',  float64),
        # electrical parameters
        ('R',   float64),
        ('Ld',  float64),
        ('Lq',  float64),
        ('KE',  float64),
        ('KA', float64),
        ('Rreq',float64),
        # mechanical parameters
        ('Js',  float64),
        ('Js_inv', float64),
        # states
        ('NS',    int32),
        ('x',   float64[:]),
        # inputs
        ('uab',   float64[:]),
        ('udq',   float64[:]),
        ('TLoad', float64),
        # output
        ('omega_slip', float64),
        ('omega_r_elec', float64),
        ('omega_r_mech', float64),
        ('omega_syn', float64),
        ('theta_d', float64),
        ('theta_d_mech', float64),
        ('iD', float64),
        ('iQ', float64),
        ('iAlfa', float64),
        ('iBeta', float64),
        ('ia', float64),
        ('ib', float64),
        ('ic', float64),
        ('Tem', float64),
        ('cosT', float64),
        ('sinT', float64),
        # simulation settings
        ('MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD', int32),
        ('bool_apply_load_model', int32)
    ])
class The_AC_Machine:
    def __init__(self, CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1):
        # name plate data
        self.npp = CTRL.npp
        self.npp_inv = 1.0/self.npp
        self.IN  = CTRL.IN
        # electrical parameters
        self.R   = CTRL.R
        self.Ld  = CTRL.Ld
        self.Lq  = CTRL.Lq
        self.Rreq  = CTRL.Rreq
        if self.Rreq >0:
            if CTRL.KE>0:
                raise Exception('KE should be zero for induction machine.')
        self.KA  = CTRL.KA
        # mechanical parameters
        self.Js  = CTRL.Js # kg.m^2
        self.Js_inv = 1.0/self.Js
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

@jitclass(
    spec=[
        ('Kp', float64),
        ('Ki', float64),
        ('Err', float64),
        ('setpoint', float64),
        ('measurement', float64),
        ('Out', float64),
        ('OutLimit', float64),
        ('ErrPrev', float64),
        ('OutPrev', float64),
    ])
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

@jitclass(
    spec=[
        ('Kp', float64),
        ('Ki', float64),
        ('Kd', float64),
        ('tau', float64),
        ('OutLimit', float64),
        ('IntLimit', float64),
        ('T', float64),
        ('integrator', float64),
        ('prevError', float64),
        ('differentiator', float64),
        ('prevMeasurement', float64),
        ('Out', float64),
        ('setpoint', float64),
        ('measurement', float64),
    ])
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
        self.measurement = 0.0;

@jitclass(
    spec=[
        ('Ualfa', float64),
        ('Ubeta', float64),
        ('Unot', float64),
        ('Ta', float64),
        ('Tb', float64),
        ('Tc', float64),
        ('SYSTEM_MAX_PWM_DUTY_LIMATATION', float64),
        ('SYSTEM_MIN_PWM_DUTY_LIMATATION', float64),
        # Those variables are only needed in simulation
        ('bool_interupt_event', float64),
        ('bool_counting_down', float64),
        ('bool_RisingEdgeDelay_is_active', float64[:]),
        ('bool_FallingEdgeDelay_is_active', float64[:]),
        ('carrier_counter', float64),
        ('deadtime_counter', float64[:]),
        ('S1', int32),
        ('S2', int32),
        ('S3', int32),
        ('S4', int32),
        ('S5', int32),
        ('S6', int32),
        ('EPwm1Regs_CMPA_bit_CMPA', float64),
        ('EPwm2Regs_CMPA_bit_CMPA', float64),
        ('EPwm3Regs_CMPA_bit_CMPA', float64),
        ('phase_U_gate_signal', float64),
        ('phase_V_gate_signal', float64),
        ('phase_W_gate_signal', float64),
        ('voltage_potential_at_terminal', float64[:]),
        ('line_to_line_voltage_AC', float64),
        ('line_to_line_voltage_BC', float64),
        ('line_to_line_voltage_AB', float64),
    ])
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

############################################# OBSERVERS SECTION
@njit(nogil=True)
def DYNAMICS_SpeedObserver(x, CTRL):
    fx = np.zeros(6)

    # [rad]
    # output_error = np.sin(CTRL.theta_d - x[0])
    output_error = angle_diff(CTRL.theta_d, x[0]) # OE version 2
        # CTRL.output_error = np.sin(CTRL.theta_d - CTRL.xS[0]) # OE version 1 simple and silly
        # CTRL.output_error = angle_diff(CTRL.theta_d - CTRL.xS[0]) # OE version 2
        # CTRL.output_error = q-axis component # OE version 3 Boldea
    CTRL.speed_observer_output_error = output_error

    # 机械子系统 (omega_r_elec, theta_d, theta_r_mech)
    fx[0] = CTRL.ell1*output_error + x[1]
    fx[1] = CTRL.ell2*output_error + (CTRL.Tem + x[2]) * CTRL.npp/CTRL.Js # elec. angular rotor speed
    fx[2] = CTRL.ell3*output_error + x[3]
    fx[3] = CTRL.ell4*output_error + 0.0
    return fx

@njit(nogil=True)
def RK4_ObserverSolver_CJH_Style(THE_DYNAMICS, x, hs, CTRL):
    NS = CTRL.NS # THIS SHOULD BE A CONSTANT THROUGHOUT THE CODES!!!
    k1, k2, k3, k4 = np.zeros(NS), np.zeros(NS), np.zeros(NS), np.zeros(NS) # incrementals at 4 stages
    xk, fx = np.zeros(NS), np.zeros(NS) # state x for stage 2/3/4, state derivative

    CTRL.uab[0] = CTRL.cmd_uab[0]
    CTRL.uab[1] = CTRL.cmd_uab[1]
    CTRL.iab[0] = CTRL.iab_prev[0]
    CTRL.iab[1] = CTRL.iab_prev[1]
    fx = THE_DYNAMICS(x, CTRL)
    for i in range(0, NS):
        k1[i] = fx[i] * hs
        xk[i] = x[i] + k1[i]*0.5

    CTRL.iab[0] = 0.5*(CTRL.iab_prev[0]+CTRL.iab_curr[0])
    CTRL.iab[1] = 0.5*(CTRL.iab_prev[1]+CTRL.iab_curr[1])
    fx = THE_DYNAMICS(xk, CTRL)
    for i in range(0, NS):
        k2[i] = fx[i] * hs
        xk[i] = x[i] + k2[i]*0.5

    fx = THE_DYNAMICS(xk, CTRL)
    for i in range(0, NS):
        k3[i] = fx[i] * hs
        xk[i] = x[i] + k3[i]

    CTRL.iab[0] = CTRL.iab_curr[0]
    CTRL.iab[1] = CTRL.iab_curr[1]
    fx = THE_DYNAMICS(xk, CTRL)
    for i in range(0, NS):
        k4[i] = fx[i] * hs
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i]) * CTRL.one_over_six

@njit(nogil=True)
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

############################################# MACHINE SIMULATION SECTION
@njit(nogil=True)
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

@njit(nogil=True)
def RK4_MACHINE(t, ACM, hs): # 四阶龙格库塔法
    NS = ACM.NS
    k1, k2, k3, k4 = np.zeros(NS), np.zeros(NS), np.zeros(NS), np.zeros(NS) # incrementals at 4 stages
    xk, fx = np.zeros(NS), np.zeros(NS) # state x for stage 2/3/4, state derivative

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
        for i in range(NS):
            k1[i] = fx[i] * hs
            xk[i] = ACM.x[i] + k1[i]*0.5

        fx = DYNAMICS_MACHINE(t, xk, ACM)  # @t+hs/2
        for i in range(NS):
            k2[i] = fx[i] * hs
            xk[i] = ACM.x[i] + k2[i]*0.5

        fx = DYNAMICS_MACHINE(t, xk, ACM)  # @t+hs/2
        for i in range(NS):
            k3[i] = fx[i] * hs
            xk[i] = ACM.x[i] + k3[i]

        fx = DYNAMICS_MACHINE(t, xk, ACM)  # @t+hs
        for i in range(NS):
            k4[i] = fx[i] * hs
            # ACM.x_dot[i] = (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0 / hs # derivatives
            ACM.x[i] = ACM.x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0

############################################# BASIC FOC SECTION
@njit(nogil=True)
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

@njit(nogil=True)
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

    # Return controller output */
    return reg.Out

@njit(nogil=True)
def FOC(CTRL, reg_speed, reg_id, reg_iq):

    """ Park Transformation Essentials """
    # do this once per control interrupt
    CTRL.cosT = np.cos(CTRL.theta_d)
    CTRL.sinT = np.sin(CTRL.theta_d)
    # Park transformation
    CTRL.idq[0] = CTRL.iab[0] * CTRL.cosT + CTRL.iab[1] * CTRL.sinT
    CTRL.idq[1] = CTRL.iab[0] *-CTRL.sinT + CTRL.iab[1] * CTRL.cosT

    # now we are ready to calculate torque using dq-currents
    CTRL.KA = (CTRL.Ld - CTRL.Lq) * CTRL.idq[0] + CTRL.KE # 有功磁链计算
    CTRL.Tem =     1.5 * CTRL.npp * CTRL.idq[1] * CTRL.KA # 电磁转矩计算

    # Speed regulation (FOC)
    reg_speed.setpoint = CTRL.cmd_rpm / 60 * 2*np.pi * CTRL.npp # [elec.rad/s]
    reg_speed.measurement = CTRL.omega_r_elec # [elec.rad/s]
    CTRL.velocity_loop_counter += 1
    if CTRL.velocity_loop_counter >= CTRL.velocity_loop_ceiling:
        CTRL.velocity_loop_counter = 0
        # incremental_pi(reg_speed)
        tustin_pid(reg_speed)

    # dq-frame current commands
    if CTRL.bool_apply_speed_closed_loop_control == True:
        CTRL.cmd_idq[1] = reg_speed.Out
        # CTRL.cmd_idq[0] = 0.0 # for user specifying

    # slip and syn frequencies
    if CTRL.Rreq>0: # IM
        CTRL.cmd_idq[0] = CTRL.cmd_KA / (CTRL.Ld - CTRL.Lq) # [Wb] / [H]
        CTRL.omega_slip = CTRL.Rreq * CTRL.cmd_idq[1] / CTRL.KA # Use commands for calculation (base off Harnefors recommendations)
    else:
        CTRL.omega_slip = 0.0

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

# @njit(nogil=True)
# def get_angular_frequency_of_psi_ab(CTRL):
#     x = CTRL.psi_stator_ab_fb[0]
#     y = CTRL.psi_stator_ab_fb[1]
#     # xdot = (x - CTRL.lastcomplex[0])/CTRL.CLTS #需要低通滤波
#     # ydot = (y - CTRL.lastcomplex[1])/CTRL.CLTS
#     xdot = CTRL.cmd_uab[0] - CTRL.iab[0]*CTRL.R 
#     ydot = CTRL.cmd_uab[1]] - CTRL.iab[1]*CTRL.R
#     den = np.sqrt(x**2 - y**2)
#     CTRL.omega_syn = (ydot*x - xdot*y) /den

#     # CTRL.lastcomplex[0] = CTRL.psi_stator_ab_fb[0]
#     # CTRL.lastcomplex[1] = CTRL.psi_stator_ab_fb[1]

@njit(nogil=True)
def SFOC_Dynamic(CTRL, reg_speed, reg_id, reg_iq):

    """ Park Transformation Essentials """
    # Park transformation
    CTRL.iMT[0] = CTRL.iab[0] * CTRL.cosT_MT + CTRL.iab[1] * CTRL.sinT_MT
    CTRL.iMT[1] = CTRL.iab[0] *-CTRL.sinT_MT + CTRL.iab[1] * CTRL.cosT_MT

    # Speed regulation (SFOC)
    reg_speed.setpoint = CTRL.cmd_rpm / 60 * 2*np.pi * CTRL.npp # [elec.rad/s]
    reg_speed.measurement = CTRL.omega_r_elec # [elec.rad/s]
    CTRL.velocity_loop_counter += 1
    if CTRL.velocity_loop_counter >= CTRL.velocity_loop_ceiling:
        CTRL.velocity_loop_counter = 0
        # incremental_pi(reg_speed)
        tustin_pid(reg_speed)

    # dq-frame current commands
    if CTRL.bool_apply_speed_closed_loop_control == True:
        CTRL.cmd_iMT[1] = reg_speed.Out
        # CTRL.cmd_idq[0] = 0.0 # for user specifying
    else:
        CTRL.cmd_iMT[1] = 0.5


    if True:
        # constant stator flux linkage
        # CTRL.cmd_psi_Ms = 0.0
        pass
    else:
        CTRL.cmd_idq[0] = 1.5
        # CTRL.cmd_psi_Ms = np.sqrt( (CTRL.Ld*CTRL.cmd_i??[0])**2 + (CTRL.Lq * CTRL.cmd_i??[1])**2 )

    # CTRL.cmd_uMT[0] = CTRL.kPFL*( CTRL.cmd_psi_Ms - CTRL.psi_stator_MT_fb[0] ) + 0.8*CTRL.R*CTRL.iMT[0] + 0
    CTRL.cmd_uMT[0] = CTRL.kPFL*( CTRL.cmd_psi_Ms - CTRL.psi_stator_MT_fb[0] ) + 1.0*CTRL.R*CTRL.iMT[0] + 0

    # CTRL.omega_slip = CTRL.omega_syn - CTRL.omega_r_elec
    if True:
        CTRL.emf_stator[0] = CTRL.cmd_uab[0] - CTRL.R * CTRL.iab[0]
        CTRL.emf_stator[1] = CTRL.cmd_uab[1] - CTRL.R * CTRL.iab[1]
        # CTRL.omega_syn = - (CTRL.psi_stator_ab_fb[0]*-CTRL.emf_stator[1] + CTRL.psi_stator_ab_fb[1]*CTRL.emf_stator[0]) / (CTRL.psi_stator_MT_fb[0]*CTRL.psi_stator_MT_fb[0])
        CTRL.omega_slip = - (CTRL.psi_stator_ab_fb[0]*-CTRL.emf_stator[1] + CTRL.psi_stator_ab_fb[1]*CTRL.emf_stator[0]) / (CTRL.psi_stator_MT_fb[0]*CTRL.psi_stator_MT_fb[0])
    else:
        x = CTRL.psi_stator_ab_fb[0]
        y = CTRL.psi_stator_ab_fb[1]
        xdot = CTRL.cmd_uab[0] - CTRL.iab[0]*CTRL.R
        ydot = CTRL.cmd_uab[1] - CTRL.iab[1]*CTRL.R
        # CTRL.omega_syn = (ydot*x - xdot*y) / (x**2 + y**2)
        CTRL.omega_slip = (ydot*x - xdot*y) / (x**2 + y**2)

    CTRL.Intergral_of_iT_error += CTRL.CL_TS * CTRL.kPCL * (CTRL.cmd_iMT[1] - CTRL.iMT[1])
    CTRL.cmd_uMT[1] = CTRL.R * (CTRL.Intergral_of_iT_error - CTRL.cmd_iMT[1]) + CTRL.omega_syn*CTRL.psi_stator_MT_fb[0]

############################################# DSP SECTION
@njit(nogil=True)
def DSP(ACM, CTRL, reg_speed, reg_id, reg_iq):
    CTRL.timebase += CTRL.CL_TS

    """ Measurement """
    CTRL.iab[0] = ACM.iAlfa
    CTRL.iab[1] = ACM.iBeta
    CTRL.theta_d = ACM.theta_d

    # STATOR FLXU MEASUREMENT FOR SFOC
    CTRL.psi_stator_ab_fb[0] = ACM.KA*np.cos(ACM.theta_d) + CTRL.Lq*CTRL.iab[0]
    CTRL.psi_stator_ab_fb[1] = ACM.KA*np.sin(ACM.theta_d) + CTRL.Lq*CTRL.iab[1]
    # CTRL.theta_m = np.atan2(CTRL.psi_stator_ab_fb[1], CTRL.psi_stator_ab_fb[0])
    CTRL.psi_stator_MT_fb[0] = np.sqrt(CTRL.psi_stator_ab_fb[0]**2 + CTRL.psi_stator_ab_fb[1]**2)
    CTRL.psi_stator_MT_fb[1] = 0.0
    if CTRL.psi_stator_MT_fb[0] > 0:
        inverse = 1.0 / CTRL.psi_stator_MT_fb[0]
        CTRL.cosT_MT = CTRL.psi_stator_ab_fb[0] * inverse
        CTRL.sinT_MT = CTRL.psi_stator_ab_fb[1] * inverse
    else:
        CTRL.cosT_MT = 1.0
        CTRL.sinT_MT = 0.0
    CTRL.omega_syn = ACM.omega_syn
    CTRL.omega_slip = ACM.omega_syn

    """ Speed Estimation """
    if CTRL.index_separate_speed_estimation == 0:
        #TODO simulate the encoder
        CTRL.omega_r_elec = ACM.omega_r_elec
    elif CTRL.index_separate_speed_estimation == 1:
        RK4_ObserverSolver_CJH_Style(DYNAMICS_SpeedObserver, CTRL.xS, CTRL.CL_TS, CTRL)
        while CTRL.xS[0]> np.pi: CTRL.xS[0] -= 2*np.pi
        while CTRL.xS[0]<-np.pi: CTRL.xS[0] += 2*np.pi
        CTRL.iab_prev[0] = CTRL.iab_curr[0]
        CTRL.iab_prev[1] = CTRL.iab_curr[1]
        # CTRL.uab_prev[0] = CTRL.uab_curr[0] # This is needed only if voltage is measured, e.g., by eCAP. Remember to update the code below marked by [$].
        # CTRL.uab_prev[1] = CTRL.uab_curr[1] # This is needed only if voltage is measured, e.g., by eCAP. Remember to update the code below marked by [$].

        """ Speed Observer Outputs """
        CTRL.vartheta_d = CTRL.xS[0]
        CTRL.omega_r_elec = CTRL.xS[1]
        if CTRL.use_disturbance_feedforward_rejection == 0:
            CTRL.total_disrubance_feedforward = 0.0
        if CTRL.use_disturbance_feedforward_rejection == 1:
            CTRL.total_disrubance_feedforward = CTRL.xS[2]
        elif CTRL.use_disturbance_feedforward_rejection == 2:
            CTRL.total_disrubance_feedforward = CTRL.xS[2] + CTRL.ell2*CTRL.speed_observer_output_error

    """ (Optional) Do Park transformation again using the position estimate from the speed observer """

    """ Speed and Current Controller (two cascaded closed loops) """
    if CTRL.bool_use_FOC_or_SFOC == True:
        FOC(CTRL, reg_speed, reg_id, reg_iq)

        # [$] Inverse Park transformation: get voltage commands in alpha-beta frame as SVPWM input
        CTRL.cmd_uab[0] = CTRL.cmd_udq[0] * CTRL.cosT + CTRL.cmd_udq[1] *-CTRL.sinT
        CTRL.cmd_uab[1] = CTRL.cmd_udq[0] * CTRL.sinT + CTRL.cmd_udq[1] * CTRL.cosT

    else:
        SFOC_Dynamic(CTRL, reg_speed, reg_id, reg_iq)

        # [$] Inverse Park transformation: get voltage commands in alpha-beta frame as SVPWM input
        CTRL.cmd_uab[0] = CTRL.cmd_uMT[0] * CTRL.cosT_MT + CTRL.cmd_uMT[1] *-CTRL.sinT_MT
        CTRL.cmd_uab[1] = CTRL.cmd_uMT[0] * CTRL.sinT_MT + CTRL.cmd_uMT[1] * CTRL.cosT_MT


############################################# Inverter and PWM
@njit(nogil=True)
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

@njit(nogil=True)
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



############################################# Wrapper level 1 (Main simulation | Incremental Edition)
""" MAIN for Real-time simulation """
@njit(nogil=True)
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

@njit(nogil=True)
def ACMSimPyIncremental(t0, TIME, ACM=None, CTRL=None, reg_id=None, reg_iq=None, reg_speed=None):

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
    watch_data    = np.zeros( (100, len(machine_times)) ) # new
    # control_times = np.arange(t0, t0+TIME, CTRL.CL_TS)
    # watch_data = np.zeros( (40, len(control_times)) ) # old

    # Main loop
    jj = controller_down_sampling_ceiling # run controller at step 1
    watch_index = 0
    for ii in range(len(machine_times)):

        t = machine_times[ii]

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

            """ Console @ CL_TS """
            if CTRL.bool_overwrite_speed_commands == False:
                if t < 1.0:
                    CTRL.cmd_rpm = 50
                elif t < 1.5:
                    ACM.TLoad = 2
                elif t < 2.0:
                    CTRL.cmd_rpm = 200
                elif t < 3.0:
                    CTRL.cmd_rpm = -200
                elif t < 4.0:
                    CTRL.cmd_rpm = 0
                elif t < 4.5:
                    CTRL.cmd_rpm = 2000
                elif t < 5:
                    CTRL.cmd_idq[0] = 2
                elif t < 5.5:
                    ACM.TLoad = 0.0
                elif t < 6: 
                    CTRL.CMD_SPEED_SINE_RPM = 500
                # else: # don't implement else to receive commands from IPython console

                # if CTRL.CMD_SPEED_SINE_RPM!=0:
                #     CTRL.cmd_rpm = CTRL.CMD_SPEED_SINE_RPM * np.sin(2*np.pi*CTRL.CMD_SPEED_SINE_HZ*t)
                pass

            if CTRL.bool_apply_sweeping_frequency_excitation == True:

                if CTRL.timebase > CTRL.CMD_SPEED_SINE_END_TIME:
                    # next frequency
                    CTRL.CMD_SPEED_SINE_HZ += CTRL.CMD_SPEED_SINE_STEP_SIZE
                    # next end time
                    CTRL.CMD_SPEED_SINE_LAST_END_TIME = CTRL.CMD_SPEED_SINE_END_TIME
                    CTRL.CMD_SPEED_SINE_END_TIME += 1.0/CTRL.CMD_SPEED_SINE_HZ # 1.0 Duration for each frequency

                if CTRL.CMD_SPEED_SINE_HZ > CTRL.CMD_SPEED_SINE_HZ_CEILING:
                    # stop
                    CTRL.cmd_rpm = 0.0
                    CTRL.cmd_idq[1] = 0.0
                else:
                    # speed control - closed-loop sweep
                    CTRL.cmd_rpm    = CTRL.CMD_SPEED_SINE_RPM      * np.sin(2*np.pi*CTRL.CMD_SPEED_SINE_HZ*(CTRL.timebase - CTRL.CMD_SPEED_SINE_LAST_END_TIME))

                    # speed control - open-loop sweep
                    CTRL.cmd_idq[1] = CTRL.CMD_CURRENT_SINE_AMPERE * np.sin(2*np.pi*CTRL.CMD_SPEED_SINE_HZ*(CTRL.timebase - CTRL.CMD_SPEED_SINE_LAST_END_TIME))

            """ DSP @ CL_TS """
            # print(ii+1)
            DSP(ACM=ACM,
                CTRL=CTRL,
                reg_speed=reg_speed,
                reg_id=reg_id,
                reg_iq=reg_iq)

            # DEBUG
            # CTRL.cmd_uab[0] = 10*np.cos(5*2*np.pi*CTRL.timebase)
            # CTRL.cmd_uab[1] = 10*np.sin(5*2*np.pi*CTRL.timebase)

            # SVPWM for voltage source inverter
            svgen1.Ualfa = CTRL.cmd_uab[0]
            svgen1.Ubeta = CTRL.cmd_uab[1]
            SVGEN_DQ(svgen1, one_over_Vdc)
            # 高低有效逻辑翻转（仿真里得马上反回来，否则输出就反相了）
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

        """ Watch @ MACHINE_TS """
        watch_data[ 0][watch_index] = divmod(ACM.theta_d, 2*np.pi)[1]
        watch_data[ 1][watch_index] = ACM.omega_r_mech / (2*np.pi) * 60 # omega_r_mech
        watch_data[ 2][watch_index] = ACM.KA
        watch_data[ 3][watch_index] = ACM.iD
        watch_data[ 4][watch_index] = ACM.iQ
        watch_data[ 5][watch_index] = ACM.Tem
        watch_data[ 6][watch_index] =   CTRL.iab[0]
        watch_data[ 7][watch_index] =   CTRL.iab[1]
        watch_data[ 8][watch_index] = CTRL.idq[0]
        watch_data[ 9][watch_index] = CTRL.idq[1]
        watch_data[10][watch_index] = divmod(CTRL.theta_d, 2*np.pi)[1]
        watch_data[11][watch_index] = CTRL.omega_r_elec / (2*np.pi*ACM.npp) * 60
        watch_data[12][watch_index] = CTRL.cmd_rpm
        watch_data[13][watch_index] = CTRL.cmd_idq[0]
        watch_data[14][watch_index] = CTRL.cmd_idq[1]
        watch_data[15][watch_index] = CTRL.xS[0] # theta_d
        watch_data[16][watch_index] = CTRL.xS[1] / (2*np.pi*ACM.npp) * 60 # omega_r_elec
        watch_data[17][watch_index] = CTRL.xS[2] # TL
        watch_data[18][watch_index] = CTRL.xS[3] # pT
        watch_data[19][watch_index] = CTRL.KA
        watch_data[20][watch_index] = CTRL.KE
        watch_data[21][watch_index] = CTRL.xT[0] # stator flux[0]
        watch_data[22][watch_index] = CTRL.xT[1] # stator flux[1]
        watch_data[23][watch_index] = CTRL.xT[2] # I term
        watch_data[24][watch_index] = CTRL.xT[3] # I term
        watch_data[25][watch_index] = 0.0 # CTRL.active_flux[0] # active flux[0]
        watch_data[26][watch_index] = 0.0 # CTRL.active_flux[1] # active flux[1]
        watch_data[27][watch_index] = CTRL.Tem

        watch_data[28][watch_index] = CTRL.cmd_uab[0] # -svgen1.line_to_line_voltage_AC # ACM.uab[0] # CTRL.cmd_uab[0]
        watch_data[29][watch_index] = CTRL.cmd_uab[1] # svgen1.line_to_line_voltage_BC # svgen1.carrier_counter # CTRL.cmd_uab[1]

        watch_data[30][watch_index] = 30+svgen1.voltage_potential_at_terminal[0] # -svgen1.line_to_line_voltage_AC # ACM.uab[0]
        watch_data[31][watch_index] = svgen1.voltage_potential_at_terminal[1] # svgen1.line_to_line_voltage_BC # ACM.uab[1]
        watch_data[32][watch_index] = -30+svgen1.voltage_potential_at_terminal[2] # svgen1.line_to_line_voltage_AB # svgen1.line_to_line_voltage_BC
        watch_data[33][watch_index] = 0.0 # svgen1.voltage_potential_at_terminal[0] # svgen1.deadtime_counter[0] # svgen1.voltage_potential_at_terminal[0]
        watch_data[34][watch_index] = ACM.uab[0] # svgen1.deadtime_counter[1] # svgen1.voltage_potential_at_terminal[1]
        watch_data[35][watch_index] = ACM.uab[1] # svgen1.deadtime_counter[2] # svgen1.voltage_potential_at_terminal[2]

        watch_data[36][watch_index] = ACM.udq[0]
        watch_data[37][watch_index] = ACM.udq[1]
        watch_data[38][watch_index] = CTRL.cmd_udq[0]
        watch_data[39][watch_index] = CTRL.cmd_udq[1]

        watch_data[40][watch_index] = reg_speed.OutLimit

        watch_data[41][watch_index] = ACM.TLoad

        watch_data[42][watch_index] = CTRL.omega_syn
        watch_data[43][watch_index] = CTRL.cmd_uMT[0]
        watch_data[44][watch_index] = CTRL.cmd_uMT[1]

        watch_data[45][watch_index] = CTRL.cmd_iMT[1]
        watch_data[46][watch_index] = CTRL.iMT[1]

        watch_data[47][watch_index] = CTRL.cmd_psi_Ms
        watch_data[48][watch_index] = CTRL.psi_stator_MT_fb[0]
        watch_data[49][watch_index] = CTRL.omega_slip

        watch_data[50][watch_index] = CTRL.iMT[0]

        watch_index += 1

    # return machine_times, watch_data # old
    return machine_times, watch_data # new



############################################# Wrapper level 2 (Collect waveforms data based off user specified names)
# TODO: need to make this globally shared between the simulation and the GUI.
_Unit_Watch_Mapping = [
    '[rad]=ACM.theta_d',
    '[rad/s]=ACM.omega_r_mech',
    '[Wb]=ACM.KA',
    '[A]=ACM.iD',
    '[A]=ACM.iQ',
    '[Nm]=ACM.Tem',
    '[A]=CTRL.iab[0]',
    '[A]=CTRL.iab[1]',
    '[A]=CTRL.idq[0]',
    '[A]=CTRL.idq[1]',
    '[rad]=CTRL.theta_d',
    '[rpm]=CTRL.omega_r_mech',
    '[rpm]=CTRL.cmd_rpm',
    '[A]=CTRL.cmd_idq[0]',
    '[A]=CTRL.cmd_idq[1]',
    '[rad]=CTRL.xS[0]',  # theta_d
    '[rpm]=CTRL.xS[1]',  # omega_r_elec
    '[Nm]=CTRL.xS[2]',   # -TL
    '[Nm/s]=CTRL.xS[3]', # DL
    '[Wb]=CTRL.KA',
    '[Wb]=CTRL.KE',
    '[Wb]=CTRL.xT[0]', # stator flux[0]
    '[Wb]=CTRL.xT[1]', # stator flux[1]
    '[V]=CTRL.xT[2]', # I term
    '[V]=CTRL.xT[3]', # I term
    '[Wb]=CTRL.active_flux[0]', # active flux[0]
    '[Wb]=CTRL.active_flux[1]', # active flux[1]
    '[Nm]=CTRL.Tem',
    '[V]=CTRL.cmd_uab[0]',
    '[V]=CTRL.cmd_uab[1]',
    '[1]=svgen1.S1',
    '[1]=svgen1.S2',
    '[1]=svgen1.S3',
    '[1]=svgen1.S4',
    '[1]=svgen1.S5',
    '[1]=svgen1.S6',
    '[V]=ACM.udq[0]',
    '[V]=ACM.udq[1]',
    '[V]=CTRL.cmd_udq[0]',
    '[V]=CTRL.cmd_udq[1]',
    '[A]=reg_speed.OutLimit',
    '[Nm]=ACM.TLoad',
    '[rad/s]=CTRL.omega_syn',
    '[V]=CTRL.cmd_uMT[0]',
    '[V]=CTRL.cmd_uMT[1]',
    '[A]=CTRL.cmd_iMT[1]',
    '[A]=CTRL.iMT[1]',
    '[Wb]=CTRL.cmd_psi_Ms',
    '[Wb]=CTRL.psi_stator_MT_fb[0]',
    '[rad/s]=CTRL.omega_slip',
    '[A]=CTRL.iMT[0]',
]
Watch_Mapping = [el[el.find('=')+1:] for el in _Unit_Watch_Mapping] # remove units before "="

def ACMSimPyWrapper(numba__scope_dict, *arg, **kwarg):

    # Do Numerical Integrations (that do not care about numba__scope_dict at all and return watch_data whatsoever)
    machine_times, watch_data = ACMSimPyIncremental(*arg, **kwarg)
    # print(f'{len(watch_data[0])=}。 end_time', machine_times[-1])
    watch_data_as_dict = dict(zip(Watch_Mapping, watch_data))
    # print(watch_data_as_dict.keys())

    # Post-processing
    numba__waveforms_dict = dict()
    if True:
        # option 1 (with exec)
        for key, expressions in numba__scope_dict.items():
            # key = r'Error Speed [rpm]',
            # expressions = ('CTRL.cmd_rpm-ACM.omega_r_mech', 'CTRL.idq[1]'),
            waveforms = []
            for expression in expressions:
                # expression = 'CTRL.cmd_rpm-ACM.omega_r_mech'
                translated_expression = ''
                for word in expression.split():
                    if 'CTRL' in word or 'ACM' in word or 'reg_' in word or 'svgen1' in word: # 这里逻辑有点怪，每增加一个新的结构体（比如svgen1）都要在修改这一行？
                        translated_expression += f'watch_data_as_dict["{word}"]'
                    else:
                        translated_expression += word
                # print('DEBUG', translated_expression)
                waveforms.append(eval(translated_expression))
            numba__waveforms_dict[key] = waveforms
        # for key, val in numba__waveforms_dict.items():
        #     print(key, len(val), len(val[0]))
        # quit()
    else:
        # option 2 (without using exec)
        for key, values in numba__scope_dict.items():
            # key = '$\alpha\beta$ current [A]'
            # values = ('CTRL.iab', 'CTRL.idq[1]'),
            waveforms = []
            for val in values:
                # val = 'CTRL.iab'
                for index, mapping in enumerate(Watch_Mapping):
                    # 'CTRL.iab' in '[A]=CTRL.iab[0]'
                    # 'CTRL.iab' in '[A]=CTRL.iab[1]'
                    if val in mapping:
                        waveforms.append(watch_data[index])
                        # print('\t', key, val, 'in', mapping)
                        if len(val) == 1:
                            raise Exception('Invalid numba__scope_dict, make sure it is a dict of tuples of strings.')

            numba__waveforms_dict[key] = waveforms
    # quit()
    return machine_times, numba__waveforms_dict



############################################# Wrapper level 3 (User Interface)
from collections import OrderedDict as OD
class Simulation_Benchmark:
    def __init__(self, d, tuner=None, bool_start_simulation=True, CTRL_execute_codes=''):

        self.CTRL_execute_codes = CTRL_execute_codes
        self.d = d
        print('Simulation_Benchmark')
        # Auto-tuning PI
        if d['CL_SERIES_KP'] is None:
            if tuner is None:
                # sys.path.append(os.path.join(os.path.dirname(__file__), "tuner"))
                import tuner
            tuner.tunner_wrapper(d)
            print('\tAuto tuning...')
            print(f'\t{d=}\n')
        else:
            print('\tSkip tuning.')

        # 这个字典决定你要画的波形是哪些信号，具体能用的信号见：_Watch_Mapping
        # 允许信号之间进行简单的加减乘除，比如：'CTRL.cmd_rpm - CTRL.omega_r_mech'
        numba__scope_dict = OD([
            # Y Labels                        Signal Name of Traces
            (r'Speed Out Limit [A]',          ( 'reg_speed.OutLimit'                                ,) ),
            (r'$q$-axis voltage [V]',         ( 'ACM.udq[1]', 'CTRL.cmd_udq[1]'                     ,) ),
            (r'$d$-axis voltage [V]',         ( 'ACM.udq[0]', 'CTRL.cmd_udq[0]'                     ,) ),
            (r'Torque [Nm]',                  ( 'ACM.Tem', 'CTRL.Tem'                               ,) ),
            (r'Speed [rpm]',                  ( 'CTRL.cmd_rpm', 'CTRL.omega_r_mech', 'CTRL.xS[1]'   ,) ),
            (r'Speed Error [rpm]',            ( 'CTRL.cmd_rpm - CTRL.omega_r_mech'                  ,) ),
            (r'Position [rad]',               ( 'ACM.theta_d', 'CTRL.theta_d', 'CTRL.xS[0]'         ,) ),
            (r'Position mech [rad]',          ( 'ACM.theta_d'                                       ,) ),
            (r'$q$-axis current [A]',         ( 'ACM.iQ', 'CTRL.cmd_idq[1]'                         ,) ),
            (r'$d$-axis current [A]',         ( 'ACM.iD', 'CTRL.cmd_idq[0]'                         ,) ),
            (r'K_{\rm Active} [A]',           ( 'ACM.KA', 'CTRL.KA'                                 ,) ),
            (r'Load torque [Nm]',             ( 'ACM.TLoad', 'CTRL.xS[2]'                           ,) ),
            (r'CTRL.iD [A]',                  ( 'CTRL.cmd_idq[0]', 'CTRL.idq[0]'                    ,) ),
            (r'CTRL.iQ [A]',                  ( 'CTRL.cmd_idq[1]', 'CTRL.idq[1]'                    ,) ),
            (r'CTRL.uab [V]',                 ( 'CTRL.cmd_uab[0]', 'CTRL.cmd_uab[1]'                ,) ),
            (r'S [1]',                        ( 'svgen1.S1', 'svgen1.S2', 'svgen1.S3', 'svgen1.S4', 'svgen1.S5', 'svgen1.S6' ,) ),
            (r'$M$-axis voltage [V]',         ( 'CTRL.cmd_uMT[0]'                                   ,) ),
            (r'$T$-axis voltage [V]',         ( 'CTRL.cmd_uMT[1]'                                   ,) ),
            (r'omega_syn [rad/s]',            ( 'CTRL.omega_syn', 'CTRL.omega_slip'                 ,) ),
            (r'Current [A]',                  ( 'CTRL.cmd_iMT[1]', 'CTRL.iMT[1]', 'CTRL.iMT[0]'     ,) ),
            (r'M-axis Flux [Wb]',             ( 'CTRL.cmd_psi_Ms', 'CTRL.psi_stator_MT_fb[0]'       ,) ),
        ])

        if bool_start_simulation:
            self.start_simulation_slices(d, numba__scope_dict)

    def get_global_objects(self):
        d = self.d
        # init
        CTRL = The_Motor_Controller(CL_TS = d['CL_TS'],
                                    VL_TS = d['VL_EXE_PER_CL_EXE']*d['CL_TS'],
                                    init_npp = d['init_npp'],
                                    init_IN = d['init_IN'],
                                    init_R = d['init_R'],
                                    init_Ld = d['init_Ld'],
                                    init_Lq = d['init_Lq'],
                                    init_KE = d['init_KE'],
                                    init_KA = d['init_KA'],
                                    init_Rreq = d['init_Rreq'],
                                    init_Js = d['init_Js'],
                                    DC_BUS_VOLTAGE = d['DC_BUS_VOLTAGE'])
        CTRL.bool_apply_decoupling_voltages_to_current_regulation = d['CTRL.bool_apply_decoupling_voltages_to_current_regulation']
        CTRL.bool_apply_sweeping_frequency_excitation = d['CTRL.bool_apply_sweeping_frequency_excitation']
        CTRL.bool_apply_speed_closed_loop_control = d['CTRL.bool_apply_speed_closed_loop_control']
        CTRL.bool_overwrite_speed_commands = d['CTRL.bool_overwrite_speed_commands']
        CTRL.bool_zero_id_control = d['CTRL.bool_zero_id_control']
        CTRL.bool_use_FOC_or_SFOC = d['CTRL.bool_use_FOC_or_SFOC']
        CTRL.omega_syn = 50

        ACM       = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=d['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])

        reg_dispX = The_PID_Regulator(d['disp.Kp'], d['disp.Ki'], d['disp.Kd'], d['disp.tau'], d['disp.OutLimit'], d['disp.IntLimit'], d['CL_TS'])
        reg_dispY = The_PID_Regulator(d['disp.Kp'], d['disp.Ki'], d['disp.Kd'], d['disp.tau'], d['disp.OutLimit'], d['disp.IntLimit'], d['CL_TS'])

        if False:
            # Use incremental_pi codes
            reg_id    = The_PI_Regulator(d['CL_SERIES_KP'], d['CL_SERIES_KP']*d['CL_SERIES_KI']*CTRL.CL_TS, d['DC_BUS_VOLTAGE']/1.732) # 我们假设调制方式是SVPWM，所以母线电压就是输出电压的线电压最大值，而我们用的是恒相幅值变换，所以限幅是相电压。
            reg_iq    = The_PI_Regulator(d['CL_SERIES_KP'], d['CL_SERIES_KP']*d['CL_SERIES_KI']*CTRL.CL_TS, d['DC_BUS_VOLTAGE']/1.732) # 我们假设调制方式是SVPWM，所以母线电压就是输出电压的线电压最大值，而我们用的是恒相幅值变换，所以限幅是相电压。
            reg_speed = The_PI_Regulator(d['VL_SERIES_KP'], d['VL_SERIES_KP']*d['VL_SERIES_KI']*CTRL.VL_TS, d['VL_LIMIT_OVERLOAD_FACTOR']*1.414*d['init_IN']) # IN 是线电流有效值，我们这边限幅是用的电流幅值。
        else:
            # Use tustin_pi codes
            local_Kp = d['CL_SERIES_KP']
            if d['CTRL.bool_apply_decoupling_voltages_to_current_regulation'] == False:
                local_Ki = d['CL_SERIES_KP']*d['CL_SERIES_KI'] * d['FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False']
                print('\tNote bool_apply_decoupling_voltages_to_current_regulation is False, to improve the current regulator performance a factor of %g has been multiplied to CL KI.' % (d['FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False']))
            else:
                local_Ki = d['CL_SERIES_KP']*d['CL_SERIES_KI']
            local_Kd = 0.0
            local_tau = 0.0
            local_OutLimit =     d['DC_BUS_VOLTAGE']/1.732
            local_IntLimit = 1.0*d['DC_BUS_VOLTAGE']/1.732 # Integrator having a lower output limit makes no sense. For example, the q-axis current regulator needs to cancel back emf using the integrator output for almost full dc bus voltage at maximum speed.
            reg_id    = The_PID_Regulator(local_Kp, local_Ki, local_Kd, local_tau, local_OutLimit, local_IntLimit, d['CL_TS'])
            reg_iq    = The_PID_Regulator(local_Kp, local_Ki, local_Kd, local_tau, local_OutLimit, local_IntLimit, d['CL_TS'])
            print('\t', reg_id.OutLimit, 'V')

            local_Kp = d['VL_SERIES_KP']
            local_Ki = d['VL_SERIES_KP']*d['VL_SERIES_KI']
            local_Kd = 0.0
            local_tau = 0.0
            local_OutLimit =     d['VL_LIMIT_OVERLOAD_FACTOR']*1.414*d['init_IN']
            local_IntLimit = 1.0*d['VL_LIMIT_OVERLOAD_FACTOR']*1.414*d['init_IN']
            reg_speed = The_PID_Regulator(local_Kp, local_Ki, local_Kd, local_tau, local_OutLimit, local_IntLimit, CTRL.VL_TS)
            print('\t', reg_speed.OutLimit, 'A')

        exec(self.CTRL_execute_codes)
        print(CTRL.cmd_psi_Ms)

        return CTRL, ACM, reg_id, reg_iq, reg_speed, reg_dispX, reg_dispY

    def start_simulation_slices(self, d, numba__scope_dict):

        global_objects = self.get_global_objects()
        self.CTRL, self.ACM, self.reg_id, self.reg_iq, self.reg_speed, self.reg_dispX, self.reg_dispY = CTRL, ACM, reg_id, reg_iq, reg_speed, reg_dispX, reg_dispY = global_objects

        global_trace_names = []
        max_number_of_traces = 0
        for ylabel, trace_names in numba__scope_dict.items():
            for name in trace_names:
                max_number_of_traces += 1
            for trace_index, name in enumerate(trace_names):
                global_trace_names.append(name)
        print(f'\t{max_number_of_traces=}')

        # init global data arrays for plotting
        global_arrays = [None] * max_number_of_traces
        global_machine_times = None

        def save_to_global(_global, _local):
            return _local if _global is None else np.append(_global, _local)

        # simulate to generate NUMBER_OF_SLICES*TIME_SLICE sec of data
        for ii in range(d['NUMBER_OF_SLICES']):

            exec(d['user_system_input_code']) # 和 CONSOLE.user_controller_commands 功能相同
            # if ii < 5:
            #     CTRL.cmd_rpm = 50
            # else:
            #     ACM.TLoad = 5

            # perform animation step
            machine_times, numba__waveforms_dict = \
                ACMSimPyWrapper(numba__scope_dict,
                            t0=ii*d['TIME_SLICE'], TIME=d['TIME_SLICE'], 
                            ACM=ACM,
                            CTRL=CTRL,
                            reg_id=reg_id,
                            reg_iq=reg_iq,
                            reg_speed=reg_speed)

            # and save slice data to global data variables
            global_machine_times = save_to_global(global_machine_times, machine_times)
            global_index = 0
            for ylabel in numba__scope_dict.keys():
                for trace_index, local_trace_data in enumerate(numba__waveforms_dict[ylabel]):
                    # trace data
                    global_arrays[global_index] = save_to_global(global_arrays[global_index], local_trace_data)

                    # next
                    global_index += 1

        # map global data to global names
        gdd = global_data_dict = OD()
        for name, array in zip(global_trace_names, global_arrays):
            global_data_dict[name] = array

        print(f'\t{gdd.keys()=}')

        self.global_machine_times = global_machine_times
        self.gdd = gdd

def lpf1_inverter(array):
    y_tminus1 = 0.0
    new_array = []
    for x in array:
        new_x = y_tminus1 + 5* 0.00020828993959591752 * (x - y_tminus1)
        y_tminus1 = new_x
        new_array.append(y_tminus1)
    return new_array

if __name__ == '__main__':
    # User input:
    d = d_user_input_motor_dict = {
        # Timing
        'CL_TS': 1e-4,
        'VL_EXE_PER_CL_EXE': 1,
        'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1, # 500,
        'TIME_SLICE': 0.2,
        'NUMBER_OF_SLICES': 6,
        # Motor data
        # 'init_npp': 2,
        # 'init_IN': 4.6,
        # 'init_R': 5.5,
        # 'init_Ld': 602*1e-3,
        # 'init_Lq':  22*1e-3,
        # 'init_KE': 0,
        # 'init_KA': 580*1e-3*1.5, # Lmu * ids_rated
        # 'init_Rreq': 2.1,
        'init_Js': 0.044,
        # Motor data
        'init_npp': 22,
        'init_IN': 1.3*6/1.414,
        'init_R': 0.035,
        'init_Ld': 1*0.036*1e-3,
        'init_Lq': 1*0.036*1e-3,
        'init_KE': 0.0125,
        'init_KA': 0.0125,
        'init_Rreq': 0.0,
        # 'init_Js': 0.44*1e-4,
        ##########################
        'DC_BUS_VOLTAGE': 48,
        'user_system_input_code': '''if ii < 1: CTRL.cmd_idq[0] = 0.0; CTRL.cmd_rpm = 50 \nelif ii <5: ACM.TLoad = 5 \nelif ii <100: CTRL.cmd_rpm = -50''',
        # Controller config
        'CTRL.bool_apply_speed_closed_loop_control': True,
        'CTRL.bool_apply_decoupling_voltages_to_current_regulation': False,
        'CTRL.bool_apply_sweeping_frequency_excitation': False,
        'CTRL.bool_overwrite_speed_commands': True,
        'CTRL.bool_zero_id_control': True,
        'CTRL.bool_use_FOC_or_SFOC': False,
        'FOC_delta': 10, # 25, # 6.5
        'FOC_desired_VLBW_HZ': 20, # 60
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
    print(f'最大电流上升率 {d["DC_BUS_VOLTAGE"]/1.732/d["init_Lq"]*1e-3} A/ms、最大转速上升率 rpm/s')

    图 = 1 # 空载加速、加载、反转
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

        # ax = axes[1]
        # ax.plot(global_machine_times, gdd['CTRL.cmd_idq[0]'], label=r'$i_d^*$')
        # ax.plot(global_machine_times, gdd['CTRL.idq[0]'], label=r'$i_d$')
        # # ax.plot(global_machine_times, gdd['ACM.iD'])
        # ax.set_ylabel(r'$i_d$ [A]', multialignment='center') #, fontdict=font)

        # ax = axes[2]
        # ax.plot(global_machine_times, gdd['CTRL.cmd_idq[1]'], label=r'$i_q^*$')
        # ax.plot(global_machine_times, gdd['CTRL.idq[1]'], label=r'$i_q$')
        # ax.set_ylabel(r'$i_q$ [A]', multialignment='center') #, fontdict=font)

        ax = axes[1]
        ax.plot(global_machine_times, gdd['CTRL.cmd_uMT[0]'], label=r'$u_{M}$')
        ax.plot(global_machine_times, gdd['CTRL.cmd_uMT[1]'], label=r'$u_{T}$')
        ax.plot(global_machine_times, (gdd['ACM.udq[1]']), label=r'$u_q$')
        ax.set_ylabel(r'$u$ [V]', multialignment='center') #, fontdict=font)

        ax = axes[2]
        ax.plot(global_machine_times, gdd['CTRL.omega_syn'], label=r'$\omega_{\rm syn}$')
        ax.plot(global_machine_times, gdd['CTRL.omega_slip'], label=r'$\hat\omega_{\rm syn}$')
        ax.set_ylabel(r'$i_q$ [A]', multialignment='center') #, fontdict=font)

        ax = axes[3]
        ax.plot(global_machine_times, gdd['ACM.Tem'], label=r'ACM.$T_{\rm em}$')
        # ax.plot(global_machine_times, gdd['CTRL.Tem'], label=r'CTRL.$T_{\rm em}$')
        ax.set_ylabel(r'$T_{\rm em}$ [Nm]', multialignment='center') #, fontdict=font)

        ax = axes[4]
        ax.plot(global_machine_times, (gdd['CTRL.cmd_iMT[1]']), label=r'$iT*$')
        ax.plot(global_machine_times, (gdd['CTRL.iMT[1]']), label=r'$iT$')
        ax.plot(global_machine_times, (gdd['CTRL.iMT[0]']), label=r'$iM$')
        ax.set_ylabel(r'Current [A]', multialignment='center')

        # ax = axes[4]
        # ax.plot(global_machine_times, (gdd['ACM.udq[0]']), label=r'$u_d$') # lpf1_inverter
        # ax.plot(global_machine_times, gdd['CTRL.cmd_udq[0]'], label=r'$u_d^*$')
        # ax.set_ylabel(r'$u_d$ [V]', multialignment='center') #, fontdict=font)

        ax = axes[5]
        # ax.plot(global_machine_times, (gdd['CTRL.cmd_psi_Ms']), label=r'$\psi_{Ms}$^*')
        ax.plot(global_machine_times, gdd['CTRL.psi_stator_MT_fb[0]'], label=r'$\psi_{Ms}$')
        ax.set_ylabel(r'$\psi_{Ms}$ [Wb]', multialignment='center') #, fontdict=font)

        ax = axes[6]
        ax.plot(global_machine_times, gdd['CTRL.cmd_uab[0]'], label=r'$u_\alpha$')
        ax.plot(global_machine_times, gdd['CTRL.cmd_uab[1]'], label=r'$u_\beta$')
        ax.set_ylabel(r'$u_{\alpha,\beta}$ [V]', multialignment='center') #, fontdict=font)

        for ax in axes:
            ax.grid(True)
            ax.legend(loc=1)
            # for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
            #     tick.label.set_font(font)
        axes[-1].set_xlabel('Time [s]') #, fontdict=font)
        return fig


    temp_a = d['TIME_SLICE']
    temp_b = d['NUMBER_OF_SLICES']
    d['TIME_SLICE'] = 0.01
    d['NUMBER_OF_SLICES'] = 1
    sim1 = Simulation_Benchmark(d); gdd, global_machine_times = sim1.gdd, sim1.global_machine_times; #fig = 图1画图代码(); # fig.savefig(f'SliceFSPM-fig-{图}.pdf', dpi=400, bbox_inches='tight', pad_inches=0)
    d['TIME_SLICE'] = temp_a
    d['NUMBER_OF_SLICES'] = temp_b

    # plt.show()

#%%

CTRL_execute_codes = '''
CTRL.kPFL = 20
CTRL.kPCL = 600
reg_speed.Kp = 0.7
reg_speed.Ki = 20
CTRL.cmd_psi_Ms = 0.9
'''
CTRL_execute_codes = '''
CTRL.kPFL = 0
CTRL.kPCL = 300
reg_speed.Kp = 0.5
reg_speed.Ki = 20
CTRL.cmd_psi_Ms = 0.0125
'''

sim1 = Simulation_Benchmark(d, CTRL_execute_codes=CTRL_execute_codes); gdd, global_machine_times = sim1.gdd, sim1.global_machine_times; fig = 图1画图代码(); # fig.savefig(f'SliceFSPM-fig-{图}.pdf', dpi=400, bbox_inches='tight', pad_inches=0)
plt.show()

# %%
