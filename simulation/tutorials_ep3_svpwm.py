# %%
import enum
from numba.experimental import jitclass
from numba import njit, int32, float64
from pylab import np, plt
import re
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
            ('omega_r_elec', float64),
            ('omega_syn', float64),
            ('uab', float64[:]),
            # ('uab_prev', float64[:]),
            # ('uab_curr', float64[:]),
            ('iab', float64[:]),
            ('iab_prev', float64[:]),
            ('iab_curr', float64[:]),
            # states
            ('timebase', float64),
            ('KA', float64),
            ('Tem', float64),
            ('cosT', float64),
            ('sinT', float64),
            # commands
            ('cmd_idq', float64[:]),
            ('cmd_udq', float64[:]),
            ('cmd_uab', float64[:]),
            ('cmd_rpm', float64),
            ('index_separate_speed_estimation', int32),
            ('use_disturbance_feedforward_rejection', int32),
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
            # electrical parameters
            ('R',   float64),
            ('Ld',  float64),
            ('Lq',  float64),
            ('KE',  float64),
            ('Rreq',float64),
            # mechanical parameters
            ('Js',  float64),
        # OBSERVER
            # feedback / inputs
            ('idq', float64[:]),
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
        init_Rreq = -1, # note division by 0 is equal to infinity
        init_Js = 0.0006168,
    ):
        ''' CONTROL '''
        # constants
        self.CL_TS = CL_TS
        self.VL_TS = VL_TS
        self.velocity_loop_ceiling = VL_TS / CL_TS
        self.velocity_loop_counter = self.velocity_loop_ceiling - 1
        print('CTRL.velocity_loop_ceiling =', self.velocity_loop_ceiling)
        # feedback / input
        self.theta_d = 0.0
        self.omega_r_elec = 0.0
        self.omega_syn = 0.0
        self.uab      = np.zeros(2, dtype=np.float64)
        # self.uab_prev = np.zeros(2, dtype=np.float64)
        # self.uab_curr = np.zeros(2, dtype=np.float64)
        self.iab      = np.zeros(2, dtype=np.float64)
        self.iab_prev = np.zeros(2, dtype=np.float64)
        self.iab_curr = np.zeros(2, dtype=np.float64)
        # states
        self.timebase = 0.0
        self.KA = init_KE
        self.Tem = 0.0
        self.cosT = 1.0
        self.sinT = 0.0
        # commands
        self.cmd_idq = np.zeros(2, dtype=np.float64)
        self.cmd_udq = np.zeros(2, dtype=np.float64)
        self.cmd_uab = np.zeros(2, dtype=np.float64)
        self.cmd_rpm = 0.0
        self.index_separate_speed_estimation = 0
        self.use_disturbance_feedforward_rejection = 0
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
        self.Rreq = init_Rreq
        self.Js   = init_Js

        ''' OBSERVER '''
        # feedback / input
        self.idq = np.zeros(2, dtype=np.float64)
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
        ('KA', float64),
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
        self.KE  = CTRL.KE
        self.Rreq  = CTRL.Rreq
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

@jitclass(
    spec=[
        ('Kp', float64),
        ('Ki', float64),
        ('Err', float64),
        ('Ref', float64),
        ('Fbk', float64),
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
        self.Ref      = 0.0
        self.Fbk      = 0.0
        self.Out      = 0.0
        self.OutLimit = OUTPUT_LIMIT
        self.ErrPrev  = 0.0
        self.OutPrev  = 0.0

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

############################################# MACHINE SIMULATION SECTION

@njit(nogil=True)
def DYNAMICS_MACHINE(t, x, ACM, CLARKE_TRANS_TORQUE_GAIN=1.5):
    fx = np.zeros(ACM.NS)

    # omega_d_mech = x[0]
    # omega_r_mech = x[1]
    KA    = x[2]
    iD    = x[3]
    iQ    = x[4]
    # ACM.theta_d = x[0]*ACM.npp
    # ACM.omega_r = x[1]*ACM.npp
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
    reg.Err = reg.Ref - reg.Fbk
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
def FOC(CTRL, reg_speed, reg_id, reg_iq):
    reg_speed.Ref = CTRL.cmd_rpm / 60 * 2*np.pi * CTRL.npp # [elec.rad]
    reg_speed.Fbk = CTRL.omega_r_elec # [elec.rad]
    CTRL.velocity_loop_counter += 1
    if CTRL.velocity_loop_counter >= CTRL.velocity_loop_ceiling:
        CTRL.velocity_loop_counter = 0
        incremental_pi(reg_speed)

    # dq-frame current commands
    CTRL.cmd_idq[1] = reg_speed.Out
    CTRL.cmd_idq[0] = 0.0
    if CTRL.Rreq>0: # IM
        CTRL.cmd_idq[0] = 0.9 / (CTRL.Ld - CTRL.Lq) # [Wb] / [H]

    # d-axis
    reg_id.Ref = CTRL.cmd_idq[0]
    reg_id.Fbk = CTRL.idq[0]
    incremental_pi(reg_id)
    CTRL.cmd_udq[0] = reg_id.Out

        # if HUMAN.use_disturbance_feedforward_rejection == 0:
        #     CTRL.cmd_idq[1] = reg_speed.Out
        # else:
        #     CTRL.cmd_idq[1] = HUMAN.KP*(reg_speed.Ref-reg_speed.Fbk) + OB.total_disrubance_feedforward

    # q-axis
    reg_iq.Ref = CTRL.cmd_idq[1]
    reg_iq.Fbk = CTRL.idq[1]
    incremental_pi(reg_iq)
    CTRL.cmd_udq[1] = reg_iq.Out

    # return CTRL.cmd_udq


############################################# DSP SECTION

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

# print(180/np.pi*angle_diff(7, -7))
# print(180/np.pi*angle_diff(7, -6))
# print(180/np.pi*angle_diff(725/180*np.pi, 190/180*np.pi))
# print(180/np.pi*angle_diff(725/180*np.pi, 175/180*np.pi))
# print(180/np.pi*angle_diff((720-25)/180*np.pi, 190/180*np.pi))
# print(180/np.pi*angle_diff((720-25)/180*np.pi, 175/180*np.pi))
# quit()

""" DSP """
@njit(nogil=True)
def DSP(ACM, CTRL, reg_speed, reg_id, reg_iq):
    CTRL.timebase += CTRL.CL_TS

    """ Measurement """
    CTRL.iab[0] = ACM.iAlfa
    CTRL.iab[1] = ACM.iBeta
    CTRL.theta_d = ACM.theta_d

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
    FOC(CTRL, reg_speed, reg_id, reg_iq)

    # [$] Inverse Park transformation: get voltage commands in alpha-beta frame as SVPWM input
    CTRL.cmd_uab[0] = CTRL.cmd_udq[0] * CTRL.cosT + CTRL.cmd_udq[1] *-CTRL.sinT
    CTRL.cmd_uab[1] = CTRL.cmd_udq[0] * CTRL.sinT + CTRL.cmd_udq[1] * CTRL.cosT

############################################# Inverter and PWM
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
                v.bool_RisingEdgeDelay_is_active[0] = True
                v.S1 = False
            else:
                pass # False
        if v.carrier_counter >= v.EPwm2Regs_CMPA_bit_CMPA:
            v.deadtime_counter[1] += 1
            if v.deadtime_counter[1] <= DEAD_TIME_AS_COUNT:
                v.bool_RisingEdgeDelay_is_active[1] = True
                v.S2 = False
            else:
                pass # False
        if v.carrier_counter >= v.EPwm3Regs_CMPA_bit_CMPA:
            v.deadtime_counter[2] += 1
            if v.deadtime_counter[2] <= DEAD_TIME_AS_COUNT:
                v.bool_RisingEdgeDelay_is_active[2] = True
                v.S3 = False
            else:
                pass # False
    elif v.bool_counting_down == True:

        if v.carrier_counter < v.EPwm1Regs_CMPA_bit_CMPA:
            v.deadtime_counter[0] += 1
            if v.deadtime_counter[0] < DEAD_TIME_AS_COUNT:
                v.bool_FallingEdgeDelay_is_active[0] = True
                v.S4 = False
            else:
                pass # False
        if v.carrier_counter < v.EPwm2Regs_CMPA_bit_CMPA:
            v.deadtime_counter[1] += 1
            if v.deadtime_counter[1] < DEAD_TIME_AS_COUNT:
                v.bool_FallingEdgeDelay_is_active[1] = True
                v.S5 = False
            else:
                pass # False
        if v.carrier_counter < v.EPwm3Regs_CMPA_bit_CMPA:
            v.deadtime_counter[2] += 1
            if v.deadtime_counter[2] < DEAD_TIME_AS_COUNT:
                v.bool_FallingEdgeDelay_is_active[2] = True
                v.S6 = False
            else:
                pass # False
    # if v.bool_RisingEdgeDelay_is_active[0]:   v.S1 = False # AHC: upper leg implements RED
    # if v.bool_RisingEdgeDelay_is_active[1]:   v.S2 = False # AHC: upper leg implements RED
    # if v.bool_RisingEdgeDelay_is_active[2]:   v.S3 = False # AHC: upper leg implements RED
    # if v.bool_FallingEdgeDelay_is_active[0]:  v.S4 = False # AHC: lower leg implements FED
    # if v.bool_FallingEdgeDelay_is_active[1]:  v.S5 = False # AHC: lower leg implements FED
    # if v.bool_FallingEdgeDelay_is_active[2]:  v.S6 = False # AHC: lower leg implements FED

############################################# Wrapper level 1
""" MAIN for Real-time simulation """
@njit(nogil=True)
def ACMSimPyIncremental(
        t0, TIME,
        ACM=None,
        CTRL=None,
        reg_id=None,
        reg_iq=None,
        reg_speed=None,
    ):

    # RK4 simulation and controller execution relative freuqencies
    MACHINE_TS = CTRL.CL_TS / ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    controller_down_sampling_ceiling = int(CTRL.CL_TS / MACHINE_TS)

    # SVPWM
    CPU_TICK_PER_SAMPLING_PERIOD = ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    DEAD_TIME_AS_COUNT = int(200*0.5e-4*CPU_TICK_PER_SAMPLING_PERIOD) # 200 count for 0--5000--0 counting sequence
    print('DEAD_TIME_AS_COUNT =', DEAD_TIME_AS_COUNT)
    Vdc = 150 # Vdc is assumed measured and known
    one_over_Vdc = 1/Vdc
    # print('Vdc, CPU_TICK_PER_SAMPLING_PERIOD, controller_down_sampling_ceiling', Vdc, CPU_TICK_PER_SAMPLING_PERIOD, controller_down_sampling_ceiling)
    svgen1 = SVgen_Object(CPU_TICK_PER_SAMPLING_PERIOD)

    # watch variabels
    machine_times = np.arange(t0, t0+TIME, MACHINE_TS)
    # control_times = np.arange(t0, t0+TIME, CTRL.CL_TS)
    watch_data    = np.zeros( (40, len(machine_times)) ) # new
    # watch_data = np.zeros( (40, len(control_times)) ) # old

    # Main loop
    # print('\tt0 =', t0)
    jj = controller_down_sampling_ceiling # run controller at step 1
    watch_index = 0
    for ii in range(len(machine_times)):

        t = machine_times[ii]

        """ Machine Simulation @ MACHINE_TS """
        # Numerical Integration (ode4) with 5 states
        RK4_MACHINE(t, ACM, hs=MACHINE_TS)

        """ Machine Simulation Output @ MACHINE_TS """
        # Generate output variables for easy access
        # ACM.x[0] = ACM.x[0] - ACM.x[0]//(2*np.pi)*(2*np.pi)
        ACM.theta_d_mech = ACM.x[0]
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

                if (CTRL.CMD_SPEED_SINE_HZ > CTRL.CMD_SPEED_SINE_HZ_CEILING):
                    # stop
                    CTRL.cmd_rpm = 0.0
                    CTRL.cmd_idq[1] = 0.0
                else:
                    # speed control - closed-cloop sweep
                    CTRL.cmd_rpm = CTRL.CMD_SPEED_SINE_RPM * np.sin(2*np.pi*CTRL.CMD_SPEED_SINE_HZ*(CTRL.timebase - CTRL.CMD_SPEED_SINE_LAST_END_TIME))

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

            # """ Watch @ CL_TS """
            # watch_data[ 0][watch_index] = divmod(ACM.theta_d, 2*np.pi)[1]
            # watch_data[ 1][watch_index] = ACM.omega_r_mech / (2*np.pi) * 60 # omega_r_mech
            # watch_data[ 2][watch_index] = ACM.KA
            # watch_data[ 3][watch_index] = ACM.iD
            # watch_data[ 4][watch_index] = ACM.iQ
            # watch_data[ 5][watch_index] = ACM.Tem
            # watch_data[ 6][watch_index] =   CTRL.iab[0]
            # watch_data[ 7][watch_index] =   CTRL.iab[1]
            # watch_data[ 8][watch_index] = CTRL.idq[0]
            # watch_data[ 9][watch_index] = CTRL.idq[1]
            # watch_data[10][watch_index] = divmod(CTRL.theta_d, 2*np.pi)[1]
            # watch_data[11][watch_index] = CTRL.omega_r_elec / (2*np.pi*ACM.npp) * 60
            # watch_data[12][watch_index] = CTRL.cmd_rpm
            # watch_data[13][watch_index] = CTRL.cmd_idq[0]
            # watch_data[14][watch_index] = CTRL.cmd_idq[1]
            # watch_data[15][watch_index] = CTRL.xS[0] # theta_d
            # watch_data[16][watch_index] = CTRL.xS[1] / (2*np.pi*ACM.npp) * 60 # omega_r_elec
            # watch_data[17][watch_index] = CTRL.xS[2] # TL
            # watch_data[18][watch_index] = CTRL.xS[3] # pT
            # watch_data[19][watch_index] = CTRL.KA
            # watch_data[20][watch_index] = CTRL.KE
            # watch_data[21][watch_index] = CTRL.xT[0] # stator flux[0]
            # watch_data[22][watch_index] = CTRL.xT[1] # stator flux[1]
            # watch_data[23][watch_index] = CTRL.xT[2] # I term
            # watch_data[24][watch_index] = CTRL.xT[3] # I term
            # watch_data[25][watch_index] = 0.0 # CTRL.active_flux[0] # active flux[0]
            # watch_data[26][watch_index] = 0.0 # CTRL.active_flux[1] # active flux[1]
            # watch_data[27][watch_index] = CTRL.Tem
            # watch_data[28][watch_index] = CTRL.cmd_uab[0]
            # watch_data[29][watch_index] = CTRL.cmd_uab[1]
            # watch_data[30][watch_index] = svgen1.S1
            # watch_data[31][watch_index] = svgen1.EPwm1Regs_CMPA_bit_CMPA
            # watch_data[32][watch_index] = svgen1.Ta
            # watch_data[33][watch_index] = svgen1.S4
            # watch_data[34][watch_index] = svgen1.EPwm2Regs_CMPA_bit_CMPA
            # watch_data[35][watch_index] = svgen1.carrier_counter
            # watch_index += 1

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

        watch_data[28][watch_index] = ACM.udq[0] # -svgen1.line_to_line_voltage_AC # ACM.uab[0] # CTRL.cmd_uab[0]
        watch_data[29][watch_index] = ACM.udq[1] # svgen1.line_to_line_voltage_BC # svgen1.carrier_counter # CTRL.cmd_uab[1]

        watch_data[30][watch_index] = 30+svgen1.voltage_potential_at_terminal[0] # -svgen1.line_to_line_voltage_AC # ACM.uab[0]
        watch_data[31][watch_index] = svgen1.voltage_potential_at_terminal[1] # svgen1.line_to_line_voltage_BC # ACM.uab[1]
        watch_data[32][watch_index] = -30+svgen1.voltage_potential_at_terminal[2] # svgen1.line_to_line_voltage_AB # svgen1.line_to_line_voltage_BC
        watch_data[33][watch_index] = 0.0 # svgen1.voltage_potential_at_terminal[0] # svgen1.deadtime_counter[0] # svgen1.voltage_potential_at_terminal[0]
        watch_data[34][watch_index] = ACM.uab[0] # svgen1.deadtime_counter[1] # svgen1.voltage_potential_at_terminal[1]
        watch_data[35][watch_index] = ACM.uab[1] # svgen1.deadtime_counter[2] # svgen1.voltage_potential_at_terminal[2]
        watch_index += 1

    # return machine_times, watch_data # old
    return machine_times, watch_data # new



# TODO: need to make this globally shared between the simulation and the GUI.
_Watch_Mapping = [
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
    '[V]=CTRL.uab[0]', 
    '[V]=CTRL.uab[1]', 
    '[1]=svgen1.S1', 
    '[1]=svgen1.S2', 
    '[1]=svgen1.S3', 
    '[1]=svgen1.S4', 
    '[1]=svgen1.S5', 
    '[1]=svgen1.S6', 
]
Watch_Mapping = [el[el.find('=')+1:] for el in _Watch_Mapping]

############################################# Wrapper level 2
def ACMSimPyWrapper(numba__scope_dict, *arg, **kwarg):

    # Do Numerical Integrations (that do not care about numba__scope_dict at all and return watch_data whatsoever)
    machine_times, watch_data = ACMSimPyIncremental(*arg, **kwarg)
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



# Test incremental simulation
if __name__ == '__main__':

    # Basic settings
    CL_TS      = 1e-4 # [sec]
    TIME_SLICE = 0.2  # [sec]

    # init
    CTRL = The_Motor_Controller(CL_TS, 5*CL_TS,
                # init_npp = 4,
                # init_IN = 3,
                # init_R = 1.1,
                # init_Ld = 5e-3,
                # init_Lq = 6e-3,
                # init_KE = 0.095,
                # init_Rreq = 0.0, # note division by 0 is equal to infinity
                # init_Js = 0.0006168)
            init_npp = 21,
            init_IN = 72/1.414,
            init_R = 0.1222,
            init_Ld = 0.000502,
            init_Lq = 0.000571,
            init_KE = 0.188492, # 150 / 1.732 / (450/60*6.28*21)
            init_Rreq = 0.0, # -1.0, # PMSM
            init_Js = 0.203)
            # {'DOWN_SAMPLE': 1, 
            # 'n_pp': 21, 
            # 'IN': 72, 
            # 'PW': 14000, 
            # 'RPM': 470, 
            # 'J_s': 0.20300000000000004, 
            # 'Udc': 150, 
            # 'Rs': 0.1222, 
            # 'Ld': 0.0005020000000000001, 
            # 'Lq': 0.000571, 
            # 'KE': 0.188492, 
            # 'Ls': 0.000571, 
            # 'CL_TS': 0.0001, 
            # 'VL_TS': 0.0005, 
            # 'EndTime': 25, 
            # 'JLoadRatio': 0.0, 
            # 'Tload': 3, 
            # 'ViscousCoeff': 0.0007, 
            # 'data_file_name_prefix': 'SlessInv-OverEstimated-a2-NL'}
    # CTRL.bool_overwrite_speed_commands = False
    ACM       = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=100)

    # reg_id    = The_PI_Regulator(1*6.39955, 
    # 1*6.39955*237.845*CTRL.CL_TS, 150/1.732)
    # reg_iq    = The_PI_Regulator(1*6.39955, 1*6.39955*237.845*CTRL.CL_TS, 150/1.732)
    # # reg_speed = The_PI_Regulator(1.0*0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)
    # # reg_speed = The_PI_Regulator(0.1*0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)
    # # reg_speed = The_PI_Regulator(10 *0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)
    # # reg_speed = The_PI_Regulator(100 *0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)
    # reg_speed = The_PI_Regulator(0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)

    reg_id    = The_PI_Regulator(0.737168, 0.737168*214.011*CTRL.CL_TS, 150/1.732)
    reg_iq    = The_PI_Regulator(0.737168, 0.737168*214.011*CTRL.CL_TS, 150/1.732)
    reg_speed = The_PI_Regulator(0.323363, 0.323363*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)

    from collections import OrderedDict as OD
    numba__scope_dict = OD([
        # Y Labels                        Signal Name of Traces
        (r'Speed [rpm]',                  ( 'CTRL.cmd_rpm', 'CTRL.omega_r_mech', 'CTRL.xS[1]'   ,) ),
        (r'Speed Error [rpm]',            ( 'CTRL.cmd_rpm - CTRL.omega_r_mech'                  ,) ),
        (r'Position [rad]',               ( 'ACM.theta_d', 'CTRL.theta_d', 'CTRL.xS[0]'         ,) ),
        (r'Position mech [rad]',          ( 'ACM.theta_d'                                       ,) ),
        (r'$q$-axis current [A]',         ( 'ACM.iQ', 'CTRL.cmd_idq[1]'                         ,) ),
        (r'$d$-axis current [A]',         ( 'ACM.iD', 'CTRL.cmd_idq[0]'                         ,) ),
        (r'K_{\rm Active} [A]',           ( 'ACM.KA', 'CTRL.KA'                                 ,) ),
        (r'Load torque [Nm]',             ( 'CTRL.xS[2]'                                        ,) ),
        (r'CTRL.iD [A]',                  ( 'CTRL.cmd_idq[0]', 'CTRL.idq[0]'                    ,) ),
        (r'CTRL.iQ [A]',                  ( 'CTRL.cmd_idq[1]', 'CTRL.idq[1]'                    ,) ),
        (r'CTRL.uab [V]',                 ( 'CTRL.uab[0]', 'CTRL.uab[1]'                        ,) ),
        (r'S [1]',                        ( 'svgen1.S1', 'svgen1.S2', 'svgen1.S3', 'svgen1.S4', 'svgen1.S5', 'svgen1.S6' ,) ),
    ])

    # Global data arrays
    global_arrays = [None] * 23
    global_machine_times = None

    # simulate to generate ? sec of data
    for ii in range(0, 2):

        # controller_commands(t=ii*TIME_SLICE)

        """perform animation step"""
        machine_times, numba__waveforms_dict = \
            ACMSimPyWrapper(numba__scope_dict,
                        t0=ii*TIME_SLICE, TIME=TIME_SLICE, 
                        ACM=ACM,
                        CTRL=CTRL,
                        reg_id=reg_id,
                        reg_iq=reg_iq,
                        reg_speed=reg_speed)

        def save_to_global(_global, _local):
            return _local if _global is None else np.append(_global, _local)
        global_machine_times = save_to_global(global_machine_times, machine_times)
        print(len(global_machine_times))

        global_index = 0
        for ylabel in [ r'Speed [rpm]',
                        r'K_{\rm Active} [A]',
                        r'$d$-axis current [A]',
                        r'Position [rad]',
                        r'Load torque [Nm]',
                        r'CTRL.iD [A]',
                        r'CTRL.iQ [A]',
                        r'CTRL.uab [V]',
                        r'S [1]'
                ]:
            print(ylabel)
            # cmd_rpm, ACM_speed, OB_speed = numba__waveforms_dict[r'Speed [rpm]']
            # x5, KA                       = numba__waveforms_dict[r'K_{\rm Active} [A]']
            # ACM_id, cmd_id               = numba__waveforms_dict[r'$d$-axis current [A]']
            # ACM_theta_d, CTRL_theta_d, OB_theta_d = numba__waveforms_dict[r'Position [rad]']
            # TL,                          = numba__waveforms_dict[r'Load torque [Nm]']
            # cmd_iD, iD                   = numba__waveforms_dict[r'CTRL.iD [A]']
            # cmd_iQ, iQ                   = numba__waveforms_dict[r'CTRL.iQ [A]']
            # ualfa, ubeta                 = numba__waveforms_dict[r'CTRL.uab [V]']
            # S1, S2, S3, S4, S5, S6       = numba__waveforms_dict[r'S [1]']
            for trace_index, local_trace_data in enumerate(numba__waveforms_dict[ylabel]):
                global_arrays[global_index] = save_to_global(global_arrays[global_index], local_trace_data)
                global_index += 1
                print('\t', global_index, numba__scope_dict[ylabel][trace_index])

        # print(ii, f'{ACM.KA=} [Wb]', f'{CTRL.KA=} [Wb]', CTRL.cmd_rpm, '[rpm]')

        # for k,v in numba__waveforms_dict.items():
        #     print(k, np.shape(v))

        # print(len(global_speed), end='|')
        # print(max(global_speed), end='|')
        # print(len(speed), end='|')
        # print(max(speed))
        # print()
        # break
        pass

    global_cmd_speed   , \
    global_ACM_speed   , \
    global__OB_speed   , \
    global________x5   , \
    global___KA        , \
    global____ACM_id   , \
    global___CTRL_id   , \
    global_ACM_theta_d , \
    global_CTRL_theta_d, \
    global_OB_theta_d  , \
    global_TL          , \
    global_cmd_iD      , \
    global_iD          , \
    global_cmd_iQ      , \
    global_iQ          , \
    global_ualfa       , \
    global_ubeta       , \
    global_S1          , \
    global_S2          , \
    global_S3          , \
    global_S4          , \
    global_S5          , \
    global_S6          = global_arrays

    plt.figure(101, figsize=(15,4))
    plt.plot(global_cmd_speed); plt.plot(global_ACM_speed); plt.plot(global__OB_speed)

    plt.figure(2, figsize=(15,4)); #plt.ylim([-1e-1, 1e-1])
    plt.plot( (global_cmd_speed - global_ACM_speed) ) # [2000:5000]

    plt.figure(3, figsize=(15,4))
    plt.plot( global________x5 ); plt.plot( global___KA )

    plt.figure(4, figsize=(15,4))
    plt.plot( global____ACM_id ); plt.plot( global___CTRL_id )

    plt.figure(5, figsize=(15,4))
    plt.plot(global_ACM_theta_d); plt.plot( global_CTRL_theta_d); plt.plot( global_OB_theta_d )

    # plt.figure(figsize=(15,4))
    # plt.plot( np.sin(global_ACM_theta_d - global_OB_theta_d) )

    # plt.figure(figsize=(15,4))
    # plt.plot( global_TL )

    plt.figure(10, figsize=(15,4))
    plt.plot( global_cmd_iD); plt.plot( global_iD)
    plt.figure(11, figsize=(15,4))
    plt.plot( global_cmd_iQ); plt.plot( global_iQ)

    plt.figure(21, figsize=(15,4))
    plt.plot(global_ualfa, '--'); plt.plot(global_ubeta, '-.')
    # plt.figure(22, figsize=(15,4))
    # plt.plot(global_S1); plt.plot(global_S2)
    # plt.figure(23, figsize=(15,4))
    # plt.plot(global_S3);
    plt.figure(24, figsize=(15,4))
    plt.plot(global_S1)
    plt.plot(global_S2)
    plt.plot(global_S3)
    plt.figure(27, figsize=(15,4))
    plt.plot(global_S4); plt.plot(global_S5); plt.plot(global_S6)


    # print(CTRL.ell1, CTRL.ell2, CTRL.ell3, CTRL.ell4)

    plt.show()

