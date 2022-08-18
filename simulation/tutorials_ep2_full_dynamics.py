# %%
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
            # states
            ('timebase', float64),
            ('KA', float64),
            ('Tem', float64),
            # commands
            ('cmd_idq', float64[:]),
            ('cmd_udq', float64[:]),
            ('cmd_uab', float64[:]),
            ('cmd_rpm', float64),
            ('index_separate_speed_estimation', int32),
            ('use_disturbance_feedforward_rejection', int32),
            # commands (sweep freuqency)
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
        init_Rreq = None, # note division by 0 is equal to infinity
        init_Js = 0.0006168,
    ):
        ''' CONTROL '''
        # constants
        self.CL_TS = CL_TS
        self.VL_TS = VL_TS
        self.velocity_loop_counter = 0 # 4
        self.velocity_loop_ceiling = 5
        # feedback / input
        self.theta_d = 0.0
        self.omega_r_elec = 0.0
        self.omega_syn = 0.0
        # states
        self.timebase = 0.0
        self.KA = init_KE
        self.Tem = 0.0
        # commands 
        self.cmd_idq = np.zeros(2, dtype=np.float64)
        self.cmd_udq = np.zeros(2, dtype=np.float64)
        self.cmd_uab = np.zeros(2, dtype=np.float64)
        self.cmd_rpm = 0.0
        self.index_separate_speed_estimation = 0
        self.use_disturbance_feedforward_rejection = 0
        # sweep frequency
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
        ('Tem', float64),
        ('cosT', float64),
        ('sinT', float64),
    ])
class The_AC_Machine:
    def __init__(self, CTRL):
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
        self.Tem = 0.0
        self.cosT = 1.0
        self.sinT = 0.0

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

    # CTRL.uab[0] = CTRL.uab_prev[0]
    # CTRL.uab[1] = CTRL.uab_prev[1]
    # CTRL.iab[0] = CTRL.iab_prev[0]
    # CTRL.iab[1] = CTRL.iab_prev[1]
    fx = THE_DYNAMICS(x, CTRL)
    for i in range(0, NS):
        k1[i] = fx[i] * hs
        xk[i] = x[i] + k1[i]*0.5

    # CTRL.iab[0] = 0.5*(CTRL.iab_prev[0]+CTRL.iab_curr[0])
    # CTRL.iab[1] = 0.5*(CTRL.iab_prev[1]+CTRL.iab_curr[1])
    # CTRL.uab[0] = 0.5*(CTRL.uab_prev[0]+CTRL.uab_curr[0])
    # CTRL.uab[1] = 0.5*(CTRL.uab_prev[1]+CTRL.uab_curr[1])
    fx = THE_DYNAMICS(xk, CTRL)
    for i in range(0, NS):
        k2[i] = fx[i] * hs
        xk[i] = x[i] + k2[i]*0.5

    fx = THE_DYNAMICS(xk, CTRL)
    for i in range(0, NS):
        k3[i] = fx[i] * hs
        xk[i] = x[i] + k3[i]

    # CTRL.iab[0] = CTRL.iab_curr[0]
    # CTRL.iab[1] = CTRL.iab_curr[1]
    # CTRL.uab[0] = CTRL.uab_curr[0]
    # CTRL.uab[1] = CTRL.uab_curr[1]
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
    CTRL.idq[0] = ACM.iD
    CTRL.idq[1] = ACM.iQ
    CTRL.theta_d = ACM.theta_d

    """ Park Transformation Essentials """
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
        # CTRL.iab_prev[0] = CTRL.iab_curr[0]
        # CTRL.iab_prev[1] = CTRL.iab_curr[1]
        # # CTRL.uab_prev[0] = CTRL.uab_curr[0] # This is needed only if voltage is measured, e.g., by eCAP.
        # # CTRL.uab_prev[1] = CTRL.uab_curr[1] # This is needed only if voltage is measured, e.g., by eCAP.

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
    MACHINE_TS = CTRL.CL_TS
    down_sampling_ceiling = int(CTRL.CL_TS / MACHINE_TS) #print('\tdown sample:', down_sampling_ceiling)

    # watch variabels
    machine_times  = np.arange(t0, t0+TIME, MACHINE_TS)
    control_times  = np.arange(t0, t0+TIME, CTRL.CL_TS)
    watch_data = np.zeros( (30, len(control_times)) )

    # Main loop
    # print('\tt0 =', t0)
    jj = 0; watch_index = 0
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
        # ACM.cosT = np.cos(ACM.theta_d)
        # ACM.sinT = np.sin(ACM.theta_d)
        # ACM.iab[0] = ACM.x[0] * ACM.cosT + ACM.x[1] *-ACM.sinT
        # ACM.iab[1] = ACM.x[0] * ACM.sinT + ACM.x[1] * ACM.cosT

        jj += 1
        if jj >= down_sampling_ceiling:
            jj = 0

            if CTRL.bool_overwrite_speed_commands == True:

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
            DSP(ACM=ACM,
                CTRL=CTRL,
                reg_speed=reg_speed,
                reg_id=reg_id,
                reg_iq=reg_iq)

            """ Console @ CL_TS """
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

            """ Watch @ CL_TS """
            watch_data[ 0][watch_index] = divmod(ACM.theta_d, 2*np.pi)[1]
            watch_data[ 1][watch_index] = ACM.omega_r_mech / (2*np.pi) * 60 # omega_r_mech
            watch_data[ 2][watch_index] = ACM.KA
            watch_data[ 3][watch_index] = ACM.iD
            watch_data[ 4][watch_index] = ACM.iQ
            watch_data[ 5][watch_index] = ACM.Tem
            watch_data[ 6][watch_index] = 0.0 # CTRL.iab[0]
            watch_data[ 7][watch_index] = 0.0 # CTRL.iab[1]
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
            watch_index += 1

            """ Voltage Source Inverter (in dq-frame) """
            ACM.udq[0] = CTRL.cmd_udq[0] # + current denpendent distorted voltage
            ACM.udq[1] = CTRL.cmd_udq[1] # + current denpendent distorted voltage

    return control_times, watch_data



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
]
Watch_Mapping = [el[el.find('=')+1:] for el in _Watch_Mapping]

############################################# Wrapper level 2
def ACMSimPyWrapper(numba__scope_dict, *arg, **kwarg):

    # Do Numerical Integrations (that do not care about numba__scope_dict at all and return watch_data whatsoever)
    control_times, watch_data = ACMSimPyIncremental(*arg, **kwarg)
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
                    if 'CTRL' in word or 'ACM' in word or 'reg_' in word:
                        translated_expression += f'watch_data_as_dict["{word}"]'
                    else:
                        translated_expression += word
                # print(translated_expression)
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
    return control_times, numba__waveforms_dict



# Test incremental simulation
if __name__ == '__main__':

    # Basic settings
    CL_TS      = 1e-4 # [sec]
    TIME_SLICE = 0.5  # [sec]

    # init
    CTRL = The_Motor_Controller(CL_TS, 5*CL_TS,
                init_npp = 4,
                init_IN = 3,
                init_R = 1.1,
                init_Ld = 5e-3,
                init_Lq = 6e-3,
                init_KE = 0.095,
                init_Rreq = -1.0, # PMSM
                # init_Rreq = 1.0, # IM
                init_Js = 0.0006168)
    CTRL.bool_overwrite_speed_commands = False
    ACM       = The_AC_Machine(CTRL)
    reg_id    = The_PI_Regulator(6.39955, 6.39955*237.845*CTRL.CL_TS, 600)
    reg_iq    = The_PI_Regulator(6.39955, 6.39955*237.845*CTRL.CL_TS, 600)
    reg_speed = The_PI_Regulator(1.0*0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)
    # reg_speed = The_PI_Regulator(10 *0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)
    # reg_speed = The_PI_Regulator(0.1*0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)

    # Global arrays
    global_cmd_speed, global_ACM_speed, global__OB_speed = None, None, None
    global___KA = None
    global____ACM_id, global___CTRL_id = None, None
    global________x5 = None
    global_ACM_theta_d, global_CTRL_theta_d, global_OB_theta_d = None, None, None
    global_TL = None

    from collections import OrderedDict as OD
    numba__scope_dict = OD([
        (r'Speed [rpm]',                  ( 'CTRL.cmd_rpm', 'CTRL.omega_r_mech', 'CTRL.xS[1]'   ,) ),
        (r'Speed Error [rpm]',            ( 'CTRL.cmd_rpm - CTRL.omega_r_mech'                  ,) ),
        (r'Position [rad]',               ( 'ACM.theta_d', 'CTRL.theta_d', 'CTRL.xS[0]'         ,) ),
        (r'Position mech [rad]',          ( 'ACM.theta_d'                                          ,) ),
        (r'$q$-axis current [A]',         ( 'ACM.iQ', 'CTRL.cmd_idq[1]'                       ,) ),
        (r'$d$-axis current [A]',         ( 'ACM.iD', 'CTRL.cmd_idq[0]'                       ,) ),
        (r'K_{\rm Active} [A]',           ( 'ACM.KA', 'CTRL.KA'                               ,) ),
        (r'Load torque [Nm]',             ( 'CTRL.xS[2]'                                        ,) ),
    ])

    # def controller_commands(t):
    #     """ Console @ CL_TS """
    #     if t <= 1.0:
    #         CTRL.cmd_rpm = 50
    #     elif t <= 1.5:
    #         ACM.TLoad = 2
    #     elif t <= 2.0:
    #         CTRL.cmd_rpm = 200
    #     elif t <= 3.0:
    #         CTRL.cmd_rpm = -200
    #     elif t <= 4.0:
    #         CTRL.cmd_rpm = 0
    #     elif t <= 4.5:
    #         CTRL.cmd_rpm = 2000
    #     elif t <= 5:
    #         CTRL.cmd_idq[0] = 2
    #     elif t <= 5.5:
    #         ACM.TLoad = 0.0
    #     elif t <= 6: 
    #         CTRL.CMD_SPEED_SINE_RPM = 500
    #     # else: # don't implement else clause to receive commands from IPython console

    #     if CTRL.CMD_SPEED_SINE_RPM!=0:
    #         CTRL.cmd_rpm = CTRL.CMD_SPEED_SINE_RPM * np.sin(2*np.pi*CTRL.CMD_SPEED_SINE_HZ*t)

    # simulate to generate 10 sec of data
    for ii in range(0, 8):

        # controller_commands(t=ii*TIME_SLICE)

        """perform animation step"""
        control_times, numba__waveforms_dict = \
            ACMSimPyWrapper(numba__scope_dict,
                        t0=ii*TIME_SLICE, TIME=TIME_SLICE, 
                        ACM=ACM,
                        CTRL=CTRL,
                        reg_id=reg_id,
                        reg_iq=reg_iq,
                        reg_speed=reg_speed)
        cmd_rpm, ACM_speed, OB_speed = numba__waveforms_dict[r'Speed [rpm]']
        x5, KA                       = numba__waveforms_dict[r'K_{\rm Active} [A]']
        ACM_id, cmd_id               = numba__waveforms_dict[r'$d$-axis current [A]']
        ACM_theta_d, CTRL_theta_d, OB_theta_d = numba__waveforms_dict[r'Position [rad]']
        TL,                          = numba__waveforms_dict[r'Load torque [Nm]']

        def save_to_global(_global, _local):
            return _local if _global is None else np.append(_global, _local)    
        global_cmd_speed = save_to_global(global_cmd_speed, cmd_rpm)
        global_ACM_speed = save_to_global(global_ACM_speed, ACM_speed)
        global__OB_speed = save_to_global(global__OB_speed,  OB_speed)
        global________x5 = save_to_global(global________x5,        x5)
        global___KA = save_to_global(global___KA,   KA)
        global____ACM_id = save_to_global(global____ACM_id,   ACM_id)
        global___CTRL_id = save_to_global(global___CTRL_id,   cmd_id)
        global_ACM_theta_d  = save_to_global(global_ACM_theta_d, ACM_theta_d)
        global_CTRL_theta_d = save_to_global(global_CTRL_theta_d, CTRL_theta_d)
        global_OB_theta_d = save_to_global(global_OB_theta_d, OB_theta_d)
        global_TL        = save_to_global(global_TL, TL)

        print(ii, 'KA =', ACM.KA, CTRL.KA, 'Wb', CTRL.cmd_rpm, 'rpm')

        # for k,v in numba__waveforms_dict.items():
        #     print(k, np.shape(v))

        # print(len(global_speed), end='|')
        # print(max(global_speed), end='|')
        # print(len(speed), end='|')
        # print(max(speed))
        # print()
        # break

    plt.figure(figsize=(15,4))
    plt.plot(global_cmd_speed); plt.plot(global_ACM_speed); plt.plot(global__OB_speed)

    plt.figure(figsize=(15,4)); #plt.ylim([-1e-1, 1e-1])
    plt.plot( (global_cmd_speed - global_ACM_speed) ) # [2000:5000]

    plt.figure(figsize=(15,4))
    plt.plot( global________x5 ); plt.plot( global___KA )

    plt.figure(figsize=(15,4))
    plt.plot( global____ACM_id ); plt.plot( global___CTRL_id )

    plt.figure(figsize=(15,4))
    plt.plot(global_ACM_theta_d); plt.plot( global_CTRL_theta_d); plt.plot( global_OB_theta_d )

    plt.figure(figsize=(15,4))
    plt.plot( np.sin(global_ACM_theta_d - global_OB_theta_d) )

    # plt.figure(figsize=(15,4))
    # plt.plot( global_TL )

    print(CTRL.ell1, CTRL.ell2, CTRL.ell3, CTRL.ell4)

    plt.show()


