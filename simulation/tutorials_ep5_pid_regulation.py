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
            # feedback / input
            ('displacement', float64[:]),
            ('velocity', float64[:]),
            ('armature_current', float64[:]),
            # states
            ('timebase', float64),
            # commands
            ('cmd_force', float64[:]),
            ('cmd_armature_current', float64[:]),
        # MOTOR
    ])
class The_Motor_Controller:
    def __init__(self, CL_TS):
        ''' CONTROL '''
        # constants
        self.CL_TS = CL_TS
        # feedback / input
        self.displacement = np.zeros(2,dtype=np.float64)
        self.velocity = np.zeros(2,dtype=np.float64)
        # commands
        self.armature_current = np.zeros(2,dtype=np.float64)
        self.cmd_force = np.zeros(2, dtype=np.float64)
        self.cmd_armature_current = np.zeros(2, dtype=np.float64)

@jitclass(
    spec=[
        # states
        ('NS',    int32),
        ('x',   float64[:]),
        # inputs
        # output
        ('LIMIT_X', float64),
        ('LIMIT_Y', float64),
        ('UMP_X', float64),
        ('UMP_Y', float64),
        ('displacement', float64[:]),
        ('velocity', float64[:]),
        ('armature_current', float64[:]),
        # simulation settings
        ('MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD', int32),
    ])
class The_AC_Machine:
    def __init__(self, CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=1):
        # states
        self.NS = 6
        self.x = np.zeros(self.NS, dtype=np.float64)
        self.x[0] = 0.002
        self.x[1] = -0.002
        # inputs
        # output
        self.LIMIT_X = 0.0
        self.LIMIT_Y = 0.0
        self.UMP_X = 0.0
        self.UMP_Y = 0.0
        self.displacement = np.zeros(2, dtype=np.float64)
        self.velocity = np.zeros(2, dtype=np.float64)
        self.armature_current = np.zeros(2, dtype=np.float64)
        # simulation settings
        self.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD = MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD

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

############################################# MACHINE SIMULATION SECTION
@njit(nogil=True)
def DYNAMICS_MACHINE(t, x, ACM, CLARKE_TRANS_TORQUE_GAIN=1.5):
    fx = np.zeros(ACM.NS) # s x = f(x)

    ACM.LIMIT_X = 3.3*10e5 * (np.abs(x[0])-0.005) + 34*x[2]*(np.abs(x[0])-0.005) if np.abs(x[0])>0.005 else 0.0
    ACM.LIMIT_Y = 3.3*10e5 * (np.abs(x[1])-0.005) + 34*x[3]*(np.abs(x[1])-0.005) if np.abs(x[1])>0.005 else 0.0

    # 机械子系统
        # ACM.displacement[0] = ACM.x[0]
        # ACM.displacement[1] = ACM.x[1]
        # ACM.velocity[0]     = ACM.x[2]
        # ACM.velocity[1]     = ACM.x[3]
        # ACM.armature_current[0] = ACM.x[5]
        # ACM.armature_current[1] = ACM.x[6]
    fx[0] = x[2]
    fx[1] = x[3]
    fx[2] = (-0.02*x[2] + 9 * x[4]* 2.16e-5/(0.006-x[0])**2 - ACM.UMP_X + np.sign(x[0])*ACM.LIMIT_X) / 0.00837
    fx[3] = (-0.02*x[3] + 9 * x[5]* 2.16e-5/(0.006-x[1])**2 - ACM.UMP_Y + np.sign(x[1])*ACM.LIMIT_Y) / 0.00837
    fx[4] = 0.0
    fx[5] = 0.0

    return fx

@njit(nogil=True)
def RK4_MACHINE(t, ACM, hs): # 四阶龙格库塔法
    NS = ACM.NS
    k1, k2, k3, k4 = np.zeros(NS), np.zeros(NS), np.zeros(NS), np.zeros(NS) # incrementals at 4 stages
    xk, fx = np.zeros(NS), np.zeros(NS) # state x for stage 2/3/4, state derivative

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

############################################# DSP SECTION
@njit(nogil=True)
def DSP(ACM, CTRL, reg_dispX, reg_dispY):
    CTRL.timebase += CTRL.CL_TS

    """ Measurement """
    CTRL.displacement[0] = ACM.displacement[0]
    CTRL.displacement[1] = ACM.displacement[1]

    """ Speed Estimation """
    # CTRL.velocity[0] = ACM.velocity[0]
    # CTRL.velocity[1] = ACM.velocity[1]

    """ Position Controller (PD control) """
    reg_dispX.setpoint = 0.0
    reg_dispX.measurement = CTRL.displacement[0]
    tustin_pid(reg_dispX)
    reg_dispY.setpoint = 0.0
    reg_dispY.measurement = CTRL.displacement[1]
    tustin_pid(reg_dispY)

    # [$] Inverse Park transformation: get voltage commands in alpha-beta frame as SVPWM input
    CTRL.cmd_force[0] = reg_dispX.Out
    CTRL.cmd_force[1] = reg_dispY.Out

############################################# Wrapper level 1 (Main simulation | Incremental Edition)
""" MAIN for Real-time simulation """
@njit(nogil=True)
def ACMSimPyIncremental(t0, TIME, ACM=None, CTRL=None, reg_dispX=None, reg_dispY=None):

    # RK4 simulation and controller execution relative freuqencies
    MACHINE_TS = CTRL.CL_TS / ACM.MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD
    controller_down_sampling_ceiling = int(CTRL.CL_TS / MACHINE_TS)

    # watch variabels
    machine_times = np.arange(t0, t0+TIME, MACHINE_TS)
    watch_data    = np.zeros( (50, len(machine_times)) ) # new

    # Main loop
    jj = controller_down_sampling_ceiling # run controller at step 1
    watch_index = 0
    for ii in range(len(machine_times)):

        t = machine_times[ii]

        """ Machine Simulation @ MACHINE_TS """
        # Numerical Integration (ode4) with 5 states
        RK4_MACHINE(t, ACM, hs=MACHINE_TS)

        """ Machine Simulation Output @ MACHINE_TS """
        # Generate output variables for easy access
        ACM.displacement[0] = ACM.x[0]
        ACM.displacement[1] = ACM.x[1]
        ACM.velocity[0]     = ACM.x[2]
        ACM.velocity[1]     = ACM.x[3]
        ACM.armature_current[0] = ACM.x[4]
        ACM.armature_current[1] = ACM.x[5]

        jj += 1
        if jj >= controller_down_sampling_ceiling:
            jj = 0

            """ DSP @ CL_TS """
            # print(ii+1)
            DSP(ACM=ACM,
                CTRL=CTRL,
                reg_dispX=reg_dispX,
                reg_dispY=reg_dispY)

        # cmd to actual motor
        ACM.x[4] = CTRL.cmd_force[0]
        ACM.x[5] = CTRL.cmd_force[1]

        """ Watch @ MACHINE_TS """
        watch_data[ 0][watch_index] = ACM.displacement[0]
        watch_data[ 1][watch_index] = ACM.displacement[1]
        watch_data[ 2][watch_index] = ACM.velocity[0]
        watch_data[ 3][watch_index] = ACM.velocity[1]
        watch_data[ 4][watch_index] = ACM.armature_current[0]
        watch_data[ 5][watch_index] = ACM.armature_current[1]
        watch_data[ 6][watch_index] = CTRL.cmd_force[0]
        watch_data[ 7][watch_index] = CTRL.cmd_force[1]
        watch_data[ 8][watch_index] = ACM.UMP_X
        watch_index += 1

    return machine_times, watch_data # new



############################################# Wrapper level 2 (Collect waveforms data based off user specified names)
# TODO: need to make this globally shared between the simulation and the GUI.
_Unit_Watch_Mapping = [
    '[m]=ACM.displacement[0]',
    '[m]=ACM.displacement[1]',
    '[m/s]=ACM.velocity[0]',
    '[m/s]=ACM.velocity[1]',
    '[A]=ACM.armature_current[0]',
    '[A]=ACM.armature_current[1]',
    '[N]=CTRL.cmd_force[0]',
    '[N]=CTRL.cmd_force[1]',
    '[N]=ACM.UMP_X',
]
Watch_Mapping = [el[el.find('=')+1:] for el in _Unit_Watch_Mapping] # remove units before "="

def ACMSimPyWrapper(numba__scope_dict, *arg, **kwarg):

    # Do Numerical Integrations (that do not care about numba__scope_dict at all and return watch_data whatsoever)
    machine_times, watch_data = ACMSimPyIncremental(*arg, **kwarg)
    watch_data_as_dict = dict(zip(Watch_Mapping, watch_data))

    # Post-processing
    numba__waveforms_dict = dict()
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
    return machine_times, numba__waveforms_dict



############################################# Wrapper level 3 (User Interface)
from collections import OrderedDict as OD
class Simulation_Benchmark:
    def __init__(self, d, tuner=None, bool_start_simulation=True):

        self.d = d

        # 这个字典决定你要画的波形是哪些信号，具体能用的信号见：_Watch_Mapping
        # 允许信号之间进行简单的加减乘除，比如：'CTRL.cmd_rpm - CTRL.omega_r_mech'
        numba__scope_dict = OD([
            # Y Labels                        Signal Name of Traces
            (r'disp [m]', ('ACM.displacement[0]', 'ACM.displacement[1]', ) ),
            (r'velo [m/s]', ('ACM.velocity[0]', 'ACM.velocity[1]',) ),
            (r'curr [A]', ('ACM.armature_current[0]', 'ACM.armature_current[1]',) ),
            (r'force [N]', ('CTRL.cmd_force[0]', 'CTRL.cmd_force[1]',) ),
            (r'UMP [N]', ('ACM.UMP_X',) ),
        ])

        if bool_start_simulation:
            self.start_simulation_slices(d, numba__scope_dict)

    def get_global_objects(self):
        d = self.d
        # init
        CTRL = The_Motor_Controller(d['CL_TS'])
        ACM  = The_AC_Machine(CTRL, MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD=d['MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD'])
        reg_dispX = The_PID_Regulator(d['disp.Kp'], d['disp.Ki'], d['disp.Kd'], d['disp.tau'], d['disp.OutLimit'], d['disp.IntLimit'], d['CL_TS'])
        reg_dispY = The_PID_Regulator(d['disp.Kp'], d['disp.Ki'], d['disp.Kd'], d['disp.tau'], d['disp.OutLimit'], d['disp.IntLimit'], d['CL_TS'])
        return CTRL, ACM, reg_dispX, reg_dispY

    def start_simulation_slices(self, d, numba__scope_dict):

        global_objects = self.get_global_objects()
        self.CTRL, self.ACM, self.reg_dispX, self.reg_dispY = CTRL, ACM, reg_dispX, reg_dispY = global_objects

        global_trace_names = []
        max_number_of_traces = 0
        for ylabel, trace_names in numba__scope_dict.items():
            for name in trace_names:
                max_number_of_traces += 1
            for trace_index, name in enumerate(trace_names):
                global_trace_names.append(name)
        print(f'{max_number_of_traces=}')
        # print(global_trace_names)

        # init global data arrays for plotting
        global_arrays = [None] * max_number_of_traces
        global_machine_times = None

        def save_to_global(_global, _local):
            return _local if _global is None else np.append(_global, _local)

        # simulate to generate NUMBER_OF_SLICES*TIME_SLICE sec of data
        for ii in range(d['NUMBER_OF_SLICES']):

            # exec(d['user_system_input_code']) # 和 CONSOLE.user_controller_commands 功能相同

            # perform animation step
            machine_times, numba__waveforms_dict = \
                ACMSimPyWrapper(numba__scope_dict,
                            t0=ii*d['TIME_SLICE'], TIME=d['TIME_SLICE'],
                            ACM=ACM,
                            CTRL=CTRL,
                            reg_dispX=reg_dispX,
                            reg_dispY=reg_dispY)

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

        # print(f'{gdd.keys()=}')

        self.global_machine_times = global_machine_times
        self.gdd = gdd

if __name__ == '__main__':
    # User input:
    d = d_user_input_motor_dict = {
        'CL_TS': 1e-4,
        'TIME_SLICE': 0.2,
        'NUMBER_OF_SLICES': 2,
        'VL_EXE_PER_CL_EXE': 1,
        'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
        'disp.Kp': 5,
        'disp.Ki': 0.0,
        'disp.Kd': 20,
        'disp.tau': 0.005, # 0.02
        'disp.OutLimit': 10,
        'disp.IntLimit': 10,
    }

    图 = 1 # 空载加速、加载
    if True:
        sim1 = Simulation_Benchmark(d); gdd, global_machine_times = sim1.gdd, sim1.global_machine_times
        def 图1画图代码():
            plt.style.use('bmh') # https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
            mpl.rc('font', family='Times New Roman', size=10.0)
            mpl.rc('legend', fontsize=10)
            mpl.rcParams['lines.linewidth'] = 0.75 # mpl.rc('lines', linewidth=4, linestyle='-.')
            mpl.rcParams['mathtext.fontset'] = 'stix'

            fig, axes = plt.subplots(nrows=3, ncols=1, dpi=150, facecolor='w', figsize=(8,12), sharex=True)

            ax = axes[0]
            ax.plot(global_machine_times, gdd['ACM.displacement[0]'], label=r'$x$')
            ax.plot(global_machine_times, gdd['ACM.displacement[1]'], label=r'$y$')
            ax.set_ylabel(r'Disp [m]', multialignment='center') #) #, fontdict=font)

            ax = axes[1]
            ax.plot(global_machine_times, gdd['ACM.armature_current[0]'], label=r'$i_x$')
            ax.plot(global_machine_times, gdd['ACM.armature_current[1]'], label=r'$i_y$')
            # ax.plot(global_machine_times, gdd['CTRL.cmd_force[0]'], label=r'$i_x^*$')
            # ax.plot(global_machine_times, gdd['CTRL.cmd_force[1]'], label=r'$i_y^*$')
            ax.set_ylabel(r'Current [A]', multialignment='center') #) #, fontdict=font)

            ax = axes[2]
            ax.plot(global_machine_times, gdd['ACM.velocity[0]'], label=r'$v_x$')
            ax.plot(global_machine_times, gdd['ACM.velocity[1]'], label=r'$v_y$')
            ax.set_ylabel(r'Velocity [m/s]', multialignment='center') #) #, fontdict=font)

            for ax in axes:
                ax.grid(True)
                ax.legend(loc=1)
                # for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
                #     tick.label.set_font(font)
            axes[-1].set_xlabel('Time [s]') #, fontdict=font)
            return fig
        fig = 图1画图代码(); fig.savefig(f'SliceFSPM-fig-{图}.pdf', dpi=400, bbox_inches='tight', pad_inches=0)

    plt.show()
