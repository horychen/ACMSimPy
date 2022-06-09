# 常用库

from concurrent.futures import thread
from pylab import np, mpl, plt
from collections import OrderedDict as OD
from dataclasses import dataclass
import matplotlib.animation as animation
import json, re, copy, pkg_resources
import traceback
import pandas as pd

# from gui.core.json_settings import Settings
# from simulation.tutorials import acmsimpy
import simulation.tutorials as acmsimpy

# 后端
# use cairo only for acmsimpy | use cairo for acmsimc will slow down plotting
# mpl.use('cairo') # try if this is faster for animate plotting? ['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']

# 风格
mpl.style.use('ggplot') # further customize: https://stackoverflow.com/questions/35223312/matplotlib-overriding-ggplot-default-style-properties and https://matplotlib.org/2.0.2/users/customizing.html
# mpl.style.use("dark_background")
# mpl.style.use('grayscale')
# mpl.style.use('classic')

# 字体
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['font.weight'] ='900'     # bold is 700, normal is 400, see https://matplotlib.org/2.0.2/users/customizing.html
mpl.rcParams['font.size']   =12 # 13 is too big, commenting out is too small

# 颜色设置壹/叁
# https://matplotlib.org/2.0.2/users/customizing.html from https://stackoverflow.com/questions/35223312/matplotlib-overriding-ggplot-default-style-properties
# mpl.rcParams['axes.labelsize']  =14   # 'medium' # 'large'是没用的！
hex_monokai   = '#272822'
hex_wanderson = '#3d444c'
mpl.rcParams['axes.facecolor']  =hex_wanderson # also need to change facecolor in mplwidget.py to get a consistent look
mpl.rcParams['axes.labelcolor'] ='#d5d1c7' # tint grey
mpl.rcParams['axes.axisbelow']  ='False'   # 'line'
mpl.rcParams['ytick.color']     ='#d5d1c7' #'white'
mpl.rcParams['xtick.color']     ='#d5d1c7' #'white'
mpl.rcParams['text.color']      ='#d5d1c7' #'white'
mpl.rcParams['grid.color']      ='#d5d1c7' # grid color
mpl.rcParams['grid.linestyle']  ='--'      # solid
mpl.rcParams['grid.linewidth']  =0.3       # in points
mpl.rcParams['grid.alpha']      =0.8       # transparency, between 0.0 and 1.0



@dataclass
class CJH_Plotting(object):
    cjh_linestyles = [
            '-', '--', (0, (3, 1, 1, 1)), ':', '-.', 
            '-', '--', (0, (3, 1, 1, 1)), ':', '-.', 
            '-', '--', (0, (3, 1, 1, 1)), ':', '-.', 
        ]
    # 颜色设置叁/叁
    cjh_colors = [
            '#ff6347', # tomato red
            'white',   # #d5d1c7',
            '#B3F5FF', # dark-theme-bright-blue
            '#FFF0B3', # dark-theme-bright-yellow
            '#ffabff', # pink
            '#E1C7E0', # dark-theme-bright-violet(purple)
            '#ABF5D1', # dark-theme-bright-green
            '#FFBDAD', # dark-theme-bright-red
            '#00dcff', # blue
            '#c593ff', # dark-theme-bright-color1
            '#02dac3', # dark-theme-bright-color2
            '#efc9a1', # dark-theme-bright-color3
            # 'black',
            'blue',
            '#857F72', # warmGrey
            '#616E7C', # Cool Grey
            '#B8B2A7', # warmGrey
            '#9AA5B1', # Cool Grey
            '#504A40', # warmGrey
            '#27241D', # warmGrey
        ] # https://uxdesign.cc/dark-mode-ui-design-the-definitive-guide-part-1-color-53dcfaea5129

    def __init__(self):
        # 作图用外观设置
        # generator的好处是循环，列表则是有限的 # https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/linestyles.html
        self.cjh_linestyle = self.cyclic_generator(self.cjh_linestyles)
        self.cjh_color = self.cyclic_generator(self.cjh_colors)
        # 注意，不可以做 list(self.cjh_color)，因为这是一个无止境循环的发生器，会卡住的。。。

    @staticmethod
    def cyclic_generator(list_of_things):
        N, i = len(list_of_things), 0
        while True:
            yield list_of_things[i]
            i +=1
            if i > (N-1):
                i = 0


@dataclass
class THE_CONSOLE:
    """ User control over the mpl animation """
    # offset_anime_ii : int = 0
    reset : int = False
    _pause : int = False
    SAMPLING_RATE : float = 1e4
    TIME_SLICE : float = 0.2
    NUMBER_OF_TIME_SLICE_TO_SHOW : int = 10
    bool_exit: int = False
    numba__scope_dict: dict = None
    ii_list: list = None
    def __post_init__(self):
        # self.pause_time = 0.0
        self.CL_TS = 1 / self.SAMPLING_RATE
        self.NUMBER_OF_SAMPLE_TO_SHOW = int(self.NUMBER_OF_TIME_SLICE_TO_SHOW * self.TIME_SLICE * self.SAMPLING_RATE)

# FUNCTIONS
class EmyFunctions(object):

    def prepare_canvas_on_page_3(mainWindowObject):

        """ Read scope dict from GUI """
        try:
            the_cmd = mainWindowObject.plainTextEdit_NumbaScopeDict.toPlainText()
            with open('user_input.txt', 'w') as f:
                f.write(the_cmd)
            numba__scope_dict = eval(the_cmd[the_cmd.find('OD'):])
        except Exception as err:
            print('-------------------')
            print('Decypher for NumbaScopeDict has failed. Will use the default dict.')
            numba__scope_dict = OD([
                (r'Speed [rpm]',                  ( 'CTRL.omega_elec'           ,)  ),
                (r'Position [rad]',               ( 'ACM.x[3]'                  ,)  ),
                (r'$q$-axis current [A]',         ( 'ACM.x[1]', 'CTRL.idq[1]'   ,)  ),
                (r'$d$-axis current [A]',         ( 'ACM.x[0]', 'CTRL.idq[0]'   ,)  ),
                (r'$\alpha\beta$ current [A]',    ( 'CTRL.iab[0]', 'CTRL.iab[1]',)  ),
            ])
            print('err:', err)
            traceback.print_exc()
            print('-------------------')

        # 颜色设置贰/叁
        # mainWindowObject.ui.MplWidget_ACMPlot.toolbar.setStyleSheet("background-color: #9AA5B1;")
        # mainWindowObject.ui.right_column.horizontalLayout_CBSMplToolBar.addWidget(ui.MplWidget_ACMPlot.toolbar)

        CONSOLE = THE_CONSOLE(numba__scope_dict=numba__scope_dict)
        mainWindowObject.console_push_variable({f'CONSOLE':CONSOLE})


        # print('numba__scope_dict =', numba__scope_dict)
        number_of_subplot = len(numba__scope_dict)
        numba__axes = []
        numba__line_dict = OD()
        first_ax = None

        cjh_plotting = CJH_Plotting()

        """ Canvas Behaviors """
        def onClick(event):
            print('Hint: double click to pause the simulation.', event)
            if event.dblclick == True:
                CONSOLE._pause ^= True
            # if event.button == 3:
            #     CONSOLE._pause ^= True
                # if CONSOLE._pause:
                #     # CONSOLE.anim.event_source.stop()
                #     CONSOLE.pause_time = 999
                # else:
                #     # CONSOLE.anim.event_source.start()
                #     CONSOLE.pause_time = 0.0
        mainWindowObject.ui.MplWidget_ACMPlot.canvas.figure.canvas.mpl_connect('button_press_event', onClick)
        # mainWindowObject.ui.MplWidget_ACMPlot.canvas.draw_idle() # draw before animation # https://stackoverflow.com/questions/8955869/why-is-plotting-with-matplotlib-so-slow

        """ Canvas """
        mainWindowObject.ui.MplWidget_ACMPlot.canvas.figure.clf() # clear canvas
        mainWindowObject.ui.MplWidget_ACMPlot.setMinimumSizeByNumberOfSubplots(number_of_subplot)

        for index, (ylabel, waveform_names) in enumerate(numba__scope_dict.items()):
            # get axes
            ax = mainWindowObject.ui.MplWidget_ACMPlot.canvas.figure.add_subplot(
                number_of_subplot, 
                1, 1+index, 
                autoscale_on=False, 
                sharex=first_ax
            )
            first_ax = ax if index == 0 else first_ax
            numba__axes.append(ax)

            numba__line_dict[ylabel] = []
            for jj, name in enumerate(waveform_names):
                # print(jj, name)
                # get line objects
                line, = ax.plot([], [], 
                                linestyle=cjh_plotting.cjh_linestyles[jj], 
                                color=cjh_plotting.cjh_colors[jj], 
                                lw=1.5, 
                                label=name, 
                                alpha=0.7) # zorder
                line.ax = ax
                line.ydata = np.array(0)
                line.MIN = -1e-10
                line.MAX =  1e-10
                numba__line_dict[ylabel].append(line)
            # print()

            ax.set_ylabel(ylabel)
            ax.legend(loc='lower left').set_zorder(202)

        numba__axes[-1].set_xlabel('Time [s]')
        # time_text = first_ax.text(0.02, 0.95, '', transform=first_ax.transAxes)

        mainWindowObject.CONSOLE = CONSOLE
        mainWindowObject.numba__line_dict = numba__line_dict
        mainWindowObject.first_ax = first_ax

    # Run real time simulation with ACMSymPy
    ''' Python-Numba-based Simulation '''
    # @Slot()
    # @decorator_status(status='Configuring') # this will not be updated until this function is done
    def runPyBasedSimulation(mainWindowObject, progress_callback=None):

        """ Initialization """
        EmyFunctions.prepare_canvas_on_page_3(mainWindowObject)
        CONSOLE = mainWindowObject.CONSOLE

        """ Simulation Globals """
        CTRL = mainWindowObject.CTRL = acmsimpy.The_Motor_Controller(CONSOLE.CL_TS, 5*CONSOLE.CL_TS,
                init_npp = 4,
                init_IN = 3,
                init_R = 1.1,
                # init_Ld = 5e-3, # PMSM
                # init_Lq = 6e-3, # PMSM
                # init_KE = 0.095, # PMSM
                # init_Rreq = -1.0, # PMSM
                init_Ld = 440e-3, # IM
                init_Lq = 25e-3, # IM
                init_KE = 0.0, # IM
                init_Rreq = 1.0, # IM
                init_Js = 0.0006168)
        ACM       = mainWindowObject.ACM       = acmsimpy.The_AC_Machine(CTRL)
        reg_id    = mainWindowObject.reg_id    = None
        reg_iq    = mainWindowObject.reg_iq    = None
        reg_speed = mainWindowObject.reg_speed = acmsimpy.The_PI_Regulator(0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)

        """ Simulation Globals Access from Console """
        if mainWindowObject.console_window is not None:
            mainWindowObject.console_push_variable({'CONSOLE':CONSOLE})
            mainWindowObject.console_push_variable({'CTRL':CTRL})
            mainWindowObject.console_push_variable({'ACM':ACM})
            mainWindowObject.console_push_variable({'reg_id':reg_id})
            mainWindowObject.console_push_variable({'reg_iq':reg_iq})
            mainWindowObject.console_push_variable({'reg_speed':reg_speed})

        """ Visualization of Realtime Simulation """
        print('\tJIT compile with numba...')
        ii = 0
        while ii<1000:
            ii += 1

            # start time for the present slice of simulation
            t0 = ii*CONSOLE.TIME_SLICE

            # Run one slice of simulation
            control_times, numba__waveforms_dict, = acmsimpy.ACMSimPyWrapper(
                CONSOLE.numba__scope_dict,
                t0=t0, TIME=CONSOLE.TIME_SLICE,
                ACM=ACM,
                CTRL=CTRL,
                reg_id=reg_id,
                reg_iq=reg_iq,
                reg_speed=reg_speed,
            )
            end_time = control_times[-1] # print('\tend_time:', end_time)
            # progress_callback.emit([t0, numba__waveforms_dict, end_time])

            """ update matplotlib artist """
            for key, waveforms_data in numba__waveforms_dict.items():
                ymin, ymax = None, None
                for line, ydata in zip(mainWindowObject.numba__line_dict[key], waveforms_data):
                    # print(key, line, ydata[0:3], mainWindowObject.numba__line_dict[key], len(waveforms_data))
                    local_ymin, local_ymax = EmyFunctions.update_line_data(CONSOLE, end_time, line, ydata)
                    if ymin is None: ymin=local_ymin; ymax=local_ymax
                    ymin = local_ymin if local_ymin < ymin else ymin
                    ymax = local_ymax if local_ymax > ymax else ymax
                # 完事了只对一个ax做一次
                if ymin != ymax:
                    line.ax.set_ylim([ymin, ymax]) # this .ax is manually assigned when line object is created.

            # time_text.set_text('time = %.1f' % end_time)

            mainWindowObject.first_ax.set_xlim([end_time-CONSOLE.NUMBER_OF_SAMPLE_TO_SHOW*CONSOLE.CL_TS, end_time])
            # first_ax.set_xticklabels(np.arange(0, int(round(end_time)), int(round(end_time))*0.1)) #, fontdict=font)
            # return line01, line02, time_text # for using blit=True but blit does not update xticks and yticks

            mainWindowObject.ui.MplWidget_ACMPlot.canvas.draw()
            # mainWindowObject.ui.MplWidget_ACMPlot.canvas.flush_events() # flush any pending GUI events, re-painting the screen if needed

            progress_callback.emit(t0)
            if CONSOLE.SAMPLING_RATE * CONSOLE.TIME_SLICE < 1e4*0.5:
                plt.pause(0.2) # avoid to call canvas.draw too frequently.
            while CONSOLE._pause:
                plt.pause(0.2)
                if CONSOLE.bool_exit == True:
                    break
            if CONSOLE.bool_exit == True:
                break

        # Animation cannot be used as the thread along with the anim object will be killed after exit here.
        return 'Simulation is done.' # the result to print

    def runPyBasedSimulationParallel(mainWindowObject, thread_index, extra_execution_codes, progress_callback=None):

        """ Initialization is moved to MainWindow.btn_clicked() """
        CONSOLE = mainWindowObject.CONSOLE # Note this is the shared CONSOLE object

        """ Simulation Globals """
        mainWindowObject.CTRL      = CTRL      = acmsimpy.The_Motor_Controller(CONSOLE.CL_TS, 5*CONSOLE.CL_TS,
                init_npp = 4,
                init_IN = 3,
                init_R = 1.1,
                # init_Ld = 5e-3, # PMSM
                # init_Lq = 6e-3, # PMSM
                # init_KE = 0.095, # PMSM
                # init_Rreq = -1.0, # PMSM
                init_Ld = 440e-3, # IM
                init_Lq = 25e-3, # IM
                init_KE = 0.0, # IM
                init_Rreq = 1.0, # IM
                init_Js = 0.0006168)
        mainWindowObject.ACM       = ACM       = acmsimpy.The_AC_Machine(CTRL)
        mainWindowObject.reg_id    = reg_id    = None
        mainWindowObject.reg_iq    = reg_iq    = None
        mainWindowObject.reg_speed = reg_speed = acmsimpy.The_PI_Regulator(0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)

        """ Execute extra codes for each thread here"""
        exec(extra_execution_codes)
        # exec("print(dir())")
        print('Parallel thread', thread_index, extra_execution_codes)

        """ Simulation Globals Access from Console """
        mainWindowObject.console_push_variable({f'CTRL{thread_index:d}':CTRL})
        mainWindowObject.console_push_variable({f'ACM{thread_index:d}':ACM})
        mainWindowObject.console_push_variable({f'reg_id{thread_index:d}':reg_id})
        mainWindowObject.console_push_variable({f'reg_iq{thread_index:d}':reg_iq})
        mainWindowObject.console_push_variable({f'reg_speed{thread_index:d}':reg_speed})

        """ Visualization of Realtime Simulation """
        print('\tJIT compile with numba...')
        ii = 0
        while ii<1000:
            ii += 1

            # start time for the present slice of simulation
            t0 = ii*CONSOLE.TIME_SLICE

            # Run one slice of simulation
            control_times, numba__waveforms_dict, = acmsimpy.ACMSimPyWrapper(
                CONSOLE.numba__scope_dict,
                t0=t0, TIME=CONSOLE.TIME_SLICE,
                ACM=ACM,
                CTRL=CTRL,
                reg_id=reg_id,
                reg_iq=reg_iq,
                reg_speed=reg_speed,
            )
            end_time = control_times[-1] # print('\tend_time:', end_time)
            # progress_callback.emit([t0, numba__waveforms_dict, end_time])

            """ update matplotlib artist """
            for key, waveforms_data in numba__waveforms_dict.items():
                ymin, ymax = None, None
                if True:
                    line, ydata = mainWindowObject.numba__line_dict[key][thread_index], waveforms_data[0]
                    # debug print
                    # print(thread_index, key, line, ydata[0:3], mainWindowObject.numba__line_dict[key], len(waveforms_data))
                    local_ymin, local_ymax = EmyFunctions.update_line_data(CONSOLE, end_time, line, ydata)
                    if ymin is None: ymin=local_ymin; ymax=local_ymax
                    ymin = local_ymin if local_ymin < ymin else ymin
                    ymax = local_ymax if local_ymax > ymax else ymax
                # 完事了只对一个ax做一次
                if ymin != ymax:
                    min_scale = 1.05 if ymin<0 else 0.95
                    max_scale = 1.05 if ymax>0 else 0.95
                    line.ax.set_ylim([ymin*min_scale, ymax*max_scale]) # this .ax is manually assigned when line object is created.

            # time_text.set_text('time = %.1f' % end_time)

            # Pause or Exit
            if CONSOLE.SAMPLING_RATE * CONSOLE.TIME_SLICE < 1e4*0.5:
                plt.pause(0.2) # avoid to call canvas.draw too frequently.
            while CONSOLE._pause:
                plt.pause(0.2)
                if CONSOLE.bool_exit == True:
                    break
            if CONSOLE.bool_exit == True:
                # CONSOLE = None # this causes problem for multiple threads to exit
                break

            # Synchronize the threads
            CONSOLE.ii_list[thread_index] = ii

            # Update the canvas
            if thread_index == 0:
                while sum([el == CONSOLE.ii_list[0] for el in CONSOLE.ii_list]) < len(CONSOLE.ii_list):
                    # print('Pause:', thread_index, CONSOLE.ii_list)
                    plt.pause(0.01)
                    if CONSOLE.bool_exit == True:
                        break
                # print(thread_index, CONSOLE.ii_list)
                mainWindowObject.first_ax.set_xlim([end_time-CONSOLE.NUMBER_OF_SAMPLE_TO_SHOW*CONSOLE.CL_TS, end_time])
                # first_ax.set_xticklabels(np.arange(0, int(round(end_time)), int(round(end_time))*0.1)) #, fontdict=font)
                # return line01, line02, time_text # for using blit=True but blit does not update xticks and yticks
                mainWindowObject.ui.MplWidget_ACMPlot2.canvas.draw()
                # mainWindowObject.ui.MplWidget_ACMPlot2.canvas.flush_events() # flush any pending GUI events, re-painting the screen if needed
                progress_callback.emit(t0)
            else:
                while ii > CONSOLE.ii_list[0]:
                    plt.pause(0.01)
                    if CONSOLE.bool_exit == True:
                        break
        # Animation cannot be used as the thread along with the anim object will be killed after exit here.
        return 'Simulation (parallel) is done.' # the result to print

    def prepare_canvas_on_page_4(mainWindowObject, number_of_threads):

        """ Read scope dict from GUI """
        try:
            the_cmd = mainWindowObject.ui.plainTextEdit_NumbaScopeDict_Parallel.toPlainText()
            with open('user_input_parallel.txt', 'w') as f:
                f.write(the_cmd)
            numba__scope_dict = eval(the_cmd[the_cmd.find('OD'):])
        except Exception as err:
            print('-------------------')
            print('Decypher for NumbaScopeDict has failed. Will use the default dict.')
            numba__scope_dict = OD([
                        (r'Speed [rpm]',      ( 'ACM.omega_r_mech',)),
                        (r'Position [rad]',   ( 'ACM.theta_d'   ,)),
                        (r'$dq$ current [A]', ( 'ACM.iQ'        ,)),
                        (r'Torque [Nm]',      ( 'ACM.Tem'       ,)),])
            print('err:', err)
            traceback.print_exc()
            print('-------------------')

        # 颜色设置贰/叁
        # mainWindowObject.ui.MplWidget_ACMPlot2.toolbar.setStyleSheet("background-color: #9AA5B1;")
        # mainWindowObject.ui.right_column.horizontalLayout_CBSMplToolBar.addWidget(ui.MplWidget_ACMPlot2.toolbar)

        CONSOLE = THE_CONSOLE(numba__scope_dict=numba__scope_dict, ii_list=[0 for _ in range(number_of_threads)])
        mainWindowObject.console_push_variable({f'CONSOLE':CONSOLE})

        # print('numba__scope_dict =', numba__scope_dict)
        number_of_subplot = len(numba__scope_dict)
        numba__axes = []
        numba__line_dict = OD()
        first_ax = None

        cjh_plotting = CJH_Plotting()

        """ Canvas Behaviors """
        def onClick(event):
            print('Hint: double click to pause the simulation.', event)
            if event.dblclick == True:
                CONSOLE._pause ^= True
        mainWindowObject.ui.MplWidget_ACMPlot2.canvas.figure.canvas.mpl_connect('button_press_event', onClick)
        # mainWindowObject.ui.MplWidget_ACMPlot2.canvas.draw_idle() # draw before animation # https://stackoverflow.com/questions/8955869/why-is-plotting-with-matplotlib-so-slow

        """ Canvas """
        mainWindowObject.ui.MplWidget_ACMPlot2.canvas.figure.clf() # clear canvas
        mainWindowObject.ui.MplWidget_ACMPlot2.setMinimumSizeByNumberOfSubplots(number_of_subplot)

        for index, (ylabel, waveform_names) in enumerate(numba__scope_dict.items()):
            # get axes
            ax = mainWindowObject.ui.MplWidget_ACMPlot2.canvas.figure.add_subplot(
                number_of_subplot, 
                1, 1+index, 
                autoscale_on=False, 
                sharex=first_ax
            )
            first_ax = ax if index == 0 else first_ax
            numba__axes.append(ax)

            numba__line_dict[ylabel] = []
            waveform_name = waveform_names[0] # Only the first waveform is plotted in parallel simulation
            for jj in range(number_of_threads):
                # print(jj, name)
                # get line objects
                line, = ax.plot([], [], 
                                linestyle=cjh_plotting.cjh_linestyles[jj], 
                                color=cjh_plotting.cjh_colors[jj], 
                                lw=1.5, 
                                label=waveform_name+str(jj), 
                                alpha=0.7) # zorder
                line.ax = ax
                line.ydata = np.array(0)
                line.MIN = -1e-10
                line.MAX =  1e-10
                numba__line_dict[ylabel].append(line)
            # print()

            ax.set_ylabel(ylabel)
            ax.legend(loc='lower left').set_zorder(202)

        numba__axes[-1].set_xlabel('Time [s]')
        # time_text = first_ax.text(0.02, 0.95, '', transform=first_ax.transAxes)

        mainWindowObject.CONSOLE = CONSOLE
        mainWindowObject.numba__line_dict = numba__line_dict
        mainWindowObject.first_ax = first_ax


    @staticmethod
    def update_line_data(CONSOLE, end_time, line, ydata):
        line.ydata = np.append(line.ydata, ydata)
        line.ydata = line.ydata[-CONSOLE.NUMBER_OF_SAMPLE_TO_SHOW:]
        line.MIN = min(line.ydata) if min(line.ydata) < line.MIN else line.MIN
        line.MAX = max(line.ydata) if max(line.ydata) > line.MAX else line.MAX

        xdata = np.arange(end_time*CONSOLE.SAMPLING_RATE-len(line.ydata), end_time*CONSOLE.SAMPLING_RATE, 1) * CONSOLE.CL_TS # 整数代替浮点数，精度高否则会出现xdata长3001，ydata长3000的问题。
        # print(len(xdata),len(line.ydata))
        line.set_data(xdata, line.ydata)

        # this is slower as we set ylim for each line rather than each axis
        # if line.MIN != line.MAX:
        #     line.ax.set_ylim([line.MIN, line.MAX]) # this .ax is manually assigned when line object is created.
        return line.MIN, line.MAX
