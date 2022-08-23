# GUI provided by WANDERSON M.PIMENTA (MADE WITH: Qt Designer and PySide6)

# IMPORT PACKAGES AND MODULES
# ///////////////////////////////////////////////////////////////
import enum
from gui.uis.windows.main_window.functions_main_window import *
import os, sys, numpy as np

# IMPORT QT CORE
# ///////////////////////////////////////////////////////////////
from qt_core import *

# IMPORT SETTINGS
# ///////////////////////////////////////////////////////////////
from gui.core.json_settings import Settings

# IMPORT PY ONE DARK WINDOWS
# ///////////////////////////////////////////////////////////////
# MAIN WINDOW
from gui.uis.windows.main_window import *
from gui.uis.windows.console_window.ui_console import ConsoleWindow

# IMPORT PY ONE DARK WIDGETS
# ///////////////////////////////////////////////////////////////
from gui.widgets import *

# ADJUST QT FONT DPI FOR HIGHT SCALE AN 4K MONITOR
# ///////////////////////////////////////////////////////////////
os.environ["QT_FONT_DPI"] = "96"
# os.environ["QT_SCALE_FACTOR"] = "1.25" # IF IS 4K MONITOR ENABLE 'os.environ["QT_SCALE_FACTOR"] = "2"'

# MAIN WINDOW
# ///////////////////////////////////////////////////////////////
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # SETUP MAIN WINDOw
        # Load widgets from "gui\uis\main_window\ui_main.py"
        # ///////////////////////////////////////////////////////////////
        self.ui = UI_MainWindow()
        self.ui.setup_ui(parent=self) # call parent.setCentralWidget(self.central_widget)

        # LOAD SETTINGS
        # ///////////////////////////////////////////////////////////////
        settings = Settings()
        self.settings = settings.items

        # LOAD USER INPUT (IF ANY)
        if os.path.exists('user_input_scope_dict.txt'):
            with open('user_input_scope_dict.txt', 'r') as f:
                self.STRING_SCOPE_DICT = f.read()
        else:
            self.STRING_SCOPE_DICT = r"""numba__scope_dict = OD([
# Y Labels                # Signal Name of Traces
(r'Speed Out Limit [A]',  ( 'reg_speed.OutLimit'                                ,) ),
(r'$q$-axis voltage [V]', ( 'ACM.udq[1]', 'CTRL.cmd_udq[1]'                     ,) ),
(r'$d$-axis voltage [V]', ( 'ACM.udq[0]', 'CTRL.cmd_udq[0]'                     ,) ),
(r'Torque [Nm]',          ( 'ACM.Tem', 'CTRL.Tem'                               ,) ),
(r'Speed [rpm]',          ( 'CTRL.cmd_rpm', 'CTRL.omega_r_mech', 'CTRL.xS[1]'   ,) ),
(r'Speed Error [rpm]',    ( 'CTRL.cmd_rpm - CTRL.omega_r_mech'                  ,) ),
(r'Position [rad]',       ( 'ACM.theta_d', 'CTRL.theta_d', 'CTRL.xS[0]'         ,) ),
(r'Position mech [rad]',  ( 'ACM.theta_d'                                       ,) ),
(r'$q$-axis current [A]', ( 'ACM.iQ', 'CTRL.cmd_idq[1]'                         ,) ),
(r'$d$-axis current [A]', ( 'ACM.iD', 'CTRL.cmd_idq[0]'                         ,) ),
(r'K_{\rm Active} [A]',   ( 'ACM.KA', 'CTRL.KA'                                 ,) ),
(r'Load torque [Nm]',     ( 'ACM.TLoad', 'CTRL.xS[2]'                           ,) ),
(r'CTRL.iD [A]',          ( 'CTRL.cmd_idq[0]', 'CTRL.idq[0]'                    ,) ),
(r'CTRL.iQ [A]',          ( 'CTRL.cmd_idq[1]', 'CTRL.idq[1]'                    ,) ),
(r'CTRL.uab [V]',         ( 'CTRL.uab[0]', 'CTRL.uab[1]'                        ,) ),
(r'S [1]',                ( 'svgen1.S1', 'svgen1.S2', 'svgen1.S3', 'svgen1.S4', 'svgen1.S5', 'svgen1.S6' ,) ),
])"""

        # LOAD USER INPUT for PARALLEL (IF ANY)
        if os.path.exists('user_input_parallel.txt'):
            with open('user_input_parallel.txt', 'r') as f:
                self.STRING_SCOPE_DICT_PARALLEL = f.read()
        else:
            self.STRING_SCOPE_DICT_PARALLEL = r"""numba__scope_dict = OD([
            (r'Speed [rpm]',      ( 'ACM.omega_r_mech',)),
            (r'Position [rad]',   ( 'ACM.theta_d'     ,)),
            (r'$d$-axis current [A]', ( 'ACM.iD'      ,)),
            (r'$q$-axis current [A]', ( 'ACM.iQ'      ,)),
            (r'Torque [Nm]',      ( 'ACM.Tem'         ,)),])"""

        if os.path.exists('user_input_motor_dict.txt'):
            with open('user_input_motor_dict.txt', 'r') as f:
                self.STRING_USER_INPUT_MOTOR_DICT = f.read()
        else:
            self.STRING_USER_INPUT_MOTOR_DICT = r"""d = d_user_input_motor_dict = {
'CL_TS': 1e-4,
'TIME_SLICE': 0.2,
'NUMBER_OF_SLICES': 20,
'VL_EXE_PER_CL_EXE': 5,
'init_npp': 21,
'init_IN': 72/1.414,
'init_R': 0.1222,
'init_Ld': 0.0011,
'init_Lq': 0.00125,
'init_KE': 0.127,
'init_Rreq': 0.0,
'init_Js': 0.203,
'CTRL.bool_apply_speed_closed_loop_control': True,
'CTRL.bool_apply_decoupling_voltages_to_current_regulation': True,
'CTRL.bool_apply_sweeping_frequency_excitation': True,
'CTRL.bool_overwrite_speed_commands': True,
'CTRL.bool_zero_id_control': False,
'MACHINE_SIMULATIONs_PER_SAMPLING_PERIOD': 1,
'DC_BUS_VOLTAGE': 300,
'FOC_delta': 6.5,
'FOC_desired_VLBW_HZ': 60,
'CL_SERIES_KP': None, # 1.61377,
'CL_SERIES_KI': None, # 97.76
'VL_SERIES_KP': None, # 0.479932,
'VL_SERIES_KI': None, # 30.5565
'VL_LIMIT_OVERLOAD_FACTOR': 1,
'user_system_input_code': '''
''',
}"""

        # LOAD 
        self.STRING_CODES = r"""list_extra_execution_code = [
                                                            "reg_speed.Kp = 1.0*0.0380362",
                                                            "reg_speed.Kp = 2.0*0.0380362",
                                                            "reg_speed.Kp = 0.5*0.0380362",
                                                            ]"""

        # SETUP MAIN WINDOW
        # ///////////////////////////////////////////////////////////////
        self.hide_grips = False # Show/Hide resize grips
        SetupMainWindow.setup_gui(self)

        # SHOW MAIN WINDOW
        # ///////////////////////////////////////////////////////////////
        self.show()

        # Get worker/thread ready 
        # ///////////////////////////////////////////////////////////////
        self.threadpool = QThreadPool()
        print("Multithreading with maximum %d threads" % self.threadpool.maxThreadCount())

        # Get timer running
        self.counter = 0
        self.timer = QTimer()
        self.timer.setInterval(100)
        self.timer.timeout.connect(self.recurring_timer)
        self.timer.start()

        # Prepare canvas for showing simulation waveforms
        # EmyFunctions.prepare_canvas(mainWindowObject=self)
        self.console_window = ConsoleWindow()
        self.CONSOLE = None

    # LEFT MENU BTN IS CLICKED
    # Run function when btn is clicked
    # Check funtion by object name / btn_id
    # ///////////////////////////////////////////////////////////////
    def btn_clicked(self):
        # GET BT CLICKED
        btn = SetupMainWindow.setup_btns(self)

        # Remove Selection If Clicked By "btn_close_left_column"
        if btn.objectName() != "btn_settings":
            self.ui.left_menu.deselect_all_tab()

        # Get Title Bar Btn And Reset Active         
        top_settings = MainFunctions.get_title_bar_btn(self, "btn_top_settings")
        top_settings.set_active(False)

        # LEFT MENU
        # ///////////////////////////////////////////////////////////////

        # HOME BTN
        if btn.objectName() == "btn_home":
            # Select Menu
            self.ui.left_menu.select_only_one(btn.objectName())

            # Load Page 1
            MainFunctions.set_page(self, self.ui.load_pages.page_1)

        # WIDGETS BTN
        if btn.objectName() == "btn_widgets":
            # Select Menu
            self.ui.left_menu.select_only_one(btn.objectName())

            # Load Page 2
            MainFunctions.set_page(self, self.ui.load_pages.page_2)

        # LOAD USER PAGE
        if btn.objectName() == "btn_add_user":
            # Select Menu
            self.ui.left_menu.select_only_one(btn.objectName())

            # Load Page 3 
            MainFunctions.set_page(self, self.ui.load_pages.page_3)

        # LOAD INFO PAGE
        if btn.objectName() == "btn_info":
            # Select Menu
            self.ui.left_menu.select_only_one(btn.objectName())

            # Load Page 4
            MainFunctions.set_page(self, self.ui.load_pages.page_4)

        # Open console window
        if btn.objectName() == "btn_console":
            self.console_window.show()

        # Run real time simulation
        if btn.objectName() == "btn_realtime_simulation":
            # press down the button
            self.ui.left_menu.select_only_one('btn_add_user')
            # Load Page 3 
            MainFunctions.set_page(self, self.ui.load_pages.page_3)

            if False:
                # Run Numba based simulation with freezed window
                EmyFunctions.runPyBasedSimulation(self, progress_callback=None)
            else:
                if self.CONSOLE is not None:
                    self.CONSOLE.bool_exit = True # don't click this button more than once
                # Pass the function to execute in a different thread
                worker = Worker(EmyFunctions.runPyBasedSimulation, self) # Any other args, kwargs are passed to the run function
                worker.signals.result.  connect(self.worker_print_output)
                worker.signals.finished.connect(self.worker_thread_complete)
                worker.signals.progress.connect(self.worker_progress_fn)
                # Execute the thread so the GUI does not freeze
                self.threadpool.start(worker)

        # Run real time simulation (parallel)
        if btn.objectName() == "btn_realtime_simulation_parallel":
            # press down the button
            self.ui.left_menu.select_only_one('btn_info')
            # Load Page 4
            MainFunctions.set_page(self, self.ui.load_pages.page_4)

            if self.CONSOLE is not None:
                self.CONSOLE.bool_exit = True # don't click this button more than once
            # Pass the function to execute in a different thread
            the_codes = self.ui.plainTextEdit_Codes.toPlainText()
            if True:
                list_extra_execution_code = eval(the_codes[the_codes.find('['):])
                print('list_extra_execution_code:', list_extra_execution_code)
            else:
                list_extra_execution_code = [
                    # "reg_speed = acmsimpy.The_PI_Regulator(1.0*0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)",
                    # "reg_speed = acmsimpy.The_PI_Regulator(10 *0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)",
                    # "reg_speed = acmsimpy.The_PI_Regulator(0.1*0.0380362, 0.0380362*30.5565*CTRL.VL_TS, 1*1.414*ACM.IN)",
                    "reg_speed.Kp = 1.0*0.0380362",
                    "reg_speed.Kp = 2.0*0.0380362",
                    "reg_speed.Kp = 0.5*0.0380362",
                ]
            list_worker = []
            for thread_index, extra_execution_code in enumerate(list_extra_execution_code):
                worker = Worker(EmyFunctions.runPyBasedSimulationParallel, self, thread_index, extra_execution_code) # Any other args, kwargs are passed to the run function
                if thread_index == 0:
                    worker.signals.result.  connect(self.worker_print_output)
                    worker.signals.finished.connect(self.worker_thread_complete)
                    worker.signals.progress.connect(self.worker_progress_fn)
                list_worker.append(worker)

            EmyFunctions.prepare_canvas_on_page_4(mainWindowObject=self, number_of_threads=len(list_extra_execution_code))

            # Execute the thread so the GUI does not freeze
            for worker in list_worker:
                self.threadpool.start(worker)

        # BOTTOM INFORMATION
        if btn.objectName() == "btn_info":
            # CHECK IF LEFT COLUMN IS VISIBLE
            if not MainFunctions.left_column_is_visible(self):
                self.ui.left_menu.select_only_one_tab(btn.objectName())

                # Show / Hide
                MainFunctions.toggle_left_column(self)
                self.ui.left_menu.select_only_one_tab(btn.objectName())
            else:
                if btn.objectName() == "btn_close_left_column":
                    self.ui.left_menu.deselect_all_tab()
                    # Show / Hide
                    MainFunctions.toggle_left_column(self)

                self.ui.left_menu.select_only_one_tab(btn.objectName())

            # Change Left Column Menu
            if btn.objectName() != "btn_close_left_column":
                MainFunctions.set_left_column_menu(
                    self, 
                    menu = self.ui.left_column.menus.menu_2,
                    title = "Info tab",
                    icon_path = Functions.set_svg_icon("icon_info.svg")
                )

        # SETTINGS LEFT
        if btn.objectName() == "btn_settings" or btn.objectName() == "btn_close_left_column":
            # CHECK IF LEFT COLUMN IS VISIBLE
            if not MainFunctions.left_column_is_visible(self):
                # Show / Hide
                MainFunctions.toggle_left_column(self)
                self.ui.left_menu.select_only_one_tab(btn.objectName())
            else:
                if btn.objectName() == "btn_close_left_column":
                    self.ui.left_menu.deselect_all_tab()
                    # Show / Hide
                    MainFunctions.toggle_left_column(self)
                self.ui.left_menu.select_only_one_tab(btn.objectName())

            # Change Left Column Menu
            if btn.objectName() != "btn_close_left_column":
                MainFunctions.set_left_column_menu(
                    self, 
                    menu = self.ui.left_column.menus.menu_1,
                    title = "Settings Left Column",
                    icon_path = Functions.set_svg_icon("icon_settings.svg")
                )

        # TITLE BAR MENU
        # ///////////////////////////////////////////////////////////////

        # SETTINGS TITLE BAR
        if btn.objectName() == "btn_top_settings":
            # Toogle Active
            if not MainFunctions.right_column_is_visible(self):
                btn.set_active(True)

                # Show / Hide
                MainFunctions.toggle_right_column(self)
            else:
                btn.set_active(False)

                # Show / Hide
                MainFunctions.toggle_right_column(self)

            # Get Left Menu Btn
            top_settings = MainFunctions.get_left_menu_btn(self, "btn_settings")
            top_settings.set_active_tab(False)

        # DEBUG
        print(f"Button {btn.objectName()}, clicked!")

    # LEFT MENU BTN IS RELEASED
    # Run function when btn is released
    # Check funtion by object name / btn_id
    # ///////////////////////////////////////////////////////////////
    def btn_released(self):
        # GET BT CLICKED
        btn = SetupMainWindow.setup_btns(self)

        # DEBUG
        print(f"Button {btn.objectName()}, released!")

    # RESIZE EVENT
    # ///////////////////////////////////////////////////////////////
    def resizeEvent(self, event):
        SetupMainWindow.resize_grips(self)

    # CLOSE EVENT
    def closeEvent(self, event):
        print("User has clicked the red x on the main window")
        self.console_window.close()
        event.accept()

    # MOUSE CLICK EVENTS
    # ///////////////////////////////////////////////////////////////
    def mousePressEvent(self, event):
        # SET DRAG POS WINDOW

        # self.dragPos = event.globalPos() # this causes DeprecationWarning

        p = event.globalPosition()
        self.dragPos = globalPos = p.toPoint()

    # Timer EVENTS
    # ///////////////////////////////////////////////////////////////
    def recurring_timer(self):
        self.counter += 1
        self.ui.title_label_counter.setText("Counter: %d" % self.counter)

    # 把变量推到 qtconsole 中去
    def console_push_variable(self, d):
        self.console_window.ConsoleWidget_ACMPlot.push_vars(d)

    # Worker (thread) Callback Functions
    # ///////////////////////////////////////////////////////////////
    def worker_print_output(self, s):
        print('Thread print out results:', s, '(EmyFunctions)')
    def worker_thread_complete(self):
        print("THREAD COMPLETE! (EmyFunctions)")
    def worker_progress_fn(self, t0):
        self.ui.title_label_time.setText(f'Elapsed: {t0:.1f} s') # note the data type of t0 is determined in WorkerSignals class definition

        if self.console_window is not None:
            temp = eval('self.'+self.console_window.t .text()); self.console_window.l .setText(f'{temp}')
            temp = eval('self.'+self.console_window.t2.text()); self.console_window.l2.setText(f'{temp}')
            temp = eval('self.'+self.console_window.t3.text()); self.console_window.l3.setText(f'{temp}')

        # self.console_window.l.setText(f'{self.CTRL.KA}')


# SETTINGS WHEN TO START
# Set the initial class and also additional parameters of the "QApplication" class
# ///////////////////////////////////////////////////////////////
if __name__ == "__main__":
    # APPLICATION
    # ///////////////////////////////////////////////////////////////
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon("acm.ico"))
    window = MainWindow()

    # EXEC APP
    # ///////////////////////////////////////////////////////////////
    if not app.exec():
        if window.CONSOLE is not None:
            window.CONSOLE.bool_exit = True # stop the worker/thread that executes EmyFunctions.runPyBasedSimulation
        print('Exit!')
        sys.exit() # The .exec() method in Qt starts the event loop of your QApplication or dialog boxes. 
