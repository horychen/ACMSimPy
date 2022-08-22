
# IMPORT PACKAGES AND MODULES
# ///////////////////////////////////////////////////////////////
from gui.core.functions import Functions
from gui.core.emy_functions import EmyFunctions

# IMPORT QT CORE
# ///////////////////////////////////////////////////////////////
from qt_core import *

# IMPORT SETTINGS
# ///////////////////////////////////////////////////////////////
from gui.core.json_settings import Settings

# IMPORT THEME COLORS
# ///////////////////////////////////////////////////////////////
from gui.core.json_themes import Themes

# IMPORT PY ONE DARK WIDGETS
# ///////////////////////////////////////////////////////////////
from gui.widgets import *



class ConsoleWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(ConsoleWindow, self).__init__(*args, **kwargs)

        # self.ui = UI_MainWindow()
        # self.ui.setup_ui(parent=self) # call parent.setCentralWidget(self.central_widget)

        self.setWindowTitle('Console for ACMSimPy')

        self.t = PyLineEdit('CTRL.KA')
        self.l = QLabel("None")
        self.t2 = PyLineEdit('CTRL.cmd_rpm')
        self.l2 = QLabel("None")
        self.t3 = PyLineEdit('CTRL.idq[1]')
        self.l3 = QLabel("None")
        b = QPushButton("DANGER!")
        b.pressed.connect(self.oh_no)
        self.ConsoleWidget_ACMPlot = ConsoleWidget()

        # Controller Commands (User Input)
        themes = Themes()
        self.themes = themes.items
        self.plainTextEdit_ControllerCommands = PyTextEdit(
            text = r'''# Console @ CL_TS
def user_controller_commands(t, ACM, CTRL, reg_id=None, reg_iq=None, reg_speed=None):
    # Sweep frequency
    CTRL.bool_overwrite_speed_commands = True
    return

    # General speed commands
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

    if CTRL.CMD_SPEED_SINE_RPM!=0:
        CTRL.cmd_rpm = CTRL.CMD_SPEED_SINE_RPM * np.sin(2*np.pi*CTRL.CMD_SPEED_SINE_HZ*t)
''',
            # place_holder_text = STRING_WATCH_MAPPING,
            radius = 8,
            border_size = 2,
            color = self.themes["app_color"]["text_foreground"],
            selection_color = self.themes["app_color"]["white"],
            bg_color = self.themes["app_color"]["dark_one"],
            bg_color_active = self.themes["app_color"]["dark_three"],
            context_color = self.themes["app_color"]["context_color"]
        )
        self.plainTextEdit_ControllerCommands.setMinimumHeight(300)


        layout2 = QHBoxLayout()
        layout2.addWidget(self.t)
        layout2.addWidget(self.l)
        layout2.addWidget(self.t2)
        layout2.addWidget(self.l2)
        layout2.addWidget(self.t3)
        layout2.addWidget(self.l3)
        layout2.addWidget(b)

        layout = QVBoxLayout()
        layout.addLayout(layout2)
        layout.addWidget(self.ConsoleWidget_ACMPlot)
        layout.addWidget(self.plainTextEdit_ControllerCommands)

        w = QWidget()
        w.setLayout(layout)
        self.setCentralWidget(w)

        # self.show()

    def oh_no(self):
        print('Console button is clicked.')
