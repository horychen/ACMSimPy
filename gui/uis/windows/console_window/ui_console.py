
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

        w = QWidget()
        w.setLayout(layout)
        self.setCentralWidget(w)

        # self.show()

    def oh_no(self):
        print('Console button is clicked.')
