
# # from app_modules import *
# # from PySide6 import QtCore, QtGui, QtWidgets
# # from PySide6.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect, QSize, QTime, QUrl, Qt, QEvent, Signal, QEasingCurve, QTimer)
# # from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence, QLinearGradient, QPalette, QPainter, QPixmap, QRadialGradient, QPen, QPainterPath, QImage)
# from PySide2.QtWidgets import QApplication, QMainWindow, QPushButton, QSizePolicy, QWidget, QLineEdit
# from PySide2.QtWidgets import QStyleOption, QStyle
# from PySide2 import QtGui, QtCore

    # from main import *
    # if bool_use_PyQt5:
    #     from PyQt5 import QtGui, QtCore
    #     from PyQt5.QtWidgets import QApplication, QWidget, QStyleOption, QStyle
    # else:
    #     from PySide2 import QtGui, QtCore
    #     from PySide2.QtWidgets import QApplication, QWidget, QStyleOption, QStyle


# IMPORT PACKAGES AND MODULES
# ///////////////////////////////////////////////////////////////
import sys

# IMPORT QT CORE
# ///////////////////////////////////////////////////////////////
from qt_core import *

# https://stackoverflow.com/questions/11513132/embedding-ipython-qt-console-in-a-pyqt-application
# from qtconsole.qt import QtGui # This works with qtconsole-4.6.0
from qtconsole.rich_jupyter_widget import RichJupyterWidget
from qtconsole.inprocess import QtInProcessKernelManager

# import logging

class ConsoleWidget(RichJupyterWidget):

    # def __init__(self, customBanner=None, *args, **kwargs):
    def __init__(self, parent=None, *args, **kwargs):

        # print('DeBUG consolewidget:', parent)

        # QWidget.__init__(self, parent)
        r''' I have to comment this out because of following error:
            D:\DrH\Codes\emachineryTestPYPI\emachinery\guiv2>python main.py
            Traceback (most recent call last):
              File "main.py", line 289, in <module>
                window = MainWindow() # ui=QUiLoader().load("GUI_BASE.ui", None)
              File "main.py", line 25, in __init__
                self.ui.setupUi(self)
              File "D:\DrH\Codes\emachineryTestPYPI\emachinery\guiv2\ui_main.py", line 1294, in setupUi
                self.ConsoleWidget = ConsoleWidget(self.page_namePlateData)
              File "D:\DrH\Codes\emachineryTestPYPI\emachinery\guiv2\consolewidget.py", line 16, in __init__
                QWidget.__init__(self, parent)
            RuntimeError: You can't initialize an object twice!
            QBasicTimer::start: QBasicTimer can only be used with threads started with QThread
        '''

        super(ConsoleWidget, self).__init__(*args, **kwargs)

        # useless for guiv2?
        # self.setStyleSheet("background: rgb(27, 29, 35);")

            #         css = '''
            # color: blue;
            # background-color: yellow;
            # selection-color: yellow;
            # selection-background-color: blue;
            # '''
            #         self.setStyleSheet(css)

        # if customBanner is not None:
        #     self.banner = customBanner
        self.banner = "This is an embedded qtconsole.\nSuggested commands:\n\tvars(CONSOLE)\n\n" \
                    + '''You will have access to CONSOLE, ACM, and CTRL after you start python+numba based simulation.
\nSome ideas: 
\tCTRL.CMD_SPEED_SINE_HZ = 20
\treg_speed.OutLimit*=2
\t_dir(CONSOLE)
\tCONSOLE.NUMBER_OF_SAMPLE_TO_SHOW=6000
\tACM.TLoad=2
\nFinally, you can double click on the canvas to stop/start simulation.
\n
'''
                    # [attr for attr in dir(emy) if not callable(getattr(emy, attr)) and not attr.startswith('__')]
                    # https://stackoverflow.com/questions/1398022/looping-over-all-member-variables-of-a-class-in-python

        # 无法修改字体大小
        self.font_size = 6+3
        # self.font_sizeInt = 9+3
        # print(self.font)
        self.font.setPointSize(24) # https://doc.qt.io/qt-5/qfont.html#setPointSize

        self.kernel_manager = kernel_manager = QtInProcessKernelManager()
        kernel_manager.start_kernel(show_banner=False)
        # kernel_manager.kernel.log.setLevel(logging.CRITICAL) #
        kernel_manager.kernel.gui = 'qt'
        self.kernel_client = kernel_client = self._kernel_manager.client()
        kernel_client.start_channels()

        def stop():
            kernel_client.stop_channels()
            kernel_manager.shutdown_kernel()
            guisupport.get_app_qt().exit()
        self.exit_requested.connect(stop)

        # # useless to set background color for promoted widget
        # self.setAttribute(QtCore.Qt.WA_StyledBackground, True) # https://stackoverflow.com/questions/54965088/pyqt-promoted-widget-background-issues?noredirect=1&lq=1

        self.execute_command('''_dir = lambda obj: [method for method in dir(obj) if not method.startswith('__')]''')

    """ Do we need to comment out this function? """
    def paintEvent(self, evt):
        # no effect???
        # Stefan Reinhardt from https://stackoverflow.com/questions/7276330/qt-stylesheet-for-custom-widget
        super(ConsoleWidget, self).paintEvent(evt)
        opt = QStyleOption()
        opt.initFrom(self)
        # opt.init(self)
        p = QPainter(self) # QtGui.QPainter
        s = self.style()
        s.drawPrimitive(QStyle.PE_Widget, opt, p, self) 

    def push_vars(self, variableDict):
        """
        Given a dictionary containing name / value pairs, push those variables
        to the Jupyter console widget
        """
        self.kernel_manager.kernel.shell.push(variableDict)

    def clear(self):
        """
        Clears the terminal
        """
        self._control.clear()

        # self.kernel_manager

    def print_text(self, text):
        """
        Prints some plain text to the console
        """
        self._append_plain_text(text)

    def execute_command(self, command):
        """
        Execute a command in the frame of the console widget
        """
        self._execute(command, False)


if __name__ == '__main__':

    app = QApplication([])
    widget = ConsoleWidget()
    # widget.setStyleSheet("background: rgb(27, 29, 35);")
    widget.show()
    if bool_use_PyQt5:
        app.exec_()
    else:
        sys.exit(app.exec_())
