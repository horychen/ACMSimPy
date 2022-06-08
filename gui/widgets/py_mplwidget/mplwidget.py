# ------------------------------------------------------
# -------------------- mplwidget.py --------------------
# ------------------------------------------------------
# from PySide2.QtCore import Qt
# from PySide2.QtWidgets import QWidget, QVBoxLayout

# IMPORT PACKAGES AND MODULES
# ///////////////////////////////////////////////////////////////
import sys

# IMPORT QT CORE
# ///////////////////////////////////////////////////////////////
from qt_core import *

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure

from pylab import mpl
# mpl.rcParams["figure.autolayout"] = True # equivalent to tight_layout=True
# mpl.rcParams['figure.constrained_layout.use'] = True

# 手动调整图片边缘空间
# mpl.rcParams["figure.subplot.hspace"] = 0.25
# mpl.rcParams["figure.subplot.top"]    = 0.995
# mpl.rcParams["figure.subplot.right"]  = 0.995  
# mpl.rcParams["figure.subplot.left"]   = 0.06
# mpl.rcParams["figure.subplot.bottom"] = 0.05 

class MplWidget(QWidget):

    def __init__(self, parent=None):

        QWidget.__init__(self, parent)

        # 获得帆布，颜色设置
        hex_monokai = '#272822'
        hex_wanderson = '#3d444c'
        self.fig = Figure(tight_layout=True, constrained_layout=False, facecolor=hex_wanderson, edgecolor='k', figsize=(8, 16))
        # fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0)
        self.canvas = FigureCanvas(self.fig)
        # self.canvas = FigureCanvas(Figure(tight_layout=True, facecolor='#3d444c', edgecolor='k'))
        # self.canvas = FigureCanvas(Figure(tight_layout=True))

        # print('|||DEBUG:', Qt.StrongFocus, type(Qt.StrongFocus))
        # print('|||DEBUG:', Qt.StrongFocus, type(Qt.StrongFocus))
        # print('|||DEBUG:', Qt.StrongFocus, type(Qt.StrongFocus))

        # 焦点 (https://stackoverflow.com/questions/41366338/matplotlib-navigation-toolbar-shortcuts-not-working)
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        # self.canvas.setFocus() # use your mouse instead of this

        # 布局
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.setLayout(vertical_layout)

        # 绘图区域
        self.canvas.axes = None #self.canvas.figure.add_subplot(111)

        # 工具栏
        self.toolbar = NavigationToolbar(self.canvas, self)

        # add new animation feature
        def on_press(event):
            print('[mplwidget.py] pass to mpl-toolbar:', event.key)
            if event.key.isspace():
                print('Space is pressed and passed to MplWidget.')
                # if anim.running:
                #     anim.event_source.stop()
                # else:
                #     anim.event_source.start()
                # anim.running ^= True
            else:
                # key_press_handler(event, self.canvas, self.mpl_nav_toolbar)
                mpl.backend_bases.key_press_handler(event, self.canvas, self.toolbar)

                # https://matplotlib.org/3.2.2/users/navigation_toolbar.html
                # Navigation Keyboard Shortcuts
                # The following table holds all the default keys, which can be overwritten by use of your matplotlibrc (#keymap.*).
                    # Command Keyboard Shortcut(s)
                    # Home/Reset  h or r or home
                    # Back    c or left arrow or backspace
                    # Forward v or right arrow
                    # Pan/Zoom    p
                    # Zoom-to-rect    o
                    # Save    ctrl + s
                    # Toggle fullscreen   f or ctrl + f
                    # Close plot  ctrl + w
                    # Close all plots shift + w
                    # Constrain pan/zoom to x axis    hold x when panning/zooming with mouse
                    # Constrain pan/zoom to y axis    hold y when panning/zooming with mouse
                    # Preserve aspect ratio   hold CONTROL when panning/zooming with mouse
                    # Toggle major grids  g when mouse is over an axes
                    # Toggle minor grids  G when mouse is over an axes
                    # Toggle x axis scale (log/linear)    L or k when mouse is over an axes
                    # Toggle y axis scale (log/linear)    l when mouse is over an axes
                pass
        self.canvas.mpl_connect('key_press_event', on_press)

    def setMinimumSizeByNumberOfSubplots(self, number_of_subplot, height=200):
        self.setMinimumSize(QSize(500, number_of_subplot*height))
