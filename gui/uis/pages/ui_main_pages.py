# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'main_pages.ui'
##
## Created by: Qt User Interface Compiler version 6.3.0
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QFrame, QHBoxLayout, QLabel,
    QScrollArea, QSizePolicy, QStackedWidget, QVBoxLayout,
    QWidget)

class Ui_MainPages(object):
    def setupUi(self, MainPages):
        if not MainPages.objectName():
            MainPages.setObjectName(u"MainPages")
        MainPages.resize(1052, 837)
        self.horizontalLayout_MainPages = QHBoxLayout(MainPages)
        self.horizontalLayout_MainPages.setSpacing(0)
        self.horizontalLayout_MainPages.setObjectName(u"horizontalLayout_MainPages")
        self.horizontalLayout_MainPages.setContentsMargins(0, 0, 0, 0)
        self.pages = QStackedWidget(MainPages)
        self.pages.setObjectName(u"pages")
        sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pages.sizePolicy().hasHeightForWidth())
        self.pages.setSizePolicy(sizePolicy)
        self.page_1 = QWidget()
        self.page_1.setObjectName(u"page_1")
        sizePolicy1 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.page_1.sizePolicy().hasHeightForWidth())
        self.page_1.setSizePolicy(sizePolicy1)
        self.page_1.setStyleSheet(u"font-size: 14pt")
        self.page_1_layout_2 = QVBoxLayout(self.page_1)
        self.page_1_layout_2.setObjectName(u"page_1_layout_2")
        self.welcome_base = QFrame(self.page_1)
        self.welcome_base.setObjectName(u"welcome_base")
        self.welcome_base.setMinimumSize(QSize(300, 150))
        self.welcome_base.setMaximumSize(QSize(1000, 150))
        self.welcome_base.setFrameShape(QFrame.NoFrame)
        self.welcome_base.setFrameShadow(QFrame.Raised)
        self.logo = QFrame(self.welcome_base)
        self.logo.setObjectName(u"logo")
        self.logo.setGeometry(QRect(580, 0, 300, 120))
        sizePolicy2 = QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.logo.sizePolicy().hasHeightForWidth())
        self.logo.setSizePolicy(sizePolicy2)
        self.logo.setMinimumSize(QSize(300, 120))
        self.logo.setMaximumSize(QSize(300, 120))
        self.logo.setFrameShape(QFrame.NoFrame)
        self.logo.setFrameShadow(QFrame.Raised)
        self.logo_layout = QVBoxLayout(self.logo)
        self.logo_layout.setSpacing(0)
        self.logo_layout.setObjectName(u"logo_layout")
        self.logo_layout.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel(self.welcome_base)
        self.label.setObjectName(u"label")
        self.label.setGeometry(QRect(530, 126, 441, 19))
        sizePolicy1.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy1)
        self.label.setAlignment(Qt.AlignCenter)

        self.page_1_layout_2.addWidget(self.welcome_base)

        self.page_1_layout = QHBoxLayout()
        self.page_1_layout.setObjectName(u"page_1_layout")

        self.page_1_layout_2.addLayout(self.page_1_layout)

        self.pages.addWidget(self.page_1)
        self.page_2 = QWidget()
        self.page_2.setObjectName(u"page_2")
        self.page_2_layout = QVBoxLayout(self.page_2)
        self.page_2_layout.setSpacing(5)
        self.page_2_layout.setObjectName(u"page_2_layout")
        self.page_2_layout.setContentsMargins(5, 5, 5, 5)
        self.scroll_area = QScrollArea(self.page_2)
        self.scroll_area.setObjectName(u"scroll_area")
        self.scroll_area.setStyleSheet(u"background: transparent;")
        self.scroll_area.setFrameShape(QFrame.NoFrame)
        self.scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.scroll_area.setWidgetResizable(True)
        self.contents = QWidget()
        self.contents.setObjectName(u"contents")
        self.contents.setGeometry(QRect(0, 0, 219, 255))
        self.contents.setStyleSheet(u"background: transparent;")
        self.verticalLayout = QVBoxLayout(self.contents)
        self.verticalLayout.setSpacing(15)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.verticalLayout.setContentsMargins(5, 5, 5, 5)
        self.title_label = QLabel(self.contents)
        self.title_label.setObjectName(u"title_label")
        self.title_label.setMaximumSize(QSize(16777215, 40))
        font = QFont()
        font.setPointSize(16)
        self.title_label.setFont(font)
        self.title_label.setStyleSheet(u"font-size: 16pt")
        self.title_label.setAlignment(Qt.AlignCenter)

        self.verticalLayout.addWidget(self.title_label)

        self.description_label = QLabel(self.contents)
        self.description_label.setObjectName(u"description_label")
        self.description_label.setAlignment(Qt.AlignHCenter|Qt.AlignTop)
        self.description_label.setWordWrap(True)

        self.verticalLayout.addWidget(self.description_label)

        self.row_1_layout = QHBoxLayout()
        self.row_1_layout.setObjectName(u"row_1_layout")

        self.verticalLayout.addLayout(self.row_1_layout)

        self.row_2_layout = QHBoxLayout()
        self.row_2_layout.setObjectName(u"row_2_layout")

        self.verticalLayout.addLayout(self.row_2_layout)

        self.row_3_layout = QHBoxLayout()
        self.row_3_layout.setObjectName(u"row_3_layout")

        self.verticalLayout.addLayout(self.row_3_layout)

        self.row_4_layout = QVBoxLayout()
        self.row_4_layout.setObjectName(u"row_4_layout")

        self.verticalLayout.addLayout(self.row_4_layout)

        self.row_5_layout = QVBoxLayout()
        self.row_5_layout.setObjectName(u"row_5_layout")

        self.verticalLayout.addLayout(self.row_5_layout)

        self.scroll_area.setWidget(self.contents)

        self.page_2_layout.addWidget(self.scroll_area)

        self.pages.addWidget(self.page_2)
        self.page_3 = QWidget()
        self.page_3.setObjectName(u"page_3")
        self.page_3.setStyleSheet(u"QFrame {\n"
"	font-size: 16pt;\n"
"}")
        self.page_3_layout = QVBoxLayout(self.page_3)
        self.page_3_layout.setSpacing(0)
        self.page_3_layout.setObjectName(u"page_3_layout")
        self.page_3_layout.setContentsMargins(0, 0, 0, 0)
        self.scrollArea_ACMSimPyScope = QScrollArea(self.page_3)
        self.scrollArea_ACMSimPyScope.setObjectName(u"scrollArea_ACMSimPyScope")
        sizePolicy3 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.scrollArea_ACMSimPyScope.sizePolicy().hasHeightForWidth())
        self.scrollArea_ACMSimPyScope.setSizePolicy(sizePolicy3)
        self.scrollArea_ACMSimPyScope.setMinimumSize(QSize(0, 0))
        self.scrollArea_ACMSimPyScope.setFrameShape(QFrame.Box)
        self.scrollArea_ACMSimPyScope.setFrameShadow(QFrame.Raised)
        self.scrollArea_ACMSimPyScope.setLineWidth(0)
        self.scrollArea_ACMSimPyScope.setWidgetResizable(True)
        self.scrollAreaWidgetContents_ACMSimPyScope = QWidget()
        self.scrollAreaWidgetContents_ACMSimPyScope.setObjectName(u"scrollAreaWidgetContents_ACMSimPyScope")
        self.scrollAreaWidgetContents_ACMSimPyScope.setGeometry(QRect(0, 0, 1052, 837))
        self.scrollAreaWidgetContents_ACMSimPyScope.setAutoFillBackground(True)
        self.scrollAreaWigetContents_verticalLayout = QVBoxLayout(self.scrollAreaWidgetContents_ACMSimPyScope)
        self.scrollAreaWigetContents_verticalLayout.setSpacing(0)
        self.scrollAreaWigetContents_verticalLayout.setObjectName(u"scrollAreaWigetContents_verticalLayout")
        self.scrollAreaWigetContents_verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.scrollArea_ACMSimPyScope.setWidget(self.scrollAreaWidgetContents_ACMSimPyScope)

        self.page_3_layout.addWidget(self.scrollArea_ACMSimPyScope)

        self.pages.addWidget(self.page_3)
        self.page_4 = QWidget()
        self.page_4.setObjectName(u"page_4")
        self.page_4_layout = QVBoxLayout(self.page_4)
        self.page_4_layout.setObjectName(u"page_4_layout")
        self.page_4_layout.setContentsMargins(0, 0, 0, 0)
        self.scrollArea_ACMSimPyScope2 = QScrollArea(self.page_4)
        self.scrollArea_ACMSimPyScope2.setObjectName(u"scrollArea_ACMSimPyScope2")
        sizePolicy3.setHeightForWidth(self.scrollArea_ACMSimPyScope2.sizePolicy().hasHeightForWidth())
        self.scrollArea_ACMSimPyScope2.setSizePolicy(sizePolicy3)
        self.scrollArea_ACMSimPyScope2.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 1050, 835))
        self.scrollAreaWidgetContents.setAutoFillBackground(True)
        self.scrollAreaWidgetContents_verticalLayout_page4 = QVBoxLayout(self.scrollAreaWidgetContents)
        self.scrollAreaWidgetContents_verticalLayout_page4.setObjectName(u"scrollAreaWidgetContents_verticalLayout_page4")
        self.scrollAreaWidgetContents_verticalLayout_page4.setContentsMargins(0, 0, 0, 0)
        self.scrollArea_ACMSimPyScope2.setWidget(self.scrollAreaWidgetContents)

        self.page_4_layout.addWidget(self.scrollArea_ACMSimPyScope2)

        self.pages.addWidget(self.page_4)

        self.horizontalLayout_MainPages.addWidget(self.pages)


        self.retranslateUi(MainPages)

        self.pages.setCurrentIndex(3)


        QMetaObject.connectSlotsByName(MainPages)
    # setupUi

    def retranslateUi(self, MainPages):
        MainPages.setWindowTitle(QCoreApplication.translate("MainPages", u"Form", None))
        self.label.setText(QCoreApplication.translate("MainPages", u"Welcome to PyOneDark GUI tuned for ACMSimPy", None))
        self.title_label.setText(QCoreApplication.translate("MainPages", u"Custom Widgets Page", None))
        self.description_label.setText(QCoreApplication.translate("MainPages", u"Here will be all the custom widgets, they will be added over time on this page.\n"
"I will try to always record a new tutorial when adding a new Widget and updating the project on Patreon before launching on GitHub and GitHub after the public release.", None))
    # retranslateUi

