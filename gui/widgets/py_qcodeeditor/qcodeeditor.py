#!/usr/bin/python3
# QcodeEditor.py by acbetter.
# -*- coding: utf-8 -*-
# from main import *
# from app_modules import *

# if bool_use_PySide6 == True:

#     from PySide6.QtGui import QColor, QPainter, QTextFormat

# elif bool_use_PyQt5:
#     #     from PyQt5 import QtCore, QtGui
#     #     from PyQt5.QtCore import Qt, QRect, QSize
#     #     from PyQt5.QtWidgets import QWidget, QPlainTextEdit, QTextEdit
#     from PyQt5.QtGui import QColor, QPainter, QTextFormat
#     #     from PyQt5.QtWidgets import QStyleOption, QStyle
# else:
#     #     from PySide2 import QtCore, QtGui
#     #     from PySide2.QtCore import Qt, QRect, QSize
#     #     from PySide2.QtWidgets import QWidget, QPlainTextEdit, QTextEdit
#     from PySide2.QtGui import QColor, QPainter, QTextFormat
#     #     from PySide2.QtWidgets import QStyleOption, QStyle


# IMPORT PACKAGES AND MODULES
# ///////////////////////////////////////////////////////////////
import sys

# IMPORT QT CORE
# ///////////////////////////////////////////////////////////////
from qt_core import *


class QLineNumberArea(QWidget):
    def __init__(self, editor):
        super().__init__(editor)
        self.codeEditor = editor

    def sizeHint(self):
        return QSize(self.editor.lineNumberAreaWidth(), 0)

    def paintEvent(self, event):
        self.codeEditor.lineNumberAreaPaintEvent(event)


class QCodeEditor(QPlainTextEdit):
    def __init__(self, parent=None):
        super(QCodeEditor, self).__init__(parent)
        self.lineNumberArea = QLineNumberArea(self)
        self.blockCountChanged.connect(self.updateLineNumberAreaWidth)
        self.updateRequest.connect(self.updateLineNumberArea)
        self.cursorPositionChanged.connect(self.highlightCurrentLine)
        self.updateLineNumberAreaWidth(0)

        # print(self.style())
        # stylesheet
        # self.setBackground(QColor(Qt.yellow).lighter(160))

        # useless for guiv2?
#         css = '''
# color: blue;
# background-color: yellow;
# selection-color: yellow;
# selection-background-color: blue;
# '''
#         self.setStyleSheet(css)

        # self.setStyleSheet("background: rgb(27, 29, 35);")
        # print('DEBUG: QcodeEditor:', self.styleSheet())

        # useless to set background color for promoted widget
        # self.setAttribute(QtCore.Qt.WA_StyledBackground, True) # https://stackoverflow.com/questions/54965088/pyqt-promoted-widget-background-issues?noredirect=1&lq=1

    # def paintEvent(self, evt):
    #     # no effect???
    #     # Stefan Reinhardt from https://stackoverflow.com/questions/7276330/qt-stylesheet-for-custom-widget
    #     super(QCodeEditor, self).paintEvent(evt)
    #     opt = QStyleOption()
    #     opt.initFrom(self)
    #     p = QtGui.QPainter(self)
    #     s = self.style()
    #     s.drawPrimitive(QStyle.PE_Widget, opt, p, self) 

    def lineNumberAreaWidth(self):
        digits = 1
        max_value = max(1, self.blockCount())
        while max_value >= 10:
            max_value /= 10
            digits += 1
        space = 3 + self.fontMetrics().width('9') * digits
        return space

    def updateLineNumberAreaWidth(self, _):
        self.setViewportMargins(self.lineNumberAreaWidth(), 0, 0, 0)

    def updateLineNumberArea(self, rect, dy):
        if dy:
            self.lineNumberArea.scroll(0, dy)
        else:
            self.lineNumberArea.update(0, rect.y(), self.lineNumberArea.width(), rect.height())
        if rect.contains(self.viewport().rect()):
            self.updateLineNumberAreaWidth(0)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        cr = self.contentsRect()
        self.lineNumberArea.setGeometry(QRect(cr.left(), cr.top(), self.lineNumberAreaWidth(), cr.height()))

    def highlightCurrentLine(self):
        extraSelections = []
        if not self.isReadOnly():
            selection = QTextEdit.ExtraSelection()
            # lineColor = QColor(Qt.yellow).lighter(160) 
            # lineColor = QColor(Qt.gray).lighter(160) # TODO: change to different color for dark stylesheet
            # lineColor = QColor('#659B91')
            lineColor = QColor('#615F4E')
            selection.format.setBackground(lineColor)
            selection.format.setProperty(QTextFormat.FullWidthSelection, True)
            selection.cursor = self.textCursor()
            selection.cursor.clearSelection()
            extraSelections.append(selection)
        self.setExtraSelections(extraSelections)

    def lineNumberAreaPaintEvent(self, event):
        painter = QPainter(self.lineNumberArea)

        painter.fillRect(event.rect(), Qt.lightGray)

        block = self.firstVisibleBlock()
        blockNumber = block.blockNumber()
        top = self.blockBoundingGeometry(block).translated(self.contentOffset()).top()
        bottom = top + self.blockBoundingRect(block).height()

        # Just to make sure I use the right font
        height = self.fontMetrics().height()
        while block.isValid() and (top <= event.rect().bottom()):
            if block.isVisible() and (bottom >= event.rect().top()):
                number = str(blockNumber + 1)
                painter.setPen(Qt.black)
                painter.drawText(0, top, self.lineNumberArea.width(), height, Qt.AlignRight, number)

            block = block.next()
            top = bottom
            bottom = top + self.blockBoundingRect(block).height()
            blockNumber += 1


if __name__ == '__main__':
    import sys
    from PyQt5.QtWidgets import QApplication

    app = QApplication(sys.argv)
    codeEditor = QCodeEditor()
    codeEditor.setStyleSheet("background: rgb(27, 29, 35);")
    # codeEditor.setStyleSheet("background: transparent;")
    codeEditor.show()
    sys.exit(app.exec_())
