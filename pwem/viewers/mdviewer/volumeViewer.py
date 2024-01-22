# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os.path
import sys
import numpy as np
import logging

from PIL import Image, ImageOps, ImageFilter
import matplotlib.pyplot as plt

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon, QPixmap, QImage, QKeySequence
from PyQt5.QtWidgets import (QApplication, QMainWindow, QTableWidget, QWidget,
                             QVBoxLayout, QComboBox, QLabel,
                             QSpinBox, QScrollBar, QAction, QMenu, QMenuBar,
                             QFileDialog)
import pyworkflow as pw
from pwem.emlib.image.image_readers import ImageReadersRegistry

logger = logging.getLogger()


class CustomWidget(QWidget):
    """Class to  custom the table cell widget"""
    def __init__(self, data, text='',  autocontrast=False,  gaussianBlurFilter=False):
        super().__init__()
        self._data = data
        self._layout = QVBoxLayout()
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._label = QLabel()
        self._adicinalText = QLabel(text)

        try:
            if autocontrast:
                data = ImageOps.autocontrast(data)
            if gaussianBlurFilter:
                data = data.filter(ImageFilter.GaussianBlur(radius=0.5))
            im = data.convert("RGBA")
            data = im.tobytes("raw", "RGBA")
            qimage = QImage(data, im.size[0], im.size[1],
                                  QImage.Format_RGBA8888)

            pixmap = QPixmap.fromImage(qimage)
            self._label.setPixmap(pixmap)
            self._layout.addSpacing(5)
            self._layout.addWidget(self._label, alignment=Qt.AlignCenter)
            self._layout.addWidget(self._adicinalText, alignment=Qt.AlignCenter)

        except Exception as e:
            logger.error("Error loading the image:", e)

        self.setLayout(self._layout)

    def widgetType(self):
        """Return the type of the widget content"""
        return self._type

    def getWidgetContent(self):
        """Return the widget content"""
        return self._label

    def getData(self):
        """Return the content data"""
        return self._data

    def getValue(self):
        """Return the image path"""
        return self._value

    def getId(self):
        """Return the widget id"""
        return self._id


class CustomScrollBar(QScrollBar):
    """Class to custom the scrollbar widget"""
    def __init__(self):
        super().__init__()
        self.setSingleStep(1)
        self.dragging = False
        self.start_value = 0
        self.start_pos = None

    def wheelEvent(self, event):
        """Handle the mouse wheel event"""
        step = 1
        delta = event.angleDelta().y() / 120  # Number of mouse wheel steps
        # If moved up, decreases the value of scroll
        if delta > 0:
            self.setValue(self.value() - step)
        # If it moves down, it increases the scroll value.
        elif delta < 0:
            self.setValue(self.value() + step)

        event.accept()


class VolumeViewer(QMainWindow):
    """The VolumeViewer show the data as a table"""
    def __init__(self, filePath):
        super().__init__()
        self._filePath = filePath
        self.setWindowTitle("Basic em-file viewer")
        self.setGeometry(100, 100, 800, 600)
        self.setMinimumSize(900, 500)
        self.centralWidget = QWidget(self)
        self.setCentralWidget(self.centralWidget)
        self.applyImageAutocontrast = False
        self.applyImageGaussianBlurFilter = False

        self._createMenu()
        self._createToolBar()
        self._createTableView()
        self._defineStyles()
        self.loadVolume()

    def _createMenu(self):
        """Create a menu bar with "File" options (open and exit)."""
        #  File menu
        menu_bar = QMenuBar(self)
        fileMenu = QMenu("&File", self)
        menu_bar.addMenu(fileMenu)

        # File actions
        self.openAction = QAction(self)
        self.openAction.setText("&Open...")
        self.openAction.setShortcut(QKeySequence("Ctrl+O"))
        # self.openAction.setIcon(QIcon(getImage(FOLDER)))
        self.openAction.triggered.connect(self.openFile)

        self.exitAction = QAction(self)
        self.exitAction.setText("E&xit")
        self.exitAction.setShortcut(QKeySequence("Ctrl+X"))
        # self.exitAction.setIcon(QIcon(getImage(EXIT)))
        self.exitAction.triggered.connect(sys.exit)

        fileMenu.addAction(self.openAction)
        fileMenu.addSeparator()
        fileMenu.addAction(self.exitAction)

        imageMenu = QMenu('&Image', self)
        self.autocontrastAction = QAction('Autocontrast', self)
        self.autocontrastAction.setCheckable(True)
        self.autocontrastAction.setShortcut(QKeySequence("Ctrl+A"))
        self.autocontrastAction.toggled.connect(self.applyAutocontrast)

        self.gaussianBlurAction = QAction('Gaussian blur filter', self)
        self.gaussianBlurAction.setCheckable(True)
        self.gaussianBlurAction.setShortcut(QKeySequence("Ctrl+G"))
        self.gaussianBlurAction.toggled.connect(self.gaussianBlurFilter)

        imageMenu.addAction(self.autocontrastAction)
        imageMenu.addAction(self.gaussianBlurAction)

        menu_bar.addMenu(imageMenu)


        self.setMenuBar(menu_bar)

    def _createToolBar(self):
        """Create the Tool bar """
        toolBar = self.addToolBar("")

        # Adding axis selector
        blockLabelIcon = QLabel('   ')
        toolBar.addWidget(blockLabelIcon)
        axisLabel = QLabel('Select an axis')
        axisLabel.setToolTip('Select an axis')
        iconZoom = pw.findResource('file_vol.png')
        icon = QIcon(iconZoom)
        axisLabel.setPixmap(icon.pixmap(20, 20))
        self.axisSelector = QComboBox(self)
        self.axisSelector.addItems(["Z", "Y", "X"])
        self.axisSelector.setToolTip('Axis')
        self.axisSelector.currentTextChanged.connect(self.showVolume)
        toolBar.addWidget(axisLabel)
        toolBar.addWidget(self.axisSelector)
        blockLabelIcon = QLabel('   ')
        toolBar.addWidget(blockLabelIcon)

        # Adding zoom action
        self.zoomLabel = QLabel('Size in pixels')
        iconZoom = pw.findResource('fa-search.png')
        icon = QIcon(iconZoom)
        self.zoomLabel.setPixmap(icon.pixmap(20, 20))
        self.zoomLabel.setToolTip('Size in pixels')
        self.zoomLabel.setEnabled(False)
        self.zoom = QSpinBox()
        self.zoom.setMaximum(2000)
        self.zoom.setMinimum(50)
        self.zoom.setValue(150)
        self.zoom.setToolTip('Size in pixels')
        self.zoom.setFixedWidth(70)
        self.zoom.setAlignment(Qt.AlignRight)
        self.zoom.setEnabled(True)
        self.zoom.valueChanged.connect(self.showVolume)

        toolBar.addWidget(self.zoomLabel)
        toolBar.addWidget(self.zoom)

    def _createTableView(self):
        """Initializes a table view for displaying slices."""
        # Adding the table
        self.layout = QVBoxLayout()
        self.tableWidget = QTableWidget()

        self.vScrollBar = CustomScrollBar()
        self.tableWidget.setVerticalScrollBar(self.vScrollBar)
        self.vScrollBar.setValue(0)
        self.vScrollBar.valueChanged.connect(self.showVolume)

        self.layout.addWidget(self.tableWidget)
        self.centralWidget.setLayout(self.layout)

        self._triggeredResize = False
        self._isFirstImage = True

    def _defineStyles(self):
        """Define the gallery styles. Create a border to the selected image"""
        self.tableWidget.setStyleSheet(""" QTableWidget::item:selected { border: 1.5px dashed red; }"""
                                       """ QTableWidget::item:focus { border: 1.5px dashed blue;} """)

    def openFile(self):
        """This method uses QFileDialog to open a volume file (.stk or .mrc)."""
        filepath, _ = QFileDialog.getOpenFileNames(self, 'Open volume or stack',
                                                   '', '(*.stk *.mrc)')
        if filepath and os.path.exists(filepath[0]):
            self._filePath = filepath[0]
            self.loadVolume()

    def getZoom(self):
        return self.zoom.value()

    def _calculateVisibleColumns(self):
        """Method that calculate how many columns are visible"""
        viewportWidth = self.tableWidget.parent().parent().width() if self.tableWidget.parent().parent() else self.tableWidget.viewport().width()
        visibleCols = 0
        columnaX = 0
        viewportWidth = viewportWidth - 100  # Ensuring that horizontal scrolling does not appear
        while columnaX < viewportWidth:
            columnaX += self.getZoom() + 5  # Making the cells a little larger than the contents
            visibleCols += 1
        return visibleCols - 1

    def _calculateVisibleRows(self):
        """Method that calculate how many rows are visible"""
        viewportHeight = self.tableWidget.parent().parent().height() if self.tableWidget.parent().parent() else self.tableWidget.viewport().height()
        visibleRows = viewportHeight // self.getZoom() + 1
        return visibleRows

    def loadVolume(self):
        """This method reads the volume data from the selected file and calls showVolume to display it."""
        if self._filePath:
            ext = self._filePath.split('.')[-1]
            self.imageReader = ImageReadersRegistry._readers[ext]
            self.volumeData = self.imageReader.getArray(self._filePath)
            self.showVolume()
        else:
            logger.error("Unable to upload the file %s. Make sure the path is correct." % self._filePath)

    def showVolume(self):
        """This method sets up the table view based on the volume data"""
        self.tableWidget.setColumnCount(0)
        self._columnCount = self._calculateVisibleColumns()
        self._rowsCount = int(self.volumeData.shape[0] / self._columnCount)
        if self.volumeData.shape[0] % self._columnCount > 0:
            self._rowsCount += 1
        selectedAxis = self.axisSelector.currentIndex()
        self.tableWidget.setRowCount(self._rowsCount)
        self.tableWidget.setColumnCount(self._columnCount)
        currentValue = self.vScrollBar.value() - 1 if self.vScrollBar.value() else 0
        visibleRows = self._calculateVisibleRows()
        visibleColumn = self._calculateVisibleColumns()

        self._loadImages(visibleRows, visibleColumn, selectedAxis, currentValue)

    def _loadImages(self, visibleRows, visibleColumn, selectedAxis,
                    currentValue):
        """This method loads and displays the volume slices in the table view"""
        index = currentValue * visibleColumn
        nslices = self.volumeData.shape[0]
        self.iMax = self.volumeData.max()
        self.iMin = self.volumeData.min()

        for row in range(visibleRows):
            for col in range(visibleColumn):
                if index >= nslices:
                    break
                sliceData = np.take(self.volumeData, axis=selectedAxis,
                                    indices=index)
                self.addSlice(sliceData, currentValue + row, col, index)
                self.tableWidget.setColumnWidth(col, self.getZoom() + 5)
                index += 1
            self.tableWidget.setRowHeight(currentValue + row, self.getZoom() + 5)

    def addSlice(self, sliceData, row, col, index):
        """This method converts the volume slice to a PIL image and displays it in
        the table view using a custom widget (CustomWidget)"""

        # Convert the slice narray into PIL image
        im255 = ((sliceData - self.iMin) / (self.iMax - self.iMin) * 255).astype(np.uint8)
        image = Image.fromarray(im255)

        # Create a thumbnail
        size = self.getZoom()
        sizeX, sizeY = image.size
        height = int(size * sizeY / sizeX)
        imageR = image.resize((size, height))
        imageR.thumbnail((size, height))
        widget = CustomWidget(imageR, text='slice %s' % (index+1), autocontrast=self.applyImageAutocontrast,
                              gaussianBlurFilter=self.applyImageGaussianBlurFilter)

        self.tableWidget.setCellWidget(row, col, widget)

    def resizeEvent(self, event):
        """Handle the resize event and trigger a redraw of the volume slices
        when resizing the window"""
        if self._triggeredResize:
            self.showVolume()
        self._triggeredResize = True

    def applyAutocontrast(self):
        """Apply autocontrast to images"""
        if self.autocontrastAction.isChecked():
            self.applyImageAutocontrast = True
        else:
            self.applyImageAutocontrast = False

        self.showVolume()

    def gaussianBlurFilter(self):
        """Apply a gaussian blur filter to the images"""
        if self.gaussianBlurAction.isChecked():
            self.applyImageGaussianBlurFilter = True
        else:
            self.applyImageGaussianBlurFilter = False

        self.showVolume()


def launchViewer(filePath):
    app = QApplication([])
    viewer = VolumeViewer(filePath)
    viewer.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    filePath = sys.argv[-1]
    launchViewer(filePath)
