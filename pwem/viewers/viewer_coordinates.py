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
import os
import numpy as np
import tkinter as tk
import tkcolorpicker
from tkinter import Menu, ttk
from tkinter import filedialog, messagebox

from PIL import ImageOps, ImageFilter, ImageEnhance, ImageTk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import logging

from matplotlib.widgets import RangeSlider

from pwem.emlib.image.image_readers import ImageReadersRegistry
from pwem.objects import EMSet, SetOfCoordinates, SetOfMicrographs
from pyworkflow.gui import getImage, ToolTip
from pyworkflow.utils import Icon
from pyworkflow.viewer import View, Viewer

logger = logging.getLogger()


class MainWindow:
    """Coordinate viewer"""
    def __init__(self, root, setOfCoordinate, protocol):
        # Initialize the main window
        self.root = root
        self.setOfCoordinate = setOfCoordinate
        self.protocol = protocol
        self.root.title("Scipion coordinates viewer")
        self.root.resizable(False, False)
        self._initGui()

    def _initGui(self):
        self.micId = None
        self.boxSize = self.setOfCoordinate.getBoxSize() if self.setOfCoordinate.getBoxSize() else 100
        self.selectedColor = '#00FF00'
        self.circleButtonRelieve = tk.SUNKEN
        self.squareButtonRelieve = tk.GROOVE
        self.zoomButtonRelieve = tk.SUNKEN
        self.dragButtonRelieve = tk.SUNKEN
        self.pointButtonRelieve = tk.GROOVE
        self.mousePress = False
        self.eraser = False
        self.picking = False
        self.filament = False
        self.zoom = True
        self.drag = True
        self.selectedCoordinate = None
        self.coordinatesDict = dict()
        self.oldCoordinatesDict = dict()
        self.micrographPathDict = dict()
        self.deletedCoordinates = dict()
        self.movedCoordinates = dict()
        self.newCoordinates = dict()
        self.imageCanvasSize = 780
        self.scale = 5
        self.drawSquares = False
        self.drawCircles = True
        self.zoomFactor = 1
        self.moveShape = False
        self.hasChanges = {}
        self.dragData = {'x': 0, 'y': 0, 'item': None}
        self.xOffset = 0
        self.yOffset = 0
        self.particlesWindowVisible = False

        # Menu setup
        menubar = Menu(self.root)
        self.root.config(menu=menubar)

        # File menu
        fileMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=fileMenu)

        fileMenu.add_command(label="Open", command=self.openFile,
                             image=getImage(Icon.ACTION_BROWSE), compound=tk.LEFT, state='disabled')
        fileMenu.add_separator()
        fileMenu.add_command(label="Exit", command=self.root.destroy,
                             image=getImage(Icon.ACTION_OUT), compound=tk.LEFT)

        # Filters menu
        self.filterMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Filters", menu=self.filterMenu)
        self.enhanceContrastVar = tk.BooleanVar(value=True)
        self.gaussianBlurVar = tk.BooleanVar()
        self.bandpassFilterVar = tk.BooleanVar()
        self.invertContrastVar = tk.BooleanVar()
        self.filterMenu.add_command(label="Enhance contrast", command=self.applyContrast, image=getImage(Icon.CHECKED),
                                    compound=tk.LEFT)
        self.filterMenu.add_command(label="Gaussian blur",  command=self.applyGaussianBlur, image=getImage(Icon.UNCHECKED),
                                    compound=tk.LEFT)
        self.filterMenu.add_command(label="Invert contrast", command=self.applyInvertContrast,
                                    image=getImage(Icon.UNCHECKED),
                                    compound=tk.LEFT)
        # Tools menu
        self.toolsMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Tools", menu=self.toolsMenu)
        self.toolsMenu.add_command(label="Power histogram", command=self.applyPowerHistogram,
                                   image=getImage(Icon.ACTION_STATS), compound=tk.LEFT)
        self.toolsMenu.add_command(label="Reset micrograph", command=self.resetMicrograph,
                                   image=getImage(Icon.BROOM),compound=tk.LEFT)
        self.toolsMenu.add_command(label="Restore micrograph", command=self.restoreMicrograph,
                                   image=getImage(Icon.BACKWARD), compound=tk.LEFT)

        # Windows menu
        self.windowsMenu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Window", menu=self.windowsMenu)
        self.windowsMenu.add_command(label="Particles", command=self.showParticles,
                                     image=getImage(Icon.ACTION_GRID), compound=tk.LEFT)

        # Main frame setup
        self.mainFrame = ttk.Frame(self.root)
        self.mainFrame.grid(row=0, column=0, sticky="news")
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        # Toolbar setup
        self._createToolbar(self.mainFrame)

        # Frame divided in two
        self.contentFrame = ttk.Frame(self.mainFrame)
        self.contentFrame.grid(row=2, column=0, sticky="news")

        # Table on the left
        self._createTable()

        # Micrograph on the right
        self._createImageFrame()

        self.onTableClick(None)

    def _createToolbar(self, parent):
        """Create the toolbar with a slider for adjusting the box size"""
        self.toolbar = ttk.Frame(parent, relief=tk.SUNKEN)
        self.toolbar.bind("<Enter>", self.recoveryPointer)

        sizeLabel = ttk.Label(self.toolbar, text="Box size(px):", image=getImage(Icon.FILE_VOL), compound=tk.LEFT)
        sizeLabel.grid(row=0, column=0, padx=10, pady=5, sticky="e")

        self.sizeVar = tk.IntVar(parent, self.boxSize)
        self.sizeSlider = ttk.Scale(self.toolbar, from_=0, to=self.boxSize*4, orient=tk.HORIZONTAL, length=320,
                                    variable=self.sizeVar,  command=self.updateSize)
        self.sizeSlider.grid(row=0, column=1, padx=5, pady=5)
        self.sizeSlider.bind("<B1-Motion>", self.onBoxSizeSlider)  # Bind event for box size slider motion
        self.sizeSlider.bind("<ButtonRelease-1>", self.onSliderRelease)

        self.sizeValueLabel = ttk.Label(self.toolbar, text=f"({int(self.sizeVar.get())})")
        self.sizeValueLabel.grid(row=0, column=2, padx=5, pady=5, sticky="w")

        self.totalMicrographButton = tk.Button(self.toolbar, bg='#CFE8CF', relief=tk.SUNKEN,
                                               activebackground='#CFE8CF', compound=tk.LEFT)
        self.totalMicrographButton.grid(row=0, column=3, padx=20, pady=5, sticky="e")

        self.totalCoordinates = self.setOfCoordinate.getSize()
        self.totalPickButton = tk.Button(self.toolbar, text=f"Total picks: {self.totalCoordinates}",
                                         bg='#C3D7DF', relief=tk.SUNKEN,
                                         activebackground='#C3D7DF', compound=tk.LEFT)

        self.totalPickButton.grid(row=0, column=4, padx=0, pady=5, sticky="e")

        self.toolbar2 = ttk.Frame(parent, relief=tk.SUNKEN)
        self.toolbar2.bind("<Enter>", self.recoveryPointer)
        # Shapes
        shapeLabel = ttk.Label(self.toolbar2, text="Shape:")
        shapeLabel.grid(row=0, column=0, padx=12, pady=5, sticky="e")
        # Circle shape  text='\u20DD',
        self.circleShape = tk.Button(self.toolbar2, command=self.showCircles, width=25, height=25,
                                     image=getImage(Icon.ACTION_CIRCLE),
                                     relief=self.circleButtonRelieve)
        self.circleShape.grid(row=0, column=1, padx=0, pady=5, sticky="w")

        # Square shape text='\u25a1',
        self.squareShape = tk.Button(self.toolbar2,  command=self.showSquare, width=25, height=25,
                                     image=getImage(Icon.UNCHECKED),
                                     relief=self.squareButtonRelieve)
        self.squareShape.grid(row=0, column=2,  padx=0, pady=5, sticky="w")

        # Point shape text='\u2715',
        # self.pointShape = tk.Button(toolbar2,  command=self.showPoint, width=25, height=25,
        #                              image=getImage(Icon.ACTION_CLOSE),
        #                              relief=self.pointButtonRelieve)
        # self.pointShape.grid(row=0, column=3, padx=0, pady=5, sticky="w")

        separador = ttk.Separator(self.toolbar2, orient='vertical')
        separador.grid(row=0, column=4, padx=10, sticky='ns')

        # Palette color
        colorLabel = ttk.Label(self.toolbar2, text="Color:")
        colorLabel.grid(row=0, column=5, padx=5, pady=5, sticky="e")
        self.paletteButton = tk.Button(self.toolbar2,  command=self.showPalette, bg=self.selectedColor,
                                       relief=tk.SUNKEN, activebackground=self.selectedColor)
        self.paletteButton.grid(row=0, column=6, padx=5, pady=5, sticky="e")

        separador = ttk.Separator(self.toolbar2, orient='vertical')
        separador.grid(row=0, column=7, padx=10, sticky='ns')

        # Picker Tools
        pickerLabel = ttk.Label(self.toolbar2, text="Picker tools:")
        pickerLabel.grid(row=0, column=8, padx=5, pady=5, sticky="e")

        self.pickerAction = tk.Button(self.toolbar2, command=self.pickingActivate, width=25, height=25,
                                      image=getImage(Icon.ACTION_PICKING),
                                      relief=self.squareButtonRelieve)
        self.pickerAction.grid(row=0, column=9, pady=5, sticky="e")
        tooltip = "Particle picker"
        ToolTip(self.pickerAction, tooltip, delay=150)

        self.filamentPickerAction = tk.Button(self.toolbar2, command=self.filamentActivate, width=25, height=25,
                                              image=getImage(Icon.ACTION_FILAMENT_PICKING),
                                              relief=self.squareButtonRelieve)
        self.filamentPickerAction.grid(row=0, column=10,  pady=5, sticky="e")
        tooltip = "Filament picker"
        ToolTip(self.filamentPickerAction, tooltip, delay=150)

        self.eraserAction = tk.Button(self.toolbar2, command=self.eraserActivate, width=25, height=25,
                                      image=getImage(Icon.BROOM),
                                      relief=self.squareButtonRelieve)
        self.eraserAction.grid(row=0, column=11, pady=5, sticky="e")
        tooltip = "Eraser"
        ToolTip(self.eraserAction, tooltip, delay=150)

        separador = ttk.Separator(self.toolbar2, orient='vertical')
        separador.grid(row=0, column=12, padx=10, sticky='ns')

        # Image tools

        # Picker Tools
        imageToolsPanel = ttk.Frame(self.toolbar2)
        imageToolsPanel.grid(row=0, column=13, padx=5, pady=0, sticky="e")

        pickerLabel = ttk.Label(imageToolsPanel, text="Navigate:")
        pickerLabel.grid(row=0, column=0, padx=5, pady=5, sticky="e")

        self.fitDisplay = tk.Button(imageToolsPanel, width=25, height=25, command=self.onFitActivate,
                                    image=getImage(Icon.ACTION_EXPAND),
                                    relief=self.squareButtonRelieve)
        self.fitDisplay.grid(row=0, column=1, padx=0, pady=0, sticky="ns")
        tooltip = "Fit to display area"
        ToolTip(self.fitDisplay, tooltip, delay=150)

        self.zoomOnScroll = tk.Button(imageToolsPanel, width=25, height=25, command=self.onZoomActivate,
                                    image=getImage(Icon.ACTION_ZOOM),
                                    relief=self.zoomButtonRelieve)
        self.zoomOnScroll.grid(row=0, column=2, padx=0, pady=0, sticky="ns")
        tooltip = "Zoom on scroll"
        ToolTip(self.zoomOnScroll, tooltip, delay=150)

        self.dragButton = tk.Button(imageToolsPanel, width=25, height=25, command=self.onDragActivate,
                                      image=getImage(Icon.ACTION_HAND),
                                      relief=self.dragButtonRelieve)
        self.dragButton.grid(row=0, column=15, padx=0, pady=0, sticky="ns")
        tooltip = "Click and drag to move"
        ToolTip(self.dragButton, tooltip, delay=150)

        # Label to display coordinates and pixel value
        self.infoLabel = ttk.Label(imageToolsPanel, text="")
        self.infoLabel.grid(row=0, column=17, sticky="w", padx=250)

        self.toolbar.grid(row=0, column=0, sticky="ew")
        self.toolbar2.grid(row=1, column=0, sticky="ew")

    def recoveryPointer(self, event):
        """Restoring the default cursor"""
        self.root.config(cursor='')
        self.infoLabel.config(text='')

    def pickingActivate(self):
        """ Recovering some actions when picking button is pressed """
        if self.picking:
            self.pickerAction.configure(relief=tk.GROOVE)
            self.picking = False
        else:
            self.pickerAction.configure(relief=tk.SUNKEN)
            self.eraserAction.configure(relief=tk.GROOVE)
            self.dragButton.configure(relief=tk.GROOVE)
            self.filamentPickerAction.configure(relief=tk.GROOVE)
            self.root.config(cursor='')
            self.picking = True
            self.eraser = False
            self.filament = False
            self.drag = True

    def eraserActivate(self):
        """ Recovering some actions when eraser button is pressed """
        if self.eraser:
            self.eraserAction.configure(relief=tk.GROOVE)
            self.eraser = False
        else:
            self.pickerAction.configure(relief=tk.GROOVE)
            self.filamentPickerAction.configure(relief=tk.GROOVE)
            self.dragButton.configure(relief=tk.GROOVE)
            self.eraserAction.configure(relief=tk.SUNKEN)
            self.picking = False
            self.filament = False
            self.eraser = True
            self.drag = True

    def filamentActivate(self):
        """ Recovering some actions when filament button is pressed """
        if self.filament:
            self.filamentPickerAction.configure(relief=tk.GROOVE)
            self.filament = False
        else:
            self.pickerAction.configure(relief=tk.GROOVE)
            self.eraserAction.configure(relief=tk.GROOVE)
            self.dragButton.configure(relief=tk.GROOVE)
            self.filamentPickerAction.configure(relief=tk.SUNKEN)
            self.root.config(cursor='')
            self.filament = True
            self.picking = False
            self.eraser = False
            self.drag = False

    def onBoxSizeSlider(self, event):
        """Handle box size slider motion"""
        self.boxSize = int(self.sizeVar.get())
        self.resizeRadious()

    def resizeRadious(self):
        """Calculate the shape radius taking into account the box size value"""
        self.shapeRadius = self.boxSize / self.scale / 2 * self.zoomFactor
        self.auxCoordinatesDict = dict()
        for index, coord in self.shapes.items():
            self.imageCanvas.coords(index, [coord[0] + self.xOffset - self.shapeRadius,
                                            coord[1] + self.yOffset - self.shapeRadius,
                                            coord[0] + self.xOffset + self.shapeRadius,
                                            coord[1] + self.yOffset + self.shapeRadius])

    def onSliderRelease(self, event):
        """Handle box size slider release"""
        if self.particlesWindowVisible:
            self.extractImages()

    def showCircles(self):
        """Draw a circle on the particle"""
        if self.circleButtonRelieve == tk.GROOVE:
            self.circleButtonRelieve = tk.SUNKEN
            self.drawCircles = True
        else:
            self.circleButtonRelieve = tk.GROOVE
            self.drawCircles = False

        self.circleShape.configure(relief=self.circleButtonRelieve)
        self.showCoordinates()
        if self.particlesWindowVisible:
            self.extractImages()

    def showSquare(self):
        """Draw a square on the particle"""
        if self.squareButtonRelieve == tk.GROOVE:
            self.squareButtonRelieve = tk.SUNKEN
            self.drawSquares = True
        else:
            self.squareButtonRelieve = tk.GROOVE
            self.drawSquares = False

        self.squareShape.configure(relief=self.squareButtonRelieve)
        self.showCoordinates()
        if self.particlesWindowVisible:
            self.extractImages()

    def showPoint(self):
        """Draw a point on the particle"""
        if self.pointButtonRelieve == tk.GROOVE:
            self.pointButtonRelieve = tk.SUNKEN
        else:
            self.pointButtonRelieve = tk.GROOVE

        self.pointShape.configure(relief=self.pointButtonRelieve)
        self.showCoordinates()
        if self.particlesWindowVisible:
            self.extractImages()

    def showPalette(self):
        """Color picker dialog"""
        color = tkcolorpicker.ColorPicker(self.root)
        self.root.wait_window(color)
        if color.get_color():
            self.selectedColor = color.get_color()[2]
            self.paletteButton.configure(bg=self.selectedColor, activebackground=self.selectedColor)
            self.showCoordinates()
            if self.particlesWindowVisible:
                self.extractImages()

    def _createTable(self):
        """Create a table on the left side to display information"""
        columns = ('Index', 'File', 'Particles', 'Updated')
        data = {}

        setOfMicrographs = self.setOfCoordinate.getMicrographs()
        for micrograph in setOfMicrographs.iterItems():
            micId = str(micrograph.getObjId())
            self.coordinatesDict[micId] = {}
            self.oldCoordinatesDict[micId] = {}
            self.deletedCoordinates[micId] = {}
            self.movedCoordinates[micId] = {}
            self.newCoordinates[micId] = {}
            self.micrographPathDict[micId] = (micrograph.getFileName(), micrograph.getObjId())
            data[micId] = (micrograph.getObjId(), micrograph.getMicName(), 0, 'No')

        label = '_micName' if isinstance(setOfMicrographs, SetOfMicrographs) else '_micId'
        for micAgg in self.setOfCoordinate.aggregate(["count"], label, [label, "_micId"]):
            micName = micAgg[label]
            data[str(micAgg['_micId'])] = (micAgg['_micId'], micName, micAgg['count'], 'No')

        self.table = ttk.Treeview(self.contentFrame, columns=columns, show="headings")
        self.table.bind("<Enter>", self.recoveryPointer)
        self.table.bind("<ButtonRelease-1>", self.onTableClick)  # Bind event for table click
        for col, width in zip(columns, (50, 300, 100, 100)):
            self.table.heading(col, text=col, anchor="center")
            self.table.column(col, anchor="center", width=width)

        for row in data.values():
            self.table.insert("", tk.END, values=row)

        y_scrollbar = ttk.Scrollbar(self.contentFrame, orient="vertical", command=self.table.yview)
        self.table.configure(yscrollcommand=y_scrollbar.set)

        self.table.grid(row=0, column=0, sticky='nsew')
        y_scrollbar.grid(row=0, column=1, sticky="ns")

        firstItem = self.table.get_children()[0]
        self.table.selection_set(firstItem)  # Select the first item by default
        self.contentFrame.grid_rowconfigure(0, weight=1)
        self.contentFrame.grid_columnconfigure(0, weight=1)

        self.totalMicrographButton.configure(text=f"Total micrograph: {len(data)}")

    def onTableClick(self, event):
        """Handle the event when a row in the table is clicked to select a micrograph"""
        selected_item = self.table.selection()
        if selected_item:
            micrograph = self.table.item(selected_item, "values")
            if len(micrograph) > 1 and micrograph[0] != self.micId:
                self.micId = micrograph[0]
                self.zoomFactor = 1
                self.xOffset = 0
                self.yOffset = 0
                self.dragData = {'x': 0, 'y': 0, 'item': None}
                self.selectedCoordinate = None
                self.loadCoordinates(micrograph)
                self.showCoordinates()
                if self.particlesWindowVisible:
                    self.extractImages()

    def loadCoordinates(self, micrograph):
        """Load the coordinate for a given micrograph """
        micId = micrograph[0]
        if not self.coordinatesDict[self.micId] and micrograph[2] != 0:
            for index, coordinate in enumerate(self.setOfCoordinate.iterCoordinates(int(micId))):
                self.coordinatesDict[self.micId][index+2] = (coordinate.getX(), coordinate.getY(), coordinate.getObjId())
                self.oldCoordinatesDict[self.micId][index+2] = (coordinate.getX(), coordinate.getY(), coordinate.getObjId())

    def _createImageFrame(self):
        """Create the image frame on the right side to display micrographs and associated coordinates"""
        self.imageFrame = ttk.Frame(self.contentFrame)
        self.imageFrame.grid(row=0, column=2, sticky="nw")
        self.contentFrame.grid_rowconfigure(0, weight=1)
        self.contentFrame.grid_columnconfigure(1, weight=1)

        self.imageCanvas = tk.Canvas(self.imageFrame, width=self.imageCanvasSize, height=self.imageCanvasSize)  #, borderwidth=5, highlightthickness=1, highlightbackground='red')
        self.imageCanvas.grid(row=0, column=0, sticky="nw", pady=10)
        self.imageCanvas.bind("<Motion>", self.onMotion)
        self.imageCanvas.bind("<B1-Motion>", self.onMotion)
        self.imageCanvas.bind("<Button-1>", self.onClickPress)
        self.imageCanvas.bind('<ButtonRelease-1>', self.onClickRelease)
        self.imageCanvas.bind("<Button-4>", self.zoomerP)
        self.imageCanvas.bind("<Button-5>", self.zoomerM)

        buttonsFrame = tk.Frame(self.imageFrame)
        buttonsFrame.grid(row=1, column=0, sticky="ns", padx=510)
        buttonsFrame.bind("<Enter>", self.recoveryPointer)

        # Button to close de application
        closeButton = tk.Button(buttonsFrame, text='Close', command=self.root.destroy, font=('Helvetica', 10, 'bold'),
                                      image=getImage(Icon.BUTTON_CANCEL), compound=tk.LEFT)
        closeButton.grid(row=0, column=0, sticky="w", padx=5)

        # Button to create a new set of Coordinates
        coordinateButton = tk.Button(buttonsFrame, text='Coordinate', fg='white', font=('Helvetica', 10, 'bold'),
                                     image=getImage(Icon.ACTION_NEW), compound=tk.LEFT, bg='#B22A2A',
                                     activebackground='#B22A2A', command=self._createOutput)
        coordinateButton.grid(row=0, column=1, sticky="w", padx=5)

        self.imageFrame.grid_rowconfigure(0, weight=1)
        self.imageFrame.grid_columnconfigure(0, weight=1)

    def _createOutput(self):
        """Create a new set of coordinates"""
        result = messagebox.askquestion("Confirmation", "Are you going to generate a new set of coordinates with the new changes?",
                                        icon='warning', **{'parent': self.root})
        if result == messagebox.YES:
            micSet = self.setOfCoordinate.getMicrographs()
            coordSet = self.protocol._createSetOfCoordinates(micSet, suffix=str(self.protocol.getOutputsSize()))
            coordSet.copyInfo(self.setOfCoordinate)
            coordSet.setBoxSize(self.boxSize)
            coordSet.copyItems(self.setOfCoordinate, updateItemCallback=self._removeCoordinate)

            for micId, coord in self.newCoordinates.items():
                for newCoord in self.newCoordinates[micId].values():
                    newCoordinate = self.setOfCoordinate.getFirstItem().clone()
                    newCoordinate.setObjId(None)
                    micrographs = self.setOfCoordinate.getMicrographs()
                    newCoordinate.setMicrograph(micrographs[self.micrographPathDict[micId][1]])
                    newCoordinate.setPosition(newCoord[0], newCoord[1])
                    coordSet.append(newCoordinate)

            coordSet.write()
            self.protocol._defineOutputs(**{'coordinates_' + str(self.protocol.getOutputsSize()+1): coordSet})
            self.protocol._defineSourceRelation(micSet, coordSet)

    def _removeCoordinate(self, item, row):
        """Remove the deleted coordinate"""
        for micId, coord in self.newCoordinates.items():
            if item.getObjId() in self.deletedCoordinates[micId]:
                setattr(item, "_appendItem", False)
            if item.getObjId() in self.movedCoordinates[micId]:
                item.setX(self.movedCoordinates[micId][item.getObjId()][0])
                item.setY(self.movedCoordinates[micId][item.getObjId()][1])

    def onFitActivate(self):
        """Fit to display area and reset all parameters"""
        self.zoomFactor = 1
        self.showCoordinates()
        if self.particlesWindowVisible:
            self.extractImages()
        self.dragData = {'x': 0, 'y': 0, 'item': None}
        self.xOffset = 0
        self.yOffset = 0

    def onZoomActivate(self):
        """Allow zoom on scroll"""
        if self.zoom:
            self.zoomOnScroll.configure(relief=tk.GROOVE)
            self.zoom = False
        else:
            self.zoomOnScroll.configure(relief=tk.SUNKEN)
            self.zoom = True

    def onDragActivate(self):
        """Allow drag to move"""
        if self.drag:
            self.dragButton.configure(relief=tk.GROOVE)
            self.drag = False
        else:
            self.dragButton.configure(relief=tk.SUNKEN)
            self.drag = True
            self.filamentPickerAction.configure(relief=tk.GROOVE)
            self.pickerAction.configure(relief=tk.GROOVE)
            self.eraserAction.configure(relief=tk.GROOVE)
            self.eraser = False
            self.picking = False
            self.filament = False

    def zoomerP(self, event):
        """Increase size"""
        direction = 1
        self.onZoomImage(direction, event)

    def zoomerM(self, event):
        """Decrease size"""
        direction = -1
        self.onZoomImage(direction, event)

    def onZoomImage(self, direction, event):
        """Zoom action"""
        if self.zoom:
            self.zoomFactor += direction * 0.1
            self.zoomFactor = max(0.1, min(3.0, self.zoomFactor))
            self.xOffset = 0
            self.yOffset = 0

            if self.zoomFactor >= 1:
                self.updateImage(event)
                if self.particlesWindowVisible:
                    self.extractImages()
            else:
                self.zoomFactor = 1

    def onClickPress(self, event):
        """Actions when the mouse left-click is pressed """
        self.mousePress = True
        self.coordX = self.root.winfo_pointerx() - self.root.winfo_rootx()
        self.coordY = self.root.winfo_pointery() - self.root.winfo_rooty()
        self.dragData['x'] = event.x
        self.dragData['y'] = event.y
        if self.picking and self.canPick(event):
            coordinate_count = int(self.table.item(self.table.selection(), "values")[2])
            shape = self.addCoordinate(event.x * self.scale / self.zoomFactor, event.y * self.scale / self.zoomFactor)
            new_value = coordinate_count + 1
            self.table.set(self.table.selection(), column="Particles", value=new_value)
            self.table.set(self.table.selection(), column="Updated", value='Yes')
            self.coordinatesDict[self.micId][shape] = ((event.x - self.xOffset) * self.scale / self.zoomFactor,
                                                       (event.y - self.yOffset) * self.scale / self.zoomFactor,
                                                       None)
            self.newCoordinates[self.micId][shape] = self.coordinatesDict[self.micId][shape]
            self.totalCoordinates += 1
            self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
            self.hasChanges[self.micId] = True
            if self.particlesWindowVisible:
                self.extractImages()

        findCoord = False
        if not self.eraser:
            for index, coords in self.shapes.items():
                distance = np.sqrt((event.x - self.xOffset - coords[0]) ** 2 + (event.y - self.yOffset - coords[1]) ** 2)
                if distance < self.shapeRadius:
                    findCoord = True
                    indexesToPaint = self.nearCoordinates(index)
                    if self.selectedCoordinate is not None:
                        indexesToRestore = self.nearCoordinates(self.selectedCoordinate)
                        for idx in indexesToRestore:
                            self.imageCanvas.itemconfigure(idx, outline=self.selectedColor)

                    for idx in indexesToPaint:
                        self.imageCanvas.itemconfigure(idx, outline='red')
                        if idx in self.coordinatesDict[self.micId]:
                            self.selectedCoordinate = idx

                    if self.drag:
                        self.root.config(cursor='hand2')

                    if self.particlesWindowVisible:
                        self.updateParticle(index)

            if not findCoord and self.drag:
                if self.selectedCoordinate is not None:
                    indexesToRestore = self.nearCoordinates(self.selectedCoordinate)
                    for idx in indexesToRestore:
                        self.imageCanvas.itemconfigure(idx, outline=self.selectedColor)
                    self.selectedCoordinate = None
                self.root.config(cursor='fleur')
        else:
            self.onMotion(event)

    def onClickRelease(self, event):
        """Actions when the mouse left-click is released """
        self.mousePress = False
        self.root.config(cursor='')
        self.moveShape = False
        if self.particlesWindowVisible and (self.eraser or self.filament):
            self.extractImages()

    def nearCoordinates(self, index):
        """Calculate the near coordinates to a given coordinate"""
        coord = self.imageCanvas.coords(index)
        nearIndex = [index]
        if self.imageCanvas.coords(index):
            coord1 = self.imageCanvas.coords(index - 1)
            if coord == coord1:
                nearIndex.append(index - 1)
        if self.imageCanvas.coords(index+1):
            coord2 = self.imageCanvas.coords(index + 1)
            if coord == coord2:
                nearIndex.append(index + 1)

        return nearIndex

    def onMotion(self, event):
        """Handle the eraser, picking filament and drag actions"""
        x, y = int(event.x) if event.x else None, int(event.y) if event.y else None
        if x is not None and y is not None:
            cursor_target = self.eraser and 'target' or ''
            self.root.config(cursor=cursor_target)

            if 0 < x < self.imageSize[0] and 0 < y < self.imageSize[1]:
                pixel_value = self.getPixelValue(x, y)
                self.infoLabel.config(text=f"x={int(x * self.scale / self.zoomFactor)}, "
                                           f"y={int(y * self.scale / self.zoomFactor)} : {pixel_value}")

            if self.mousePress:
                coordinateCount = int(self.table.item(self.table.selection(), "values")[2])
                if self.eraser:  # Eraser action
                    for index, coords in self.shapes.items():
                        distance = np.sqrt((x - self.xOffset - coords[0]) ** 2 + (y - self.yOffset - coords[1]) ** 2)
                        if distance < self.shapeRadius:
                            new_value = coordinateCount - 1
                            self.table.set(self.table.selection(), column='Particles', value=new_value)
                            self.table.set(self.table.selection(), column='Updated', value='Yes')
                            indexesToDelete = self.nearCoordinates(index)
                            for idx in indexesToDelete:
                                self.imageCanvas.delete(idx)
                                self.shapes.pop(idx)

                            if index in self.coordinatesDict[self.micId]:
                                coordObjId = self.coordinatesDict[self.micId][index][2]
                                self.deletedCoordinates[self.micId][coordObjId] = True
                                self.coordinatesDict[self.micId].pop(index)

                            self.totalCoordinates -= 1
                            self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
                            self.hasChanges[self.micId] = True
                            break

                elif self.filament and self.canPick(event):  # Filament picking action
                    shape = self.addCoordinate(x * self.scale / self.zoomFactor, y * self.scale / self.zoomFactor)
                    new_value = coordinateCount + 1
                    self.table.set(self.table.selection(), column="Particles", value=new_value)
                    self.table.set(self.table.selection(), column="Updated", value='Yes')
                    self.coordinatesDict[self.micId][shape] = ((x - self.xOffset) * self.scale / self.zoomFactor,
                                                               (y - self.yOffset) * self.scale / self.zoomFactor,
                                                               None)
                    self.newCoordinates[self.micId][shape] = self.coordinatesDict[self.micId][shape]
                    self.totalCoordinates += 1
                    self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
                    self.hasChanges[self.micId] = True

                elif self.drag or (self.picking and not self.canPick(event)):  # Move coordinate or drag all image
                    if self.selectedCoordinate:
                        index, coords = self.selectedCoordinate, self.shapes[self.selectedCoordinate]
                        self.root.config(cursor='hand2')
                        indexesToMove = self.nearCoordinates(index)
                        newX = self.root.winfo_pointerx() - self.root.winfo_rootx() - self.coordX
                        newY = self.root.winfo_pointery() - self.root.winfo_rooty() - self.coordY

                        if self.isMoveIn(coords[0] + newX, coords[1] + newY):
                            self.moveShape = True
                            for idx in indexesToMove:
                                self.shapes[idx] = (self.shapes[idx][0] + newX, self.shapes[idx][1] + newY)
                                # Move the shape to a new position
                                self.imageCanvas.move(idx, newX, newY)

                            # Update de coordinates
                            self.coordX += newX
                            self.coordY += newY
                            coordXY = self.coordinatesDict[self.micId][self.selectedCoordinate]
                            self.coordinatesDict[self.micId][self.selectedCoordinate] = (coordXY[0] + newX * self.scale / self.zoomFactor,
                                                                                         coordXY[1] + newY * self.scale / self.zoomFactor,
                                                                                         coordXY[2])
                            self.movedCoordinates[self.micId][coordXY[2]] = self.coordinatesDict[self.micId][self.selectedCoordinate]

                            if self.selectedCoordinate is not None:
                                indexesToPaint = self.nearCoordinates(self.selectedCoordinate)
                                for idx in indexesToPaint:
                                    self.imageCanvas.itemconfigure(idx, outline=self.selectedColor)

                            if self.particlesWindowVisible:
                                self.moveParticle(index)

                            for idx in indexesToMove:
                                self.imageCanvas.itemconfigure(idx, outline='red')

                            if newX != 0 or newY != 0:
                                self.table.set(self.table.selection(), column="Updated", value='Yes')
                                self.hasChanges[self.micId] = True

                    if not self.moveShape:  # Drag the image and the shapes
                        self.root.config(cursor='fleur')
                        self.moveShape = False
                        self.onDrag(event)
        else:
            self.root.config(cursor='')

    def canPick(self, event):
        """Returns true if picking outside of any coordinate"""
        canPick = True
        for index, coords in self.shapes.items():
            distance = np.sqrt((event.x - self.xOffset - coords[0]) ** 2 + (event.y - self.yOffset - coords[1]) ** 2)
            if distance < self.shapeRadius:
                canPick = False
                break
        return canPick

    def isMoveIn(self, x, y):
        if 0 < x < self.imageSize[0] / self.scale * self.zoomFactor and 0 < y < self.imageSize[1] / self.scale * self.zoomFactor:
            return True
        return False

    def onDrag(self, event):
        """Move the image and shapes on drag"""
        x = (event.x - self.dragData['x'])
        y = (event.y - self.dragData['y'])
        self.xOffset += x
        self.yOffset += y
        self.imageCanvas.move('all', x, y)
        self.dragData['x'] = event.x
        self.dragData['y'] = event.y

    def getPixelValue(self, x, y):
        """Return the pixel value in the given x and y coordinates"""
        value = self.imagePIL.getpixel((x, y))
        return str(value)

    def updateImage(self, event):
        """Update the image with new size"""
        self.imageCanvas.delete("all")
        width, height = self.imagePIL.size
        new_width = int(width * self.zoomFactor)
        new_height = int(height * self.zoomFactor)
        self.scaledImage = self.imagePIL.resize((int(new_width / self.scale), int(new_height / self.scale)))
        self.imageTk = ImageTk.PhotoImage(self.scaledImage)
        self.imageCanvas.delete("all")
        self.imageCanvas.config(scrollregion=self.imageCanvas.bbox("all"))
        self.image = self.imageCanvas.create_image(0, 0, anchor=tk.NW, image=self.imageTk, tags='image')
        self.drawCoordinates(self.micId)
        if self.zoomFactor >= 1:
            x, y = ((event.x - self.dragData['x']) / self.scale * self.zoomFactor,
                    (event.y - self.dragData['y']) / self.scale * self.zoomFactor)
            self.xOffset -= x
            self.yOffset -= y
            self.imageCanvas.move('all', -x, -y)
            self.dragData['x'] = event.x
            self.dragData['y'] = event.y

    def showCoordinates(self):
        """Load and display the selected micrograph and coordinates"""
        imagePath = os.path.abspath(self.micrographPathDict[self.micId][0])
        if imagePath:
            try:
                imageReader = ImageReadersRegistry.getReader(imagePath)
                self.imagePIL = imageReader.open(imagePath)
                self.imageSize = self.imagePIL.size

                self.scale = max(self.imageSize[0]/self.imageCanvasSize, self.imageSize[1]/self.imageCanvasSize)
                dpiWidth = self.imageSize[0] / self.scale
                dpiHeight = self.imageSize[1] / self.scale
                self.imageCanvas.configure(width=dpiWidth, height=dpiHeight)

                if self.enhanceContrastVar.get():
                    contrast = ImageEnhance.Contrast(self.imagePIL)
                    self.imagePIL = contrast.enhance(8)

                if self.gaussianBlurVar.get():
                    self.imagePIL = self.imagePIL.filter(ImageFilter.GaussianBlur(radius=0.5))

                if self.invertContrastVar.get():
                    self.imagePIL = ImageOps.invert(self.imagePIL)

                self.scaledImage = self.imagePIL.resize((int(dpiWidth), int(dpiHeight)))
                self.imageTk = ImageTk.PhotoImage(self.scaledImage)

                # self.quadtree = Index(bbox=[0, 0, self.imagePIL.size[0], self.imagePIL.size[1]])
                self.imageCanvas.delete("all")
                self.image = self.imageCanvas.create_image(0, 0, anchor=tk.NW, image=self.imageTk, tags='image')
                self.zoomFactor = 1
                self.drawCoordinates(self.micId)

            except Exception as e:
                logger.error(f"Error loading image '{self.micId}': {e}")
        else:
            logger.error("Unable to upload the file %s. Make sure the path is correct." % imagePath)

    def applyPowerHistogram(self):
        """ Create the histogram and the plot"""
        self.histWindow = tk.Toplevel(self.root)
        self.histWindow.title('Power histogram')
        self.histWindow.resizable(False, False)
        self.figure, self.axes = plt.subplots(figsize=(6, 4))
        self.plotHistogram(self.axes)
        histCanvas = FigureCanvasTkAgg(self.figure, master=self.histWindow)
        histCanvasWidget = histCanvas.get_tk_widget()
        histCanvasWidget.grid(row=0, column=0, sticky="news", padx=10, pady=10)

        def update(val):
            minValue, maxValue = self.rangeSlider.val
            self.drawRangeLines(self.axes, minValue, maxValue)
            plt.draw()
            self.removeCoordinates(minValue, maxValue)

        self.rangeLines = []
        ax_slider = plt.axes([0.2, 0.05, 0.65, 0.03], facecolor='lightgoldenrodyellow')
        self.rangeSlider = RangeSlider(ax_slider, '', 0, 256, valinit=(0, 256), valstep=1, valfmt='%d')

        self.rangeSlider.on_changed(update)

        buttonsFrame = tk.Frame(self.histWindow)
        buttonsFrame.grid(row=3, column=0, sticky="e", pady=10)
        # Close button
        closeButton = tk.Button(buttonsFrame, text='Close', command=self.histWindowClose, font=('Helvetica', 10, 'bold'),
                                image=getImage(Icon.BUTTON_CANCEL), compound=tk.LEFT)
        closeButton.grid(row=3, column=0, sticky="e", padx=5)

        # Save button
        saveButton = tk.Button(buttonsFrame, text='Save & Close', command=self.histWindowSaveClose, fg='white',
                               font=('Helvetica', 10, 'bold'), image=getImage(Icon.ACTION_SAVE), compound=tk.LEFT,
                               bg='#B22A2A', activebackground='#B22A2A')
        saveButton.grid(row=3, column=1, sticky="w", padx=5)

    def histWindowClose(self):
        self.histWindow.destroy()
        self.showCoordinates()
        self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
        self.table.set(self.table.selection(), column="Particles", value=len(self.coordinatesDict[self.micId]))

    def histWindowSaveClose(self):
        """Save the new coordinates taking into account the sliders selected pixel range """
        # Updating current micrograph
        shapeToDelete = []
        for index, coord in self.shapes.items():
            currentState = self.imageCanvas.itemcget(index, "state")
            if currentState == 'hidden':
                shapeToDelete.append(index)
                self.imageCanvas.delete(index)
                coordObjId = self.coordinatesDict[self.micId][index][2]
                self.deletedCoordinates[self.micId][coordObjId] = True
                self.coordinatesDict[self.micId].pop(index)
                self.totalCoordinates -= 1

        for i in range(len(shapeToDelete)):
            self.shapes.pop(shapeToDelete[i])

        self.histWindowClose()

    def removeCoordinates(self, value1, value2):
        """Remove the coordinate taking into account the pixel values"""
        pixelRange = (round(float(value1)), round(float(value2)))
        for index, coords in self.shapes.items():
            pixelValue = self.calculateAveragePixel(coords[0], coords[1])
            new_state = "normal"
            if int(pixelValue) < pixelRange[0] or int(pixelValue) > pixelRange[1]:
                new_state = "hidden"
            self.imageCanvas.itemconfigure(index, state=new_state)

    def calculateAveragePixel(self, shapeX, shapeY):
        """Calculate the shape pixel average"""
        # Take a region delimited by the shape(circle or square)
        region = self.extractRegion(shapeX, shapeY)

        pixelArray = np.array(region)
        # Calculate the mean of RGB canals
        averagePixel = np.mean(pixelArray, axis=(0, 1))

        return averagePixel

    def extractRegion(self, x, y):
        """Take the particle region"""
        return self.scaledImage.crop((x - self.shapeRadius,  y - self.shapeRadius,
                                     x + self.shapeRadius, y + self.shapeRadius))

    def plotHistogram(self, axes):
        """Plot the micrograph histogram"""
        img_array = np.array(self.imagePIL)
        # Apply the power transformation
        img_power = np.power(img_array / 255.0, 2.0) * 255.0
        hist_power, bins_power = np.histogram(img_power.flatten(), bins=256, range=[0, 256])

        # Plot histogram
        axes.clear()
        axes.plot(hist_power)
        axes.set_xlabel("Pixel Value")
        axes.set_ylabel("Frequency")
        plt.subplots_adjust(bottom=0.2)

    def drawRangeLines(self, ax, minVal, maxVal):
        """Draw vertical lines to represent the selected range"""
        # Remove old lines
        for line in self.rangeLines:
            line.remove()
        # Draw the new vertical lines
        lineMin = ax.axvline(minVal, color='red', linestyle='--', linewidth=1)
        lineMax = ax.axvline(maxVal, color='red', linestyle='--', linewidth=1)
        # Store the new vertical lines
        self.rangeLines = [lineMin, lineMax]

    def resetMicrograph(self):
        """Remove all coordinates from current micrograph"""
        result = messagebox.askquestion("Confirmation", "Are you sure you want to reset the micrograph?",
                                        icon='warning', **{'parent': self.root})
        if result == messagebox.YES:
            self.totalCoordinates -= len(self.coordinatesDict[self.micId])
            self.coordinatesDict[self.micId] = {}
            self.shapes.clear()
            self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
            self.table.set(self.table.selection(), column="Particles", value=0)
            self.imageCanvas.delete("shape")
            self.table.set(self.table.selection(), column='Updated', value='Yes')
            self.hasChanges[self.micId] = True
            if self.particlesWindowVisible:
                self.onParticlesWindowClosing()

    def restoreMicrograph(self):
        """Restore all removed coordinates from current micrograph"""
        result = messagebox.askquestion("Confirmation", "Are you sure you want to restore the micrograph?",
                                        icon='warning', **{'parent': self.root})
        if result == messagebox.YES:
            self.totalCoordinates -= len(self.coordinatesDict[self.micId])
            self.coordinatesDict[self.micId] = self.oldCoordinatesDict[self.micId]
            self.deletedCoordinates[self.micId] = {}
            self.newCoordinates[self.micId] = {}
            self.movedCoordinates[self.micId] = {}
            self.totalCoordinates += len(self.coordinatesDict[self.micId])
            self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
            self.imageCanvas.delete("shape")
            self.table.set(self.table.selection(), column="Particles", value=len(self.coordinatesDict[self.micId]))
            self.table.set(self.table.selection(), column='Updated', value='No')
            self.drawCoordinates(self.micId)
            if self.micId in self.hasChanges:
                self.hasChanges.pop(self.micId)
            if self.particlesWindowVisible:
                self.extractImages()

    def showParticles(self, geometry=None):
        """Create the particles window"""
        if not self.particlesWindowVisible:
            self.particlesWindow = tk.Toplevel(self.root)
            self.particlesWindow.title('Particles')
            self.particlesWindow.attributes('-topmost', True)
            self.particlesWindowVisible = True
            self.particlesWindow.protocol("WM_DELETE_WINDOW", self.onParticlesWindowClosing)
            self.particlesCanvas = tk.Canvas(self.particlesWindow)
            self.particlesCanvas.pack(side="left", fill="both", expand=True)

            self.scrollbar = tk.Scrollbar(self.particlesWindow, orient="vertical", command=self.particlesCanvas.yview)
            self.scrollbar.pack(side="right", fill="y")

            self.imageGrid = tk.Frame(self.particlesCanvas, bg='white')
            self.particlesCanvas.create_window((0, 0), window=self.imageGrid, anchor="nw")
            self.particlesCanvas.configure(yscrollcommand=self.scrollbar.set)
            self.particlesCanvas.bind("<Configure>", lambda e: self.particlesCanvas.configure(scrollregion=self.particlesCanvas.bbox("all")))

            self.selectedParticle = None
            self.extractImages()

    def extractImages(self):
        """Show the particles"""
        self.clearImageGrid()
        self.particlesIndex = dict()
        self.particlesWidget = dict()
        self.image_references = []
        self.selectedParticle = None
        self.particlesWindow.geometry("%sx800" % str(3*self.boxSize + 60))
        row = column = count = 0
        for index, coords in self.shapes.items():
            if index in self.coordinatesDict[self.micId]:
                region = self.extractRegion(coords[0], coords[1])
                scaledImage = region.resize((self.boxSize, self.boxSize))
                contrast = ImageEnhance.Contrast(scaledImage)
                scaledImage = contrast.enhance(1)
                imageTk = ImageTk.PhotoImage(scaledImage)
                self.image_references.append(imageTk)

                label = tk.Label(self.imageGrid, image=self.image_references[count], width=self.boxSize + 5,
                                 height=self.boxSize + 5, borderwidth=1, highlightthickness=1, highlightbackground='blue')
                label.image = self.image_references[count]
                label.bind("<Button-1>", self.selectCoordinate)
                self.particlesIndex[label] = index
                self.particlesWidget[index] = label

                label.image = self.image_references[count]
                label.grid(row=row, column=column, padx=3)
                column += 1
                if column % 3 == 0:  # always display 3 columns
                    row += 1
                    column = 0
                count += 1

        self.particlesCanvas.update_idletasks()
        scroll_height = self.particlesCanvas.bbox("all")[3]
        self.scrollbar.set(0, 1)
        self.particlesCanvas.config(scrollregion=(0, 0, 1, scroll_height))

    def clearImageGrid(self):
        for widget in self.imageGrid.winfo_children():
            widget.destroy()
        self.selectedParticle = None

    def selectCoordinate(self, event):
        index = self.particlesIndex[event.widget]

        if self.selectedParticle is not None:
            self.selectedParticle.config(highlightbackground='blue', bg=self.root.cget("bg"))

        self.selectedParticle = event.widget

        if self.selectedCoordinate is not None:
            indexesToRestore = self.nearCoordinates(self.selectedCoordinate)
            for idx in indexesToRestore:
                self.imageCanvas.itemconfigure(idx, outline=self.selectedColor)

        self.selectedCoordinate = index
        indexesToPaint = self.nearCoordinates(self.selectedCoordinate)
        for idx in indexesToPaint:
            self.imageCanvas.itemconfigure(idx, outline='red')

        self.selectedParticle.config(highlightbackground='red', bg='red')

    def onParticlesWindowClosing(self):
        """Close the particles window"""
        self.particlesWindowVisible = False
        self.selectedParticle = None
        self.particlesWindow.destroy()

    def updateParticle(self, index):
        if index not in self.particlesWidget:
            index = index - 1
        labelWidget = self.particlesWidget[index]
        if self.selectedParticle is not None:
            self.selectedParticle.config(highlightbackground='blue', bg=self.root.cget("bg"))
        self.selectedParticle = labelWidget

        yPosition = self.selectedParticle.winfo_y()
        yview_position = yPosition / self.imageGrid.winfo_height()
        self.particlesCanvas.yview_moveto(yview_position)
        self.selectedParticle.config(highlightbackground='red', bg='red')

    def moveParticle(self, index):
        newCoords = self.shapes[index]
        region = self.extractRegion(newCoords[0], newCoords[1])
        scaledImage = region.resize((self.boxSize, self.boxSize))
        contrast = ImageEnhance.Contrast(scaledImage)
        scaledImage = contrast.enhance(1)
        imageTk = ImageTk.PhotoImage(scaledImage)
        self.image_references.append(imageTk)
        self.selectedParticle.config(image=self.image_references[-1])

    def onSliderMapMove(self, value):
        self.infoLabel.config(text=f"{float(value):.0f}")

    def drawCoordinates(self, micId):
        """Draw the coordinates over the micrograph"""
        self.shapes = {}
        coordinates = self.coordinatesDict[micId]
        self.shapeRadius = self.boxSize / self.scale / 2 * self.zoomFactor
        self.auxCoordinatesDict = dict()
        for index, coord in coordinates.items():
            self.addCoordinate(coord[0], coord[1], coord[2])
        if self.auxCoordinatesDict:
            self.coordinatesDict[micId] = self.auxCoordinatesDict

    def addCoordinate(self, x, y, coordId=None):
        """Create a coordinate"""
        xTrans, yTrans = x / self.scale * self.zoomFactor, y / self.scale * self.zoomFactor
        shape = None
        if self.drawCircles:
            circle = self.imageCanvas.create_oval(xTrans - self.shapeRadius,
                                                  yTrans - self.shapeRadius,
                                                  xTrans + self.shapeRadius,
                                                  yTrans + self.shapeRadius,
                                                  outline=self.selectedColor, width=1, fill="", tags='shape')
            self.shapes[circle] = (xTrans, yTrans)
            self.auxCoordinatesDict[circle] = (x, y, coordId)
            shape = circle
            # self.quadtree.insert(circle, (x, y))

        if self.drawSquares:
            square = self.imageCanvas.create_rectangle(xTrans - self.shapeRadius,
                                                       yTrans - self.shapeRadius,
                                                       xTrans + self.shapeRadius,
                                                       yTrans + self.shapeRadius,
                                                       outline=self.selectedColor, width=1, fill="", tags='shape')
            self.shapes[square] = (xTrans, yTrans)
            if not self.drawCircles:
                self.auxCoordinatesDict[square] = (x, y, coordId)
                shape = square
        return shape

    def openFile(self):
        """Open a file dialog to load an image"""
        file_path = filedialog.askopenfilename(filetypes=[("Image files", "*.png;*.jpg;*.jpeg;*.gif")])
        if file_path:
            pass

    def updateSize(self, event=None):
        """Update the size value label when the slider is moved"""
        self.sizeValueLabel.config(text=f"({int(self.sizeVar.get())})")

    def applyContrast(self):
        self.enhanceContrastVar.set(not self.enhanceContrastVar.get())
        if self.enhanceContrastVar.get():
            self.filterMenu.entryconfigure(0, image=getImage(Icon.CHECKED))
        else:
            self.filterMenu.entryconfigure(0, image=getImage(Icon.UNCHECKED))
        self.showCoordinates()
        if self.particlesWindowVisible:
            self.extractImages()

    def applyGaussianBlur(self):
        """Apply the gaussian blur image filter"""
        self.gaussianBlurVar.set(not self.gaussianBlurVar.get())
        if self.gaussianBlurVar.get():
            self.filterMenu.entryconfigure(1, image=getImage(Icon.CHECKED))
        else:
            self.filterMenu.entryconfigure(1, image=getImage(Icon.UNCHECKED))
        self.showCoordinates()
        if self.particlesWindowVisible:
            self.extractImages()

    def applyInvertContrast(self):
        """Invert the image contrast"""
        self.invertContrastVar.set(not self.invertContrastVar.get())
        if self.invertContrastVar.get():
            self.filterMenu.entryconfigure(2, image=getImage(Icon.CHECKED))
        else:
            self.filterMenu.entryconfigure(2, image=getImage(Icon.UNCHECKED))
        self.showCoordinates()
        if self.particlesWindowVisible:
            self.extractImages()

    def applyBandpassFilter(self):
        """Apply the bandpass image filter"""
        self.bandpassFilterVar = not self.bandpassFilterVar
        self.showCoordinates()
        if self.particlesWindowVisible:
            self.extractImages()


class CoordinateView(View):

    def __init__(self, root, emSet: EMSet, protocol):
        self.root = root
        self._emSet = emSet
        self.protocol = protocol

    def show(self):
        root = tk.Toplevel(self.root)
        app = MainWindow(root, self._emSet, self.protocol)
        width = 570 + app.scaledImage.size[0]
        height = 140 + app.scaledImage.size[1]
        root.geometry(f"{width}x{height}")
        root.mainloop()


class CoordinateViewer(Viewer):
    _name = 'Scipion coordinates viewer'
    _targets = [SetOfCoordinates]

    def _visualize(self, obj, **kwargs):
        return [CoordinateView(self.getTkRoot(), obj, self.protocol)]

