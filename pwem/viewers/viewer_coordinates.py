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
from pwem.emlib.image.image_readers import ImageReadersRegistry
from pwem.objects import EMSet, SetOfCoordinates
from pyworkflow.gui import getImage, ToolTip
from pyworkflow.utils import Icon
from pyworkflow.viewer import View, Viewer

logger = logging.getLogger()


class MainWindow:
    """Coordinate viewer"""
    def __init__(self, root, setOfCoordinate):
        # Initialize the main window
        self.root = root
        self.setOfCoordinate = setOfCoordinate
        self.root.title("Scipion coordinates viewer")
        self.root.resizable(False, False)
        self.initGui()

    def initGui(self):
        self.micName = None
        self.boxSize = self.setOfCoordinate.getBoxSize()
        self.selectedColor = '#00FF00'
        self.circleButtonRelieve = tk.SUNKEN
        self.squareButtonRelieve = tk.GROOVE
        self.zoomButtonRelieve = tk.SUNKEN
        self.dragButtonReliev = tk.SUNKEN
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
        self.scale = 5
        self.drawSquares = False
        self.drawCircles = True
        self.zoomFactor = 1
        self.moveShape = False
        self.hasChanges = {}
        self.drag_data = {'x': 0, 'y': 0, 'item': None}
        self.xOffset = 0
        self.yOffset = 0

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

        # Main frame setup
        self.mainFrame = ttk.Frame(self.root)
        self.mainFrame.grid(row=0, column=0, sticky="news")
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

        # Toolbar setup
        self.createToolbar(self.mainFrame)

        # Frame divided in two
        self.contentFrame = ttk.Frame(self.mainFrame)
        self.contentFrame.grid(row=2, column=0, sticky="news")

        # Table on the left
        self.createTable()

        # Micrograph on the right
        self.createImageFrame()
        self.loadMicrograph()

    def createToolbar(self, parent):
        """Create the toolbar with a slider for adjusting the box size"""
        toolbar = ttk.Frame(parent)
        toolbar.bind("<Enter>", self.recoveryPointer)

        sizeLabel = ttk.Label(toolbar, text="Box size(px):", image=getImage(Icon.FILE_VOL), compound=tk.LEFT)
        sizeLabel.grid(row=0, column=0, padx=10, pady=5, sticky="e")

        self.sizeVar = tk.IntVar(parent, self.boxSize)
        self.sizeSlider = ttk.Scale(toolbar, from_=0, to=self.boxSize*2, orient=tk.HORIZONTAL, length=self.boxSize*2,
                                    variable=self.sizeVar,  command=self.updateSize)
        self.sizeSlider.grid(row=0, column=1, padx=5, pady=5)
        self.sizeSlider.bind("<B1-Motion>", self.onBoxSizeSlider)  # Bind event for box size slider motion
        # self.sizeSlider.bind("<ButtonRelease-1>", self.onSliderRelease)

        self.sizeValueLabel = ttk.Label(toolbar, text=f"({int(self.sizeVar.get())})")
        self.sizeValueLabel.grid(row=0, column=2, padx=5, pady=5, sticky="w")

        self.totalMicrographButton = tk.Button(toolbar, bg='#CFE8CF', relief=tk.SUNKEN,
                                               activebackground='#CFE8CF', compound=tk.LEFT)
        self.totalMicrographButton.grid(row=0, column=3, padx=20, pady=5, sticky="e")

        self.totalCoordinates = self.setOfCoordinate.getSize()
        self.totalPickButton = tk.Button(toolbar, text=f"Total picks: {self.totalCoordinates}",
                                         bg='#C3D7DF', relief=tk.SUNKEN,
                                         activebackground='#C3D7DF', compound=tk.LEFT)

        self.totalPickButton.grid(row=0, column=4, padx=0, pady=5, sticky="e")

        self.toolbar2 = ttk.Frame(parent)
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
        self.pickerLabel = ttk.Label(self.toolbar2, text="Picker tools:")
        self.pickerLabel.grid(row=0, column=8, padx=5, pady=5, sticky="e")
        self.pickerLabel.pack_forget()

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
        imageToolsPanel = ttk.Frame(self.toolbar2)
        imageToolsPanel.grid(row=0, column=13, padx=450, pady=0, sticky="e")

        self.fitDisplay = tk.Button(imageToolsPanel, width=25, height=25, command=self.onFitActivate,
                                    image=getImage(Icon.ACTION_EXPAND),
                                    relief=self.squareButtonRelieve)
        self.fitDisplay.grid(row=0, column=12, padx=0, pady=0, sticky="ns")
        tooltip = "Fit to display area"
        ToolTip(self.fitDisplay, tooltip, delay=150)

        self.zoomOnScroll = tk.Button(imageToolsPanel, width=25, height=25, command=self.onZoomActivate,
                                    image=getImage(Icon.ACTION_ZOOM),
                                    relief=self.zoomButtonRelieve)
        self.zoomOnScroll.grid(row=0, column=14, padx=0, pady=0, sticky="ns")
        tooltip = "Zoom on scroll"
        ToolTip(self.zoomOnScroll, tooltip, delay=150)

        self.dragButton = tk.Button(imageToolsPanel, width=25, height=25, command=self.onDragActivate,
                                      image=getImage(Icon.ACTION_HAND),
                                      relief=self.dragButtonReliev)
        self.dragButton.grid(row=0, column=13, padx=0, pady=0, sticky="ns")
        tooltip = "Click and drag to move"
        ToolTip(self.dragButton, tooltip, delay=150)

        # Label to display coordinates and pixel value
        self.infoLabel = ttk.Label(imageToolsPanel, text="")
        self.infoLabel.grid(row=0, column=17, sticky="w", padx=120)

        toolbar.grid(row=0, column=0, sticky="ew")
        self.toolbar2.grid(row=1, column=0, sticky="ew")

    def recoveryPointer(self, event):
        self.root.config(cursor='')
        self.infoLabel.config(text='')

    def pickingActivate(self):
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
            self.drag = False

    def eraserActivate(self):
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
            self.drag = False

    def filamentActivate(self):
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
        self.shapeRadius = self.boxSize / self.scale / 2 * self.zoomFactor
        self.auxCoordinatesDict = dict()
        for index, coord in self.shapes.items():
            self.imageCanvas.coords(index, [coord[0] - self.shapeRadius, coord[1] - self.shapeRadius,
                                            coord[0] + self.shapeRadius, coord[1] + self.shapeRadius])

    def onSliderRelease(self, event):
        """Handle box size slider release"""
        self.loadMicrograph()

    def showCircles(self):
        """Draw a circle on the particle"""
        if self.circleButtonRelieve == tk.GROOVE:
            self.circleButtonRelieve = tk.SUNKEN
            self.drawCircles = True
        else:
            self.circleButtonRelieve = tk.GROOVE
            self.drawCircles = False

        self.circleShape.configure(relief=self.circleButtonRelieve)
        self.loadMicrograph()

    def showSquare(self):
        """Draw a square on the particle"""
        if self.squareButtonRelieve == tk.GROOVE:
            self.squareButtonRelieve = tk.SUNKEN
            self.drawSquares = True
        else:
            self.squareButtonRelieve = tk.GROOVE
            self.drawSquares = False

        self.squareShape.configure(relief=self.squareButtonRelieve)
        self.loadMicrograph()

    def showPoint(self):
        """Draw a point on the particle"""
        if self.pointButtonRelieve == tk.GROOVE:
            self.pointButtonRelieve = tk.SUNKEN
        else:
            self.pointButtonRelieve = tk.GROOVE

        self.pointShape.configure(relief=self.pointButtonRelieve)
        self.loadMicrograph()

    def showPalette(self):
        color = tkcolorpicker.ColorPicker(self.root)
        self.root.wait_window(color)
        if color.get_color():
            self.selectedColor = color.get_color()[2]
            self.paletteButton.configure(bg=self.selectedColor, activebackground=self.selectedColor)
            self.loadMicrograph()

    def createTable(self):
        """Create a table on the left side to display information"""
        columns = ('Index', 'File', 'Particles', 'Updated')
        data = []

        setOfMicrographs = self.setOfCoordinate.getMicrographs()
        for micrograph in setOfMicrographs.iterItems():
            micBaseName = os.path.basename(micrograph.getFileName())
            self.coordinatesDict[micBaseName] = {}
            self.oldCoordinatesDict[micBaseName] = {}
            self.micrographPathDict[micBaseName] = micrograph.getFileName()

        for index, coordinate in enumerate(self.setOfCoordinate.iterItems()):
            micName = coordinate.getMicName()
            if len(micName.split('.')) == 1:  # ensuring that it contains an extension
                micName += '.mrc'
            self.coordinatesDict[micName][index+2] = (coordinate.getX(), coordinate.getY())
            self.oldCoordinatesDict[micName][index+2] = (coordinate.getX(), coordinate.getY())

        for mic, micrograph in enumerate(self.coordinatesDict):
            data.append((mic + 1, micrograph, len(self.coordinatesDict[micrograph]), 'No'))

        self.table = ttk.Treeview(self.contentFrame, columns=columns, show="headings")
        self.table.bind("<Enter>", self.recoveryPointer)
        self.table.bind("<ButtonRelease-1>", self.onTableClick)  # Bind event for table click
        for col, width in zip(columns, (50, 300, 100, 100)):
            self.table.heading(col, text=col, anchor="center")
            self.table.column(col, anchor="center", width=width)

        for row in data:
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
            values = self.table.item(selected_item, "values")
            if len(values) > 1 and values[1] != self.micName:
                self.micName = values[1]
                self.zoomFactor = 1
                self.selectedCoordinate = None
                self.loadMicrograph()

    def createImageFrame(self):
        """Create the image frame on the right side to display micrographs and associated coordinates"""
        self.imageFrame = ttk.Frame(self.contentFrame)
        self.imageFrame.grid(row=0, column=2, sticky="nw")
        self.contentFrame.grid_rowconfigure(0, weight=1)
        self.contentFrame.grid_columnconfigure(1, weight=1)

        self.imageCanvas = tk.Canvas(self.imageFrame, width=780, height=780)  #, borderwidth=5, highlightthickness=1, highlightbackground='red')
        self.imageCanvas.grid(row=0, column=0, sticky="nw", pady=10)
        self.imageCanvas.bind("<Motion>", self.onPickerEraserAction)
        self.imageCanvas.bind("<B1-Motion>", self.onPickerEraserAction)
        self.imageCanvas.bind("<Button-1>", self.onClickPress)
        self.imageCanvas.bind('<ButtonRelease-1>', self.onClickRelease)
        self.imageCanvas.bind("<Button-4>", self.zoomerP)
        self.imageCanvas.bind("<Button-5>", self.zoomerM)

        buttonsFrame = tk.Frame(self.imageFrame)
        buttonsFrame.grid(row=1, column=0, sticky="ns", padx=510)

        # Button to close de application
        closeButton = tk.Button(buttonsFrame, text='Close', command=self.root.destroy, font=('Helvetica', 10, 'bold'),
                                      image=getImage(Icon.BUTTON_CANCEL), compound=tk.LEFT)
        closeButton.grid(row=0, column=0, sticky="w", padx=5)

        # Button to create a new set of Coordinates
        coordinateButton = tk.Button(buttonsFrame, text='Coordinate', fg='white', font=('Helvetica', 10, 'bold'),
                                     image=getImage(Icon.ACTION_NEW), compound=tk.LEFT, bg='#B22A2A',
                                     activebackground='#B22A2A')
        coordinateButton.grid(row=0, column=1, sticky="w", padx=5)

        self.imageFrame.grid_rowconfigure(0, weight=1)
        self.imageFrame.grid_columnconfigure(0, weight=1)

    def onFitActivate(self):
        self.zoomFactor = 1
        self.loadMicrograph()
        self.drag_data = {'x': 0, 'y': 0, 'item': None}
        self.xOffset = 0
        self.yOffset = 0

    def onZoomActivate(self):
        if self.zoom:
            self.zoomOnScroll.configure(relief=tk.GROOVE)
            self.zoom = False
        else:
            self.zoomOnScroll.configure(relief=tk.SUNKEN)
            self.zoom = True

    def onDragActivate(self):
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

            if self.zoomFactor < 1:
                self.zoomFactor = 1

            self.updateImage(event)

    def onClickPress(self, event):
        self.mousePress = True
        self.coordX = self.root.winfo_pointerx() - self.root.winfo_rootx()
        self.coordY = self.root.winfo_pointery() - self.root.winfo_rooty()
        self.drag_data['x'] = event.x
        self.drag_data['y'] = event.y
        if self.picking:
            pick = True
            for index, coords in self.shapes.items():
                distance = np.sqrt((event.x - self.xOffset - coords[0]) ** 2 + (event.y - self.yOffset - coords[1]) ** 2)
                if distance < self.shapeRadius:
                    pick = False
                    break

            if pick:
                coordinate_count = int(self.table.item(self.table.selection(), "values")[2])
                shape = self.addCoordinate(event.x * self.scale / self.zoomFactor, event.y * self.scale / self.zoomFactor)
                new_value = coordinate_count + 1
                self.table.set(self.table.selection(), column="Particles", value=new_value)
                self.table.set(self.table.selection(), column="Updated", value='Yes')
                self.coordinatesDict[self.micName][shape] = ((event.x - self.xOffset) * self.scale / self.zoomFactor,
                                                              (event.y - self.yOffset) * self.scale / self.zoomFactor)
                self.totalCoordinates += 1
                self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
                self.hasChanges[self.micName] = True

        findCoord = False
        for index, coords in self.shapes.items():
            distance = np.sqrt((event.x - self.xOffset - coords[0]) ** 2 + (event.y - self.yOffset - coords[1]) ** 2)
            if distance < self.shapeRadius:
                findCoord = True
                indexesToPaint = self.nearCoordinates(index)
                if self.selectedCoordinate is not None:
                    indexesToRestore = self.nearCoordinates(self.selectedCoordinate)
                    for idx in indexesToRestore:
                        self.imageCanvas.itemconfigure(idx, outline=self.selectedColor)
                self.selectedCoordinate = index
                for idx in indexesToPaint:
                    self.imageCanvas.itemconfigure(idx, outline='red')

                if self.drag:
                    self.root.config(cursor='hand2')
        if not findCoord and self.drag:
            self.root.config(cursor='fleur')

    def onClickRelease(self, event):
        self.mousePress = False
        self.root.config(cursor='')
        self.moveShape = False

    def nearCoordinates(self, index):
        coord = self.imageCanvas.coords(index)
        coord1 = self.imageCanvas.coords(index - 1)
        coord2 = self.imageCanvas.coords(index + 1)
        nearIndex = [index]
        if coord == coord1:
            nearIndex.append(index - 1)
        if coord == coord2:
            nearIndex.append(index + 1)

        return nearIndex

    def onPickerEraserAction(self, event):
        """Handle the eraser and picking action"""
        x, y = int(event.x) if event.x else None, int(event.y) if event.y else None
        if x is not None and y is not None:
            cursor_target = self.eraser and 'target' or ''
            self.root.config(cursor=cursor_target)

            if 0 < x < self.imageSize[0] and 0 < y < self.imageSize[1]:
                pixel_value = self.getPixelValue(x, y)
                self.infoLabel.config(text=f"x={int(x * self.scale / self.zoomFactor)}, "
                                           f"y={int(y * self.scale / self.zoomFactor)} : {pixel_value}")

            if self.mousePress:
                coordinate_count = int(self.table.item(self.table.selection(), "values")[2])
                if self.eraser:  # Eraser action
                    for index, coords in self.shapes.items():
                        distance = np.sqrt((x - self.xOffset - coords[0]) ** 2 + (y - self.yOffset - coords[1]) ** 2)
                        if distance < self.shapeRadius:
                            new_value = coordinate_count - 1
                            self.table.set(self.table.selection(), column='Particles', value=new_value)
                            self.table.set(self.table.selection(), column='Updated', value='Yes')
                            indexesToDelete = self.nearCoordinates(index)
                            for idx in indexesToDelete:
                                self.imageCanvas.delete(idx)
                                self.shapes.pop(idx)

                            if index in self.coordinatesDict[self.micName]:
                                self.coordinatesDict[self.micName].pop(index)

                            self.totalCoordinates -= 1
                            self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
                            self.hasChanges[self.micName] = True
                            break
                elif self.filament:  # Filament picking action
                    shape = self.addCoordinate(x * self.scale / self.zoomFactor, y * self.scale / self.zoomFactor)
                    new_value = coordinate_count + 1
                    self.table.set(self.table.selection(), column="Particles", value=new_value)
                    self.table.set(self.table.selection(), column="Updated", value='Yes')
                    self.coordinatesDict[self.micName][shape] = ((x - self.xOffset) * self.scale / self.zoomFactor,
                                                                 (y - self.yOffset) * self.scale / self.zoomFactor)
                    self.totalCoordinates += 1
                    self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
                    self.hasChanges[self.micName] = True

                elif self.drag:  # Move coordinate or drag all image
                    for index, coords in self.shapes.items():
                        distance = np.sqrt((x - self.xOffset - coords[0]) ** 2 + (y - self.yOffset - coords[1]) ** 2)
                        if distance < self.shapeRadius:
                            self.root.config(cursor='hand2')
                            self.moveShape = True
                            indexesToMove = self.nearCoordinates(index)
                            newX = self.root.winfo_pointerx() - self.root.winfo_rootx() - self.coordX
                            newY = self.root.winfo_pointery() - self.root.winfo_rooty() - self.coordY

                            for idx in indexesToMove:
                                self.shapes[idx] = (self.shapes[idx][0] + newX, self.shapes[idx][1] + newY)
                                # Move the shape to a new position
                                self.imageCanvas.move(idx, newX, newY)

                            # Update de coordinates
                            self.coordX = self.root.winfo_pointerx() - self.root.winfo_rootx()
                            self.coordY = self.root.winfo_pointery() - self.root.winfo_rooty()
                            coordXY = self.coordinatesDict[self.micName][index]
                            self.coordinatesDict[self.micName][index] = (coordXY[0] + newX * self.scale / self.zoomFactor,
                                                                         coordXY[1] + newY * self.scale / self.zoomFactor)

                            if self.selectedCoordinate is not None:
                                indexesToPaint = self.nearCoordinates(self.selectedCoordinate)
                                for idx in indexesToPaint:
                                    self.imageCanvas.itemconfigure(idx, outline=self.selectedColor)

                            self.selectedCoordinate = index

                            for idx in indexesToMove:
                                self.imageCanvas.itemconfigure(idx, outline='red')

                            if newX != 0 or newY != 0:
                                self.table.set(self.table.selection(), column="Updated", value='Yes')
                                self.hasChanges[self.micName] = True
                            break

                    if not self.moveShape:  # Drag the image and the shapes
                        self.root.config(cursor='fleur')
                        self.moveShape = False
                        self.onDrag(event)
        else:
            self.root.config(cursor='')

    def onDrag(self, event):
        """Move the image and shapes on drag"""
        x = (event.x - self.drag_data['x'])
        y = (event.y - self.drag_data['y'])
        self.xOffset += x
        self.yOffset += y
        self.imageCanvas.move('all', x, y)
        self.drag_data['x'] = event.x
        self.drag_data['y'] = event.y

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
        scaledImage = self.imagePIL.resize((int(new_width / self.scale), int(new_height / self.scale)))
        self.imageTk = ImageTk.PhotoImage(scaledImage)
        self.imageCanvas.config(scrollregion=self.imageCanvas.bbox("all"))
        self.image = self.imageCanvas.create_image(0, 0, anchor=tk.NW, image=self.imageTk, tags='image')
        self.drawCoordinates(self.micName)
        if self.zoomFactor != 1:
            x, y = ((event.x - self.drag_data['x']) / self.scale * self.zoomFactor,
                    (event.y - self.drag_data['y']) / self.scale * self.zoomFactor)
            self.xOffset -= x
            self.yOffset -= y
            self.imageCanvas.move('all', -x, -y)
            self.drag_data['x'] = event.x
            self.drag_data['y'] = event.y

    def loadMicrograph(self):
        """Load and display the selected micrograph and coordinates"""
        if self.micName is None:
            self.micName = next(iter(self.coordinatesDict.keys()))

        imagePath = os.path.abspath(self.micrographPathDict[self.micName])
        if imagePath:
            try:
                ext = imagePath.split('.')[-1]
                imageReader = ImageReadersRegistry._readers[ext]
                self.imagePIL = imageReader.open(imagePath)
                self.imageSize = self.imagePIL.size
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

                scaledImage = self.imagePIL.resize((int(dpiWidth), int(dpiHeight)))
                self.imageTk = ImageTk.PhotoImage(scaledImage)

                # self.quadtree = Index(bbox=[0, 0, self.imagePIL.size[0], self.imagePIL.size[1]])
                self.imageCanvas.delete("all")
                self.image = self.imageCanvas.create_image(0, 0, anchor=tk.NW, image=self.imageTk, tags='image')
                self.drawCoordinates(self.micName)

            except Exception as e:
                logger.error(f"Error loading image '{self.micName}': {e}")
        else:
            logger.error("Unable to upload the file %s. Make sure the path is correct." % imagePath)

    def applyPowerHistogram(self):
        """ Create the histogram and the plot"""
        self.histWindow = tk.Toplevel(self.root)
        self.histWindow.title('Power histogram')
        figure, axes = plt.subplots(figsize=(6, 4))
        self.plotHistogram(axes)
        histCanvas = FigureCanvasTkAgg(figure, master=self.histWindow)
        histCanvasWidget = histCanvas.get_tk_widget()
        histCanvasWidget.grid(row=0, column=0, sticky="news", padx=10, pady=10)
        # Create the slider
        self.powerSlider1 = ttk.Scale(self.histWindow, from_=0, to=255, length=410,
                                      orient=tk.HORIZONTAL, command=self.filterCoordinates)
        self.powerSlider1.grid(row=1, column=0, sticky="ns", padx=5, pady=5)
        self.powerSlider2 = ttk.Scale(self.histWindow, from_=0, to=255, length=410,
                                      orient=tk.HORIZONTAL, command=self.filterCoordinates)
        self.powerSlider2.grid(row=2, column=0, sticky="ns", padx=5, pady=5)
        self.powerSlider2.set(self.powerSlider2['to'])

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
        self.loadMicrograph()
        self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
        self.table.set(self.table.selection(), column="Particles", value=len(self.coordinatesDict[self.micName]))

    def histWindowSaveClose(self):
        """Save the new coordinates taking into account the sliders selected pixel range """
        pass
        # Updating current micrograph
        # for i in range(2, len(self.shapes)):
        #     currentState = self.imageCanvas.itemcget(i, "state")
        #     if currentState == 'hidden':
        #         self.shapes.pop(i)
        #         self.imageCanvas.delete(i)
        #         self.coordinatesDict[self.micName].pop(i)
        #         self.totalCoordinates -= 1

        # Updating the rest of micrographs
        # for micName, coords in self.coordinatesDict.items():
        #     if micName != self.micName:
        #         pass

        # self.histWindowClose()

    def filterCoordinates(self, value):
        slider1Value = self.powerSlider1.get()
        slider2Value = self.powerSlider2.get()
        if slider1Value + 20 >= slider2Value:
            slider2Value += 1
            if slider2Value < 255:
                self.powerSlider2.set(slider2Value)
        if slider2Value - 20 <= slider1Value:
            slider1Value -= 1
            if slider1Value > 0:
                self.powerSlider1.set(slider1Value)

        self.removeCoordinates(slider1Value, slider2Value)

    def removeCoordinates(self, value1, value2):
        pixelRange = (round(float(value1)), round(float(value2)))
        for index, coords in self.shapes.items():
            pixelValue = self.getPixelValue(coords[0], coords[1])
            new_state = "normal"
            if int(pixelValue) < pixelRange[0] or int(pixelValue) > pixelRange[1]:
                new_state = "hidden"
            self.imageCanvas.itemconfigure(index, state=new_state)

    def plotHistogram(self, ejes):
        """Plot the micrograph histogram"""
        img_array = np.array(self.imagePIL)
        # Apply the power transformation
        img_power = np.power(img_array / 255.0, 2.0) * 255.0
        hist_power, bins_power = np.histogram(img_power.flatten(), bins=256, range=[0, 256])

        # Plot histogram
        ejes.clear()
        ejes.plot(hist_power)
        ejes.set_xlabel("Pixel Value")
        ejes.set_ylabel("Frequency")

    def resetMicrograph(self):
        """Remove all coordinates from current micrograph"""
        result = messagebox.askquestion("Confirmation", "Are you sure you want to reset the micrograph?",
                                        icon='warning', **{'parent': self.root})
        if result == 'yes':
            self.totalCoordinates -= len(self.coordinatesDict[self.micName])
            self.coordinatesDict[self.micName] = {}
            self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
            self.table.set(self.table.selection(), column="Particles", value=0)
            self.imageCanvas.delete("shape")
            self.table.set(self.table.selection(), column='Updated', value='Yes')
            self.hasChanges[self.micName] = True

    def restoreMicrograph(self):
        """Restore all removed coordinates from current micrograph"""
        result = messagebox.askquestion("Confirmation", "Are you sure you want to restore the micrograph?",
                                        icon='warning', **{'parent': self.root})
        if result == 'yes':
            self.totalCoordinates -= len(self.coordinatesDict[self.micName])
            self.coordinatesDict[self.micName] = self.oldCoordinatesDict[self.micName]
            self.totalCoordinates += len(self.coordinatesDict[self.micName])
            self.totalPickButton.configure(text=f"Total picks: {self.totalCoordinates}")
            self.imageCanvas.delete("shape")
            self.table.set(self.table.selection(), column="Particles", value=len(self.coordinatesDict[self.micName]))
            self.table.set(self.table.selection(), column='Updated', value='No')
            self.drawCoordinates(self.micName)
            if self.micName in self.hasChanges:
                self.hasChanges.pop(self.micName)

    def onSliderMapMove(self, value):
        self.infoLabel.config(text=f"{float(value):.0f}")

    def drawCoordinates(self, micName):
        """Draw the coordinates over the micrograph"""
        self.shapes = {}
        coordinates = self.coordinatesDict[micName]
        self.shapeRadius = self.boxSize / self.scale / 2 * self.zoomFactor
        self.auxCoordinatesDict = dict()
        for index, coord in coordinates.items():
            self.addCoordinate(coord[0], coord[1])
        if self.auxCoordinatesDict:
            self.coordinatesDict[micName] = self.auxCoordinatesDict

    def addCoordinate(self, x, y):
        """Create a coordinate"""
        xTrans, yTrans = x / self.scale * self.zoomFactor, y / self.scale * self.zoomFactor
        shape = None
        if self.drawCircles:
            circle = self.imageCanvas.create_oval(xTrans - self.shapeRadius, yTrans - self.shapeRadius,  xTrans + self.shapeRadius,
                                                  yTrans + self.shapeRadius, outline=self.selectedColor, width=1, fill="",
                                                  tags='shape')
            self.shapes[circle] = (xTrans, yTrans)
            self.auxCoordinatesDict[circle] = (x, y)
            shape = circle
            # self.quadtree.insert(circle, (x, y))

        if self.drawSquares:
            square = self.imageCanvas.create_rectangle(xTrans - self.shapeRadius, yTrans - self.shapeRadius,  xTrans + self.shapeRadius,
                                                       yTrans + self.shapeRadius, outline=self.selectedColor, width=1, fill="",
                                                       tags='shape')
            self.shapes[square] = (xTrans, yTrans)
            if not self.drawCircles:
                self.auxCoordinatesDict[square] = (x, y)
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
        self.loadMicrograph()

    def applyGaussianBlur(self):
        self.gaussianBlurVar.set(not self.gaussianBlurVar.get())
        if self.gaussianBlurVar.get():
            self.filterMenu.entryconfigure(1, image=getImage(Icon.CHECKED))
        else:
            self.filterMenu.entryconfigure(1, image=getImage(Icon.UNCHECKED))
        self.loadMicrograph()

    def applyInvertContrast(self):
        self.invertContrastVar.set(not self.invertContrastVar.get())
        if self.invertContrastVar.get():
            self.filterMenu.entryconfigure(2, image=getImage(Icon.CHECKED))
        else:
            self.filterMenu.entryconfigure(2, image=getImage(Icon.UNCHECKED))
        self.loadMicrograph()

    def applyBandpassFilter(self):
        self.bandpassFilterVar = not self.bandpassFilterVar
        self.loadMicrograph()


class CoordinateView(View):

    def __init__(self, root, emSet: EMSet, protocol):
        self.root = root
        self._emSet = emSet
        self.protocol = protocol

    def show(self):
        root = tk.Toplevel(self.root)
        app = MainWindow(root, self._emSet)
        root.geometry("1315x915")
        root.mainloop()


class CoordinateViewer(Viewer):
    _name = 'Scipion coordinates viewer'
    _targets = [SetOfCoordinates]

    def _visualize(self, obj, **kwargs):
        return [CoordinateView(self.getTkRoot(), obj, self.protocol)]

