# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
"""
This module implement some base classes and utils for wizards
The content of this module is not discovered at runtime by pyworkflow.
Usage of this content is though importing it
"""

import tkinter as tk
import tkinter.ttk as ttk

import pyworkflow as pw
import pyworkflow.object as pwobj
import pyworkflow.wizard as pwizard
import pyworkflow.gui.dialog as dialog
from matplotlib import cm
from pwem import convertPixToLength, splitRange
from pyworkflow.gui.plotter import getHexColorList
from pyworkflow.gui.tree import BoundTree, ListTreeProvider, ListTreeProviderString, AttributesTreeProvider
from pyworkflow.gui.widgets import LabelSlider, ExplanationText

import pwem.constants as emcts
import pwem.objects as emobj
from pwem import emlib
from pyworkflow.protocol import IntParam, StringParam, FloatParam, LEVEL_ADVANCED

# Color map wizard constants
HIGHEST_ATTR = 'highest'
LOWEST_ATTR = 'lowest'
INTERVAL_ATTR = 'intervals'
COLORMAP_ATTR = 'colorMap'


# ===============================================================================
#    Wizard EM base class
# ===============================================================================

class EmWizard(pwizard.Wizard):

    def _getMics(self, objs):
        return [mic.clone() for mic in objs]

    def _getParticles(self, objs, num=100):
        particles = []
        for i, par in enumerate(objs):
            # Cloning the particle
            particle = par.clone()
            if i == num:  # Only load up to NUM particles
                break
            particles.append(particle)

        return particles

    def _getVols(self, objs):
        vols = []
        if isinstance(objs, emobj.Volume):
            vols.append(objs)
        else:
            vols = [vol.clone() for vol in objs]

        return vols

    def _getListProvider(self, objs):
        """ This should be implemented to return the list
        of object to be displayed in the tree.
        """
        provider = None
        if objs.hasValue():
            # If objs is a PointerList currently it can only be formed of SetOfVolumes and Volume
            # (for protocol align_volume). Should this change review this part
            if isinstance(objs, pwobj.PointerList):
                vols_total = []
                for pointer in objs:
                    obj = pointer.get()
                    vols = self._getVols(obj)
                    vols_total.extend(vols)
                provider = ListTreeProvider(vols_total)
            else:
                objs = objs.get()

                if isinstance(objs, emobj.SetOfMicrographs):
                    mics = self._getMics(objs)
                    provider = ListTreeProvider(mics)

                if isinstance(objs, emobj.SetOfParticles):
                    particles = self._getParticles(objs)
                    provider = ListTreeProvider(particles)

                if isinstance(objs, emobj.SetOfVolumes) or isinstance(objs,
                                                                      emobj.Volume):
                    vols = self._getVols(objs)
                    provider = ListTreeProvider(vols)

            return provider
        return None

    def _getInputProtocol(self, targets, protocol):
        label = []
        value = []

        for k, v in targets:
            if k.__name__ == protocol.getClassName():
                label = v
                for val in v:
                    value.append(protocol.getAttributeValue(val))

        if len(label) > 1:
            return label, value
        else:
            return label[0], value[0]


# ===============================================================================
#    Wizards base classes
# ===============================================================================

class DownsampleWizard(EmWizard):

    def show(self, form, value, label, units=emcts.UNIT_PIXEL):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            args = {'unit': units,
                    'downsample': value
                    }

            d = DownsampleDialog(form.root, provider, **args)
            if d.resultYes():
                form.setVar(label, d.getDownsample())
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)

    @classmethod
    def getView(self):
        return "wiz_downsampling"


class CtfWizard(EmWizard):

    def show(self, form, value, label, units=emcts.UNIT_PIXEL):
        protocol = form.protocol
        #        form.setParamFromVar('inputMicrographs') # update selected input micrographs
        provider = self._getProvider(protocol)

        if provider is not None:
            args = {'unit': units,
                    'lf': value[0],
                    'hf': value[1]
                    }
            d = CtfDialog(form.root, provider, **args)

            if d.resultYes():
                form.setVar(label[0], d.getLowFreq())
                form.setVar(label[1], d.getHighFreq())
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)

    @classmethod
    def getView(self):
        return "wiz_ctf"


class MaskRadiusWizard(EmWizard):

    def show(self, form, value, label, units=emcts.UNIT_PIXEL):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            d = MaskPreviewDialog(form.root,
                                  provider,
                                  maskRadius=value,
                                  unit=units)
            if d.resultYes():
                if units == emcts.UNIT_ANGSTROM:
                    value = round(d.getRadiusAngstroms(d.radiusSlider))  # Must be an integer
                else:
                    value = d.getRadius(d.radiusSlider)
                self.setVar(form, label, value)
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)

    def setVar(self, form, label, value):
        form.setVar(label, value)


class MaskRadiiWizard(EmWizard):

    def show(self, form, value, label, units=emcts.UNIT_PIXEL):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            d = MaskRadiiPreviewDialog(form.root,
                                       provider,
                                       innerRadius=value[0],
                                       outerRadius=value[1],
                                       unit=units)
            if d.resultYes():
                if units == emcts.UNIT_ANGSTROM:
                    form.setVar(label[0], d.getRadiusAngstroms(d.radiusSliderIn))
                    form.setVar(label[1], d.getRadiusAngstroms(d.radiusSliderOut))
                else:
                    form.setVar(label[0], d.getRadius(d.radiusSliderIn))
                    form.setVar(label[1], d.getRadius(d.radiusSliderOut))
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)

    @classmethod
    def getView(self):
        return "wiz_volume_mask_radii"


class ParticleMaskRadiusWizard(MaskRadiusWizard):
    pass


class VolumeMaskRadiusWizard(MaskRadiusWizard):
    pass


class ParticlesMaskRadiiWizard(MaskRadiiWizard):
    pass


class VolumeMaskRadiiWizard(MaskRadiiWizard):
    pass


class FilterWizard(EmWizard):

    def show(self, form, value, label, mode, unit=emcts.UNIT_PIXEL, **args):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            self.mode = mode
            args.update({'mode': self.mode,
                         'unit': unit})
            if self.mode == emcts.FILTER_LOW_PASS:
                args['showLowFreq'] = False
                args['highFreq'] = value[1]
                args['freqDecay'] = value[2]
            elif self.mode == emcts.FILTER_HIGH_PASS:
                args['showHighFreq'] = False
                args['lowFreq'] = value[0]
                args['freqDecay'] = value[2]
            elif self.mode == emcts.FILTER_BAND_PASS:
                args['lowFreq'] = value[0]
                args['highFreq'] = value[1]
                args['freqDecay'] = value[2]
            else:
                raise Exception("Unknown mode '%s'" % self.mode)

            d = BandPassFilterDialog(form.root, provider, **args)

            if d.resultYes():
                def setFormValue(flag, value, index):
                    if args.get(flag, True):
                        if unit == emcts.UNIT_ANGSTROM:
                            value = d.samplingRate / value
                        form.setVar(label[index], value)

                setFormValue('showLowFreq', d.getLowFreq(), 0)
                setFormValue('showHighFreq', d.getHighFreq(), 1)
                setFormValue('showDecay', d.getFreqDecay(), 2)

        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)


class FilterParticlesWizard(FilterWizard):
    pass


class FilterVolumesWizard(FilterWizard):
    pass


class GaussianWizard(EmWizard):

    def show(self, form, value, label, units=emcts.UNIT_PIXEL_FOURIER):
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            args = {'freqSigma': value,
                    'unit': units
                    }

            d = GaussianFilterDialog(form.root, provider, **args)
            if d.resultYes():
                form.setVar(label, d.getFreqSigma())
        else:
            dialog.showWarning("Empty input", "Select elements first", form.root)


class GaussianParticlesWizard(GaussianWizard):
    pass


class GaussianVolumesWizard(GaussianWizard):
    pass


class PDBVolumeWizard(EmWizard):
    def show(self, form):
        if form.protocol.inputPDB.hasValue():
            pdb = form.protocol.inputPDB.get()
            print("pdb ", pdb._volume, str(pdb._volume))
            if pdb._volume:
                print("Setting ", str(pdb.getVolume()))
                ptr = pwobj.Pointer()
                ptr.copy(form.protocol.inputPDB)
                ptr.setExtended(ptr.getExtended() + "._volume")
                #                 ptr.set(pdb.getVolume())
                form.setVar('inputVol', ptr)


class ColorScaleWizardBase(EmWizard):
    """ Base wizard to edit color scale parameters
    Usage:
        1.- define new class inheriting this one
        2.- declare _targets = ColorScaleWizardBase.defineTargets(youviewer) in a class scope, right after class definition.
        3.- call ColorScaleWizardBase.defineColorScaleParams(group) in your viewer._defineParams
        4.- use attributes in your plotting method
    """

    @classmethod
    def defineTargets(cls, *viewersClass):
        """ :return targets list per each viewer class passed"""
        targets = []
        for viewer in viewersClass:
            targets.append((viewer, [HIGHEST_ATTR, LOWEST_ATTR, INTERVAL_ATTR, COLORMAP_ATTR]))

        return targets

    def show(self, form):
        viewer = form.protocol

        d = ColorScaleDialog(form.root, viewer.lowest.get(), viewer.highest.get(), viewer.intervals.get(),
                             viewer.colorMap.get())

        # If accepted
        if d.resultYes():
            form.setVar(LOWEST_ATTR, d.getLowest())
            form.setVar(HIGHEST_ATTR, d.getHighest())
            form.setVar(INTERVAL_ATTR, d.getIntervals())
            form.setVar(COLORMAP_ATTR, d.getColorPalette())

    @staticmethod
    def defineColorScaleParams(form, defaultHighest=10, defaultLowest=0, defaultIntervals=11, defaultColorMap="jet"):

        line = form.addLine("Color scale options:",
                            help="Options to define the color scale limits, intervals (advanced) and color set. Useful when you have outliers ruining your "
                                 "visualization/plot.")
        line.addParam(HIGHEST_ATTR, FloatParam, default=defaultHighest,
                      label="Highest",
                      help="Highest value for the scale")

        line.addParam(LOWEST_ATTR, FloatParam, default=defaultLowest,
                      label="Lowest",
                      help="lowest value of the scale.")

        line.addParam(INTERVAL_ATTR, IntParam, default=defaultIntervals,
                      label="Intervals",
                      help="Number of labels of the scale",
                      expertLevel=LEVEL_ADVANCED)

        line.addParam(COLORMAP_ATTR, StringParam, default=defaultColorMap,
                      label="Color set",
                      help="Combination of color for the scale")


# ===============================================================================
#  Dialogs used by wizards
# ===============================================================================

class PreviewDialog(dialog.Dialog):
    """ This will be the base class for several wizards.
    The layout of this wizard will be:
    1. Left panel(Items) that contains a list of items to preview
    2. Right-top panel (Preview) where some preview of the items will be displayed
    3. Right-bottom panel (Controls) where some controls can change the preview
    """

    def __init__(self, parent, provider, **args):
        """
        Params:
            parent: parent windows of the dialog.
            provider: the TreeProvider to populate items tree.
        """
        # Set the attributes in **args
        for k, v in args.items():
            setattr(self, k, v)

        self.provider = provider
        self.firstItem = provider.getObjects()[0]
        self.isInnerRad = False
        self.isMakingBigger = False
        self.step = 1
        self.expText = []
        buttons = [('Select', dialog.RESULT_YES),
                   ('Cancel', dialog.RESULT_CANCEL)]
        dialog.Dialog.__init__(self, parent, "Wizard",
                               buttons=buttons, default='Select', **args)

    def body(self, bodyFrame):
        bodyFrame.config()
        bodyFrame.columnconfigure(0, weight=1)
        bodyFrame.columnconfigure(1, weight=1)
        bodyFrame.columnconfigure(2, weight=1)

        # Create explanation label
        self.expText = ExplanationText(bodyFrame)
        self.expText.text.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky='news')

        # Create items frame
        itemsFrame = tk.Frame(bodyFrame, bg=pw.TK_GRAY_DEFAULT)
        itemsFrame.grid(row=1, column=0, padx=5, sticky='news')
        itemsFrame.columnconfigure(0, weight=1)
        itemsFrame.rowconfigure(0, weight=1)
        itemsTree = BoundTree(itemsFrame, self.provider)
        itemsTree.grid(row=0, column=0, padx=5, pady=5, sticky='news')
        itemsTree.itemClick = self._itemSelected

        # Create preview frame
        previewFrame = tk.Frame(bodyFrame, bg=pw.TK_GRAY_DEFAULT)
        previewFrame.grid(row=1, column=1, padx=5, pady=5)
        self._beforePreview()
        self._createPreview(previewFrame)

        # Create controls frame
        controlsFrame = tk.Frame(bodyFrame)
        controlsFrame.grid(row=2, column=0, columnspan=2, padx=5, pady=5, sticky='sew')
        controlsFrame.columnconfigure(0, weight=1)
        controlsFrame.rowconfigure(2, weight=1)

        self._createControls(controlsFrame)
        self._itemSelected(self.firstItem)
        itemsTree.selectChildByIndex(0)  # Select the first item

    def _beforePreview(self):
        """ Called just before setting the preview.
        This is the place to set data values such as: labels, constants...
        """
        pass

    def _createPreview(self, frame):
        """ Should be implemented by subclasses to
        create the items preview.
        """
        pass

    def _createControls(self, frame):
        """ Create controls to be used. """
        pass

    def _itemSelected(self, obj):
        """ This will be call when an item in the tree is selected. """
        pass


class ImagePreviewDialog(PreviewDialog):

    def _beforePreview(self):
        self.dim = 256
        self.previewLabel = ''

    def _createPreview(self, frame):
        """ Should be implemented by subclasses to
        create the items preview.
        """
        from pyworkflow.gui.matplotlib_image import ImagePreview
        self.preview = ImagePreview(frame, self.dim, label=self.previewLabel)
        self.preview.grid(row=0, column=0)

    def _itemSelected(self, obj):

        index = obj.getIndex()
        filename = emlib.image.ImageHandler.fixXmippVolumeFileName(obj)
        if index:
            filename = "%03d@%s" % (index, filename)

        self.image = emlib.image.ImageHandler()._img

        try:
            self.image.readPreview(filename, self.dim)
            if filename.endswith('.psd'):
                self.image.convertPSD()
            self.Z = self.image.getData()
        except Exception as e:
            from pyworkflow.gui.matplotlib_image import getPngData
            self.Z = getPngData(pw.findResource('no-image.gif'))
            dialog.showError("Input particles", "Error reading image <%s>"
                             % filename, self)
        self.preview.updateData(self.Z)


class DownsampleDialog(ImagePreviewDialog):

    def _beforePreview(self):
        self.expText.updateExpText(emcts.FREQ_BANDPASS_WIZ_MSG, width=75)
        ImagePreviewDialog._beforePreview(self)
        self.lastObj = None
        self.rightPreviewLabel = "PSD"
        self.message = "Computing PSD..."
        self.previewLabel = "Micrograph"
        self.rightImage = emlib.image.ImageHandler()._img

    def _createPreview(self, frame):
        """ Should be implemented by subclasses to
        create the items preview.
        """
        leftFrame = tk.Frame(frame)
        leftFrame.grid(row=0, column=0)

        rightFrame = tk.Frame(frame)
        rightFrame.grid(row=0, column=1)

        ImagePreviewDialog._createPreview(self, leftFrame)
        self.rightPreview = self._createRightPreview(rightFrame)
        self.rightPreview.grid(row=0, column=0)

    def _createRightPreview(self, rightFrame):
        from pyworkflow.gui.matplotlib_image import ImagePreview
        return ImagePreview(rightFrame, self.dim, label=self.rightPreviewLabel)

    def _createControls(self, frame):
        self.downVar = tk.StringVar()
        self.downVar.set(getattr(self, 'downsample', 1))
        downFrame = tk.Frame(frame)
        downFrame.grid(row=0, column=0, pady=5, sticky='new')
        downLabel = tk.Label(downFrame, text='Downsample')
        downLabel.grid(row=0, column=0, sticky='ew')
        downEntry = tk.Entry(downFrame, width=10, textvariable=self.downVar)
        downEntry.grid(row=0, column=1, padx=5, sticky='ew')
        downButton = tk.Button(downFrame, text='Preview', command=self._doPreview)
        downButton.grid(row=1, column=1, padx=5, pady=5, sticky='ew')

    def getDownsample(self):
        return float(self.downVar.get())

    def manageMaskVals(self):
        # To be coded by each child, if necessary
        pass

    def _itemSelected(self, obj):
        self.lastObj = obj
        ImagePreviewDialog._itemSelected(self, obj)

        dialog.FlashMessage(self, self.message, func=self._computeRightPreview)
        self.rightPreview.updateData(self.rightImage.getData())
        self.manageMaskVals()

    def _doPreview(self, e=None):
        if self.lastObj is None:
            dialog.showError("Empty selection",
                             "Select an item first before preview", self)
        else:
            self._itemSelected(self.lastObj)

    def _computeRightPreview(self):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        emlib.fastEstimateEnhancedPSD(self.rightImage,
                                      self.lastObj.getFileName(),
                                      self.getDownsample(), self.dim, 2)


class CtfDialog(DownsampleDialog):

    def _createRightPreview(self, rightFrame):
        self.isInnerRad = True
        listeners = {"<Button-4>": self.makeBigger, "<Button-5>": self.makeSmaller,
                     "<Up>": self.upKeyPress, "<Down>": self.downKeyPress}
        from pyworkflow.gui.matplotlib_image import PsdPreview
        return PsdPreview(rightFrame, self.dim, self.lf, self.hf,
                          label=self.rightPreviewLabel, listenersDict=listeners)

    def downKeyPress(self, event):
        self.isInnerRad = False
        self.highlightOuterSlider()

    def upKeyPress(self, event):
        self.isInnerRad = True
        self.highlightInnerSlider()

    def highlightOuterSlider(self):
        self.hfSlider.highlightLabel()
        self.lfSlider.removeHighlightFromLabel()

    def highlightInnerSlider(self):
        self.lfSlider.highlightLabel()
        self.hfSlider.removeHighlightFromLabel()

    def makeBigger(self, event):
        self.isMakingBigger = True
        if self.isInnerRad:
            self.lf = self.lf + self.step
        else:
            new_val = self.hf + self.step
            if new_val <= 0.5:  # Don't make the mask bigger unless the is equal or lower than the max
                self.hf = new_val
        self.manageMaskVals()

    def makeSmaller(self, event):
        self.isMakingBigger = False
        if self.isInnerRad:
            new_val = self.lf - self.step
            if new_val >= 0:
                self.lf = new_val
        else:
            self.hf = self.hf - self.step
        self.manageMaskVals()

    def manageMaskVals(self):
        if self.isMakingBigger:
            if self.isInnerRad and self.lf >= self.hf:  # Inner ring can't be bigger than outer ring
                # Subtract one step to go back to the nearest lower value
                self.lf = self.hf - self.step
        else:
            if not self.isInnerRad and self.hf <= self.lf:  # Outer ring can't be smaller than inner
                # ring
                # Add one step to go back to the nearest higher value
                self.hf = self.lf + self.step

        # Set the ring sliders in case it comes from the mouse wheel
        self.setFreq(self.lfSlider, self.lf)
        self.setFreq(self.hfSlider, self.hf)

        # Show values
        self.showValues(self.hfVar, self.hfSlider)
        self.showValues(self.lfVar, self.lfSlider)

        # Update mask
        self.rightPreview.updateFreq(self.lf, self.hf)

    def _showInAngstroms(self):
        return getattr(self, 'showInAngstroms', False)

    def _createControls(self, frame):
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(2, weight=1)

        self.step = 0.01
        self.freqFrame = ttk.LabelFrame(frame, text="Frequencies")
        self.freqFrame.columnconfigure(0, weight=3)
        self.freqFrame.grid(row=0, column=2, sticky='news')
        self.lfSlider = LabelSlider(self.freqFrame, 'Low freq',
                                    value=self.lf,
                                    from_=0.01,
                                    to=0.5,
                                    length=178,
                                    labelWidth=13,
                                    step=self.step,
                                    showvalue=0,
                                    callback=lambda a, b, c: self.updateSliderInnerRadius(),
                                    )
        self.lfSlider.grid(row=0, column=0, padx=5, pady=5, sticky='news')
        self.lfSlider.highlightLabel()
        self.hfSlider = LabelSlider(self.freqFrame, 'High freq',
                                    value=self.lf,
                                    from_=0.01,
                                    to=0.5,
                                    length=178,
                                    labelWidth=13,
                                    step=self.step,
                                    showvalue=0,
                                    callback=lambda a, b, c: self.updateSliderOuterRadius(),
                                    )
        self.hfSlider.grid(row=1, column=0, padx=5, pady=5, sticky='news')

        # Pack and configure low freq slider
        self.lfVar = tk.StringVar()
        self.lfLabel = tk.Label(self.freqFrame, textvariable=self.lfVar, width=26)
        self.lfLabel.grid(row=0, column=1, sticky='nse', padx=5, pady=5)
        # Pack and configure high freq slider
        self.hfVar = tk.StringVar()
        self.hfLabel = tk.Label(self.freqFrame, textvariable=self.hfVar, width=26)
        self.hfLabel.grid(row=1, column=1, sticky='nse', padx=5, pady=5)
        # Update both mask and sliders with the initial values
        self.manageMaskVals()

    def getDownsample(self):
        return 1.0  # Micrograph previously downsample, not taken into account here

    def updateFreqRing(self):
        self.rightPreview.updateFreq(self.getLowFreq(), self.getHighFreq())
        self.showValues(self.lfVar, self.lfSlider)
        self.showValues(self.hfVar, self.hfSlider)

    def updateSliderOuterRadius(self):
        new_val = self.getFreq(self.hfSlider)
        # Carry out the action only if the slider is clicked & dragged (prevent from minimal variations when the mouse
        # is on the slider but not clicking on it)
        if abs(new_val - self.hf) >= self.step:
            self.highlightOuterSlider()
            self.isInnerRad = False
            self.isMakingBigger = False
            # Check if the user is making the ring bigger
            if new_val > self.hf:
                self.isMakingBigger = True

            self.hf = new_val
            self.manageMaskVals()

    def updateSliderInnerRadius(self):
        new_val = self.getFreq(self.lfSlider)
        # Highlight only if the slider is clicked & dragged (prevent from minimal variations when the mouse is on the
        # slider but not clicking on it)
        if abs(new_val - self.lf) >= self.step:
            self.highlightInnerSlider()
            self.isInnerRad = True
            self.isMakingBigger = False
            # Check if the user is making the ring bigger
            if new_val > self.lf:
                self.isMakingBigger = True

            self.lf = new_val
            self.manageMaskVals()

    def showValues(self, var2set, labSlider):
        """
        Show the values selected for the inner and outer radius. If the units are angstroms (sampling_rate = 1,
        it will show only one value to avoid redundancies
        """
        sr = self.firstItem.getSamplingRate()
        freqVal = max(self.getFreq(labSlider), 0.0001)  # avoid division by zero
        if sr == 1:
            var2set.set('{:2.2f} {}'.format(labSlider.slider.get(), emcts.UNIT_PIXEL))
        else:
            var2set.set('{:2.2f} rad/Å | {:5.1f} Å'.format(labSlider.slider.get(),
                                                           self.getDownsample() * sr / freqVal))

    def getLowFreq(self):
        return self.lfSlider.get()

    def getHighFreq(self):
        return self.hfSlider.get()

    @staticmethod
    def getFreq(freqSlider):
        return freqSlider.get()

    @staticmethod
    def setFreq(freqSlider, val):
        freqSlider.slider.set(val)


class CtfDownsampleDialog(CtfDialog):

    def _createControls(self, frame):
        DownsampleDialog._createControls(self, frame)
        CtfDialog._createControls(self, frame)

    def getDownsample(self):
        return float(self.downVar.get())


class BandPassFilterDialog(DownsampleDialog):

    def _beforePreview(self):
        ImagePreviewDialog._beforePreview(self)
        self.lastObj = None
        self.rightPreviewLabel = "Filtered"
        self.message = "Computing filtered image..."
        self.previewLabel = "Image"
        self.rightImage = emlib.image.ImageHandler()._img

    def _createControls(self, frame):
        self.freqFrame = ttk.LabelFrame(frame,
                                        text="Frequencies (%s)" % self.unit,
                                        padding="5 5 5 5")
        self.freqFrame.grid(row=0, column=0)

        self.showLowFreq = getattr(self, 'showLowFreq', True)
        self.showHighFreq = getattr(self, 'showHighFreq', True)
        self.showDecay = getattr(self, 'showDecay', True)

        if (not self.showLowFreq) or (not self.showHighFreq):
            label_high = 'Freq'
            label_low = 'Freq'
        else:
            label_high = 'High freq'
            label_low = 'Low freq'

        self.samplingRate = 1.0
        self.sliTo = 0.5
        self.sliFrom = 0.01
        if self.unit == emcts.UNIT_ANGSTROM:
            self.samplingRate = self.firstItem.getSamplingRate()
            self.itemDim, _, _ = self.firstItem.getDim()
            self.sliFrom = 2. * self.samplingRate
            self.sliTo = 2. * self.itemDim * self.samplingRate

        self.step = self.sliTo / 1000

        if self.showLowFreq:
            self.lfSlider = self.addFreqSlider(label_low, self.lowFreq, col=0)
        if self.showHighFreq:
            self.hfSlider = self.addFreqSlider(label_high, self.highFreq, col=1)
        if self.showDecay:
            self.freqDecaySlider = self.addFreqSlider('Decay',
                                                      self.freqDecay, col=2)

    def addFreqSlider(self, label, value, col):
        fromValue = self.sliFrom
        toValue = self.sliTo
        if self.unit == emcts.UNIT_ANGSTROM:
            fromValue = self.sliTo
            toValue = self.sliFrom
        slider = LabelSlider(self.freqFrame, label, from_=fromValue, to=toValue,
                             step=self.step, value=value,
                             callback=lambda a, b, c: self.updateFilteredImage())
        slider.grid(row=0, column=col, padx=5, pady=5)
        return slider

    def updateFilteredImage(self):
        self._computeRightPreview()
        self.rightPreview.updateData(self.rightImage.getData())

    def _computeRightPreview(self):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        emlib.bandPassFilter(self.rightImage,
                             emlib.image.ImageHandler.locationToXmipp(self.lastObj),
                             self.getLowFreq(), self.getHighFreq(),
                             self.getFreqDecay(), self.dim)

    def getLowFreq(self):
        if self.showLowFreq:
            if self.unit == emcts.UNIT_ANGSTROM:
                return self.samplingRate / self.lfSlider.get()
            else:
                return self.lfSlider.get()
        return 0.01

    def getHighFreq(self):
        if self.showHighFreq:
            if self.unit == emcts.UNIT_ANGSTROM:
                return self.samplingRate / self.hfSlider.get()
            else:
                return self.hfSlider.get()
        return 1.0

    def getFreqDecay(self):
        if self.showDecay:
            if self.unit == emcts.UNIT_ANGSTROM:
                return self.samplingRate / self.freqDecaySlider.get()
            else:
                return self.freqDecaySlider.get()
        return 0.01


class GaussianFilterDialog(BandPassFilterDialog):

    def _createControls(self, frame):
        self.freqVar = tk.StringVar()
        self.freqVar.set(getattr(self, 'freqSigma', 1))
        freqFrame = tk.Frame(frame)
        freqFrame.grid(row=0, column=0, sticky='nw')
        downEntry = tk.Entry(freqFrame, width=10, textvariable=self.freqVar)
        downEntry.grid(row=0, column=1, padx=5, pady=5)
        downButton = tk.Button(freqFrame, text='Preview',
                               command=self._doPreview)
        downButton.grid(row=0, column=2, padx=5, pady=5)

    def getFreqSigma(self):
        return float(self.freqVar.get())

    def _computeRightPreview(self):
        """ This function should compute the right preview
        using the self.lastObj that was selected
        """
        emlib.gaussianFilter(self.rightImage,
                             emlib.image.ImageHandler.locationToXmipp(self.lastObj),
                             self.getFreqSigma(), self.dim)


class MaskPreviewDialog(ImagePreviewDialog):

    def _beforePreview(self):
        self.dim = 256
        self.unit = getattr(self, 'unit', emcts.UNIT_PIXEL)
        self.samplingRate = self.firstItem.getSamplingRate()
        self.dim_par = self.firstItem.getDim()[0]
        self.ratio = self.dim / float(self.dim_par)
        self.previewLabel = 'Central slice'

    def _createPreview(self, frame):
        """ Should be implemented by subclasses to
        create the items preview.
        """

        # Insert the corresponding explanation text
        self.expText.updateExpText(emcts.CIRCLE_MASK_WIZ_MSG)

        from pyworkflow.gui.matplotlib_image import MaskPreview

        if self.maskRadius == -1:
            self.iniRadius = self.dim_par / 2
        else:
            self.iniRadius = self.maskRadius

        if self.unit == emcts.UNIT_ANGSTROM:
            self.iniRadius = self.iniRadius / self.samplingRate

        listeners = {"<Button-4>": self.makeBigger, "<Button-5>": self.makeSmaller}

        self.preview = MaskPreview(frame, self.dim, label=self.previewLabel,
                                   listenersDict=listeners)

        self.preview.grid(row=1, column=0)

    def makeBigger(self, event):
        new_val = self.iniRadius + 1
        if new_val <= self.dim_par / 2:  # Don't make the mask bigger unless the is equal or lower than the max
            self.iniRadius = self.iniRadius + 1
            self.manageMaskVals()

    def makeSmaller(self, event):
        new_val = self.iniRadius - 1
        if new_val >= 1:  # Don't make the mask smaller unless the radius is equal or greater than 1
            self.iniRadius = self.iniRadius - 1
            self.manageMaskVals()

    def manageMaskVals(self):
        # Set the ring slider in case it comes from the mouse wheel
        self.setRadius(self.radiusSlider, self.iniRadius)
        # Show values
        self.showValues(self.hfVar, self.radiusSlider)
        # Update mask
        self.preview.updateMask(self.iniRadius * self.ratio)

    def _createControls(self, frame):
        self.addRadiusBox(frame)
        self.hfVar = tk.StringVar()
        self.hfLabel = tk.Label(frame, textvariable=self.hfVar, width=20)
        self.hfLabel.grid(row=0, column=1, sticky='NSE', padx=5, pady=5)
        self.manageMaskVals()

    def addRadiusBox(self, parent):
        self.radiusSlider = LabelSlider(parent, 'Mask radius',
                                        from_=1, to=int(self.dim_par / 2),
                                        value=self.iniRadius,
                                        step=1,
                                        length=150,
                                        showvalue=0,
                                        callback=lambda a, b, c: self.updateSliderRadius())
        self.radiusSlider.grid(row=0, column=0, sticky='NEWS')

    def updateSliderRadius(self):

        new_val = self.getRadius(self.radiusSlider)
        # Carry out the action only if the slider is clicked & dragged (prevent from minimal variations when the mouse
        # is on the slider but not clicking on it)
        if abs(new_val - self.iniRadius) >= 1:
            self.iniRadius = new_val
            self.manageMaskVals()

    def showValues(self, var2set, radiusSlider):
        """
        Show the values selected for the inner and outer radius. If the units are angstroms (sampling_rate = 1,
        it will show only one value to avoid redundancies
        """
        pixVal = self.getRadius(radiusSlider)
        if self.samplingRate == 1:
            var2set.set('{:6.1f} {}'.format(pixVal, emcts.UNIT_PIXEL))
        else:
            var2set.set('{:5.0f} pix | {:6.1f} Å'.format(pixVal,
                                                          self.getRadiusAngstroms(radiusSlider),
                                                          ))

    @staticmethod
    def getRadius(radiusSlider):
        return radiusSlider.get()

    @staticmethod
    def setRadius(radiusSlider, val):
        radiusSlider.slider.set(val)

    def getRadiusAngstroms(self, radiusSlider):
        return convertPixToLength(self.samplingRate, radiusSlider.get())


class MaskRadiiPreviewDialog(MaskPreviewDialog):

    def _createPreview(self, frame):
        """ Should be implemented by subclasses to
        create the items preview.
        """
        # Insert the corresponding explanation text
        self.expText.updateExpText(emcts.RING_MASK_WIZ_MSG)

        from pyworkflow.gui.matplotlib_image import MaskPreview
        if self.innerRadius is None:
            self.innerRadius = 1
        if self.outerRadius is None or self.outerRadius == -1 or self.outerRadius > self.dim_par / 2:
            self.outerRadius = int(self.dim_par / 2)

        if self.unit == emcts.UNIT_ANGSTROM:
            self.innerRadius = self.innerRadius / self.samplingRate
            self.outerRadius = self.innerRadius / self.samplingRate

        listeners = {"<Button-4>": self.makeBigger, "<Button-5>": self.makeSmaller,
                     "<Up>": self.upKeyPress, "<Down>": self.downKeyPress}

        self.preview = MaskPreview(frame, self.dim,
                                   label=self.previewLabel,
                                   listenersDict=listeners)

        self.preview.grid(row=1, column=0)

    def upKeyPress(self, event):
        self.isInnerRad = False
        self.highlightOuterSlider()

    def downKeyPress(self, event):
        self.isInnerRad = True
        self.highlightInnerSlider()

    def highlightOuterSlider(self):
        self.radiusSliderOut.highlightLabel()
        self.radiusSliderIn.removeHighlightFromLabel()

    def highlightInnerSlider(self):
        self.radiusSliderIn.highlightLabel()
        self.radiusSliderOut.removeHighlightFromLabel()

    def makeBigger(self, event):
        self.isMakingBigger = True
        if self.isInnerRad:
            self.innerRadius = self.innerRadius + self.step
        else:
            new_val = self.outerRadius + self.step
            if new_val <= int(self.dim_par / 2):  # Don't make the mask bigger unless the is equal or lower than the max
                self.outerRadius = new_val
        self.manageMaskVals()

    def makeSmaller(self, event):
        self.isMakingBigger = False
        if self.isInnerRad:
            new_val = self.innerRadius - self.step
            if new_val >= 0:
                self.innerRadius = new_val
        else:
            self.outerRadius = self.outerRadius - self.step
        self.manageMaskVals()

    def manageMaskVals(self):
        if self.isMakingBigger:
            if self.isInnerRad and self.innerRadius >= self.outerRadius:  # Inner ring can't be bigger than outer ring
                # Subtract one step to go back to the nearest lower value
                self.innerRadius = self.outerRadius - self.step
        else:
            if not self.isInnerRad and self.outerRadius <= self.innerRadius:  # Outer ring can't be smaller than inner
                # ring
                # Add one step to go back to the nearest higher value
                self.outerRadius = self.innerRadius + self.step

        # Set the ring sliders in case it comes from the mouse wheel
        self.setRadius(self.radiusSliderIn, self.innerRadius)
        self.setRadius(self.radiusSliderOut, self.outerRadius)

        # Show values
        self.showValues(self.orVar, self.radiusSliderOut)
        self.showValues(self.irVar, self.radiusSliderIn)

        # Update mask
        self.preview.updateMask(self.outerRadius * self.ratio,
                                self.innerRadius * self.ratio)

    def _createControls(self, frame):
        self.step = 1
        to = int(self.dim_par / 2)
        self.radiusSliderOut = LabelSlider(frame, 'Outer radius',
                                           from_=1, to=to,
                                           value=self.outerRadius, step=self.step,
                                           callback=lambda a, b, c: self.updateSliderOuterRadius(),
                                           length=150, showvalue=0)
        self.radiusSliderOut.highlightLabel()
        self.orVar = tk.StringVar()
        self.orLabel = tk.Label(frame, textvariable=self.orVar, width=20)

        self.radiusSliderIn = LabelSlider(frame, 'Inner radius',
                                          from_=1, to=to,
                                          value=self.innerRadius, step=self.step,
                                          callback=lambda a, b, c: self.updateSliderInnerRadius(),
                                          length=150, showvalue=0)
        self.irVar = tk.StringVar()
        self.irLabel = tk.Label(frame, textvariable=self.irVar, width=20)
        # Pack and configure outer radius slider
        self.radiusSliderOut.grid(row=0, column=0, sticky='NEWS')
        self.radiusSliderOut.columnconfigure(0, weight=1)
        self.orLabel.grid(row=0, column=1, sticky='NSE', padx=5, pady=5)
        self.orLabel.columnconfigure(0, weight=1)
        # Pack and configure inner radius slider
        self.radiusSliderIn.grid(row=1, column=0, sticky='NEWS')
        self.radiusSliderIn.columnconfigure(0, weight=1)
        self.irLabel.grid(row=1, column=1, sticky='NSE', padx=5, pady=5)
        self.irLabel.columnconfigure(0, weight=1)
        # Update both mask and sliders with the initial values
        self.manageMaskVals()

    def updateSliderOuterRadius(self):
        new_val = self.getRadius(self.radiusSliderOut)
        # Carry out the action only if the slider is clicked & dragged (prevent from minimal variations when the mouse
        # is on the slider but not clicking on it)
        if abs(new_val - self.outerRadius) >= self.step:
            self.highlightOuterSlider()
            self.isInnerRad = False
            self.isMakingBigger = False
            # Check if the user is making the ring bigger
            if new_val > self.outerRadius:
                self.isMakingBigger = True

            self.outerRadius = new_val
            self.manageMaskVals()

    def updateSliderInnerRadius(self):
        new_val = self.getRadius(self.radiusSliderIn)
        # Highlight only if the slider is clicked & dragged (prevent from minimal variations when the mouse is on the
        # slider but not clicking on it)
        if abs(new_val - self.innerRadius) >= self.step:
            self.highlightInnerSlider()
            self.isInnerRad = True
            self.isMakingBigger = False
            # Check if the user is making the ring bigger
            if new_val > self.innerRadius:
                self.isMakingBigger = True

            self.innerRadius = new_val
            self.manageMaskVals()


class ColorScaleDialog(dialog.Dialog):
    """ This will assist users to choose the color scale and range for
    local resolution viewers
    """

    def __init__(self, parentWindow, lowest, highest, intervals, colorPalette):
        """
            :param highest: highest resolution value for the scale
            :param lowest: lowest resolution value for the scale
            :param intervals: number of labels for the scale
            :param colorPalette: color palette to use
        """
        self.highest = tk.DoubleVar(value=highest)
        self.lowest = tk.DoubleVar(value=lowest)
        self.intervals = tk.IntVar(value=intervals)
        self.colorPalette = tk.StringVar(value=colorPalette)
        self.info = tk.StringVar()

        # GUI attributes initialization
        self.palette = None
        self.params = None

        dialog.Dialog.__init__(self, parentWindow, "Color scale wizard", default="None")

    # Getters
    def getIntervals(self):
        return int(self.intervals.get())

    def getHighest(self):
        return float(self.highest.get())

    def getLowest(self):
        return float(self.lowest.get())

    def getColorPalette(self):
        return self.colorPalette.get()

    ### ----- GUI methods ------ ###
    def body(self, master):
        """ Draws the main frame of the dialog"""
        body = tk.Frame(self)
        body.grid(row=0, column=0, sticky="news")
        body.grid_columnconfigure(1, weight=10)

        # GUI attributes
        self.palette = tk.Frame(body)
        self.palette.grid(row=0, column=0, sticky='nes',
                          padx=5, pady=5)

        # Params
        self.params = tk.Frame(body)
        self.params.grid(row=0, column=1, sticky='news',
                         padx=5, pady=5)
        self.params.bind("<Key>", self._keyPressedOnParams)

        self._drawPalette()
        self._fillParams()

    def _drawLabel(self, color, value, count):
        label = tk.Label(self.palette, text=value)
        label.config(bg=color)
        label.grid(row=count, column=0, sticky='news')

    def _fillParams(self):

        # Add highest scale value
        highValueLabel = tk.Label(self.params, text="Highest:")
        highValueLabel.grid(row=0, column=0, sticky="e")
        self._addEntry(0, 1, self.highest)

        # Add lowest value
        lowValueLabel = tk.Label(self.params, text="Lowest:")
        lowValueLabel.grid(row=1, column=0, sticky="e")
        self._addEntry(1, 1, self.lowest)

        intervalsLabel = tk.Label(self.params, text="Intervals:")
        intervalsLabel.grid(row=2, column=0, sticky="e")
        self._addEntry(2, 1, self.intervals)

        # Palettes
        paletteLabel = tk.Label(self.params, text="Color set:")
        paletteLabel.grid(row=3, column=0, sticky="e")

        availablePalettes = self.getAvailablePalettes()
        opt = ttk.Combobox(self.params, textvariable=self.colorPalette, values=availablePalettes)
        opt.grid(row=3, column=1, sticky="news")
        self.colorPalette.trace("w", self._paletteChanged)

        # Label to store information: number of color of the color map, e.g.
        infoLabel = tk.Label(self.params, textvariable=self.info, text="")
        infoLabel.grid(row=4, column=0, columnspan=2, sticky="news")

    def _addEntry(self, row, column, var):
        newEntry = tk.Entry(self.params, textvariable=var)
        newEntry.grid(row=row, column=column, sticky="news")
        newEntry.lift()
        # Bind events
        newEntry.bind("<Return>", self._paramChanged)
        newEntry.bind("<FocusIn>", self._selectAllText)
        newEntry.bind("<Button-1>", self._selectAllText)
        newEntry.bind("<FocusOut>", self._paramChanged)

    def _paletteChanged(self, *args):
        self._drawPalette()
        info = "%s has %s colors." % (self.getColorPalette(), cm.get_cmap(self.getColorPalette()).N)
        self.info.set(info)

    def _drawPalette(self):
        """ Draws the palette using current values: palette, highest, lowest and intervals"""
        # clean first
        for widget in self.palette.winfo_children():
            widget.destroy()

        # get the colors
        colorList = getHexColorList(self.getIntervals(), colorName=self.getColorPalette())

        # Get the label values
        labelValues = splitRange(self.getLowest(), self.getHighest(), self.getIntervals())

        # Draw it
        count = 0
        for color, value in zip(colorList, labelValues):
            self._drawLabel(color, value, count)
            count = count + 1

    ### Events handling
    def _paramChanged(self, event):
        # Handles change in palette Option
        self._drawPalette()

    def _keyPressedOnParams(self, event):
        # handles keypress on params frame
        print("PASA key press")
        self.colorPalette.set("jet")

    def getAvailablePalettes(self):
        """ Returns a list of all available palettes"""
        return list(cm.cmap_d.keys())

    def _selectAllText(self, event):
        """ Select all the text of the widget that triggered the event"""
        # select text
        event.widget.select_range(0, 'end')
        # move cursor to the end
        event.widget.icursor('end')
        # stop propagation
        # return 'break'


class FormulaDialog(dialog.Dialog):
    """ This will assist users to create a formula based on class attibutes.
    """

    def __init__(self, parentWindow, set, formula=""):
        """
            :param set: set with the items to use in the formulae
            :param formula: initial formulate to load
        """
        self.set = set
        self.item = set.getFirstItem()
        self.formula = tk.StringVar(value=formula)
        self.formula.trace("w", self.evaluateFormula)
        self.formulaTxt = None # To be instantiated later.
        self.info = tk.StringVar()
        dialog.Dialog.__init__(self, parentWindow, "Formula wizard", default="None")


    # Getters
    def getFormula(self):
        return self.formula.get()

    def getItem(self):
        return self.item

    ### ----- GUI methods ------ ###
    def _createTree(self, parent):
        provider = AttributesTreeProvider(self.item)
        self.tree = BoundTree(parent, provider)
        self.tree.itemDoubleClick = self.addAttributeToFormula
        self.tree.grid(row=1, column=0)

    def addAttributeToFormula(self, event):
        attr = self.tree.getSelectedObjects()[0]
        self._insertText("item." + attr.attrName , self.formulaTxt)

    def body(self, master):
        """ Draws the main frame of the dialog"""
        body = tk.Frame(self)
        body.grid(row=0, column=0, sticky="news")
        body.grid_columnconfigure(1, weight=10)

        # GUI attributes
        self.attributes = tk.Frame(body)
        self.attributes.grid(row=0, column=0, sticky='nes',
                          padx=5, pady=5)

        # Formula
        self.formulaFrame = tk.Frame(body)
        self.formulaFrame.grid(row=0, column=1, sticky='news',
                         padx=5, pady=5)
        #self.params.bind("<Key>", self._keyPressedOnParams)

        self._fillParams()


    def _fillParams(self):

        # Add the formula text
        formulaLbl = tk.Label(self.formulaFrame, text="Formula:")
        formulaLbl.grid(row=0, column=0, sticky="w")
        self.formulaTxt = self._addEntry(1, 0, self.formula,  width=100)

        self._createTree(self.attributes)

        # Label to show formula result or error
        infoLabel = tk.Label(self.formulaFrame, textvariable=self.info)
        infoLabel.grid(row=2, column=0, columnspan=2, sticky="news")

    def _addEntry(self, row, column, var, width=None):

        newEntry = tk.Entry(self.formulaFrame, textvariable=var, width=width)
        newEntry.grid(row=row, column=column, sticky="news")
        newEntry.lift()
        # Bind events
        newEntry.bind("<Return>", self._paramChanged)
        newEntry.bind("<FocusOut>", self._paramChanged)

        return newEntry

    def evaluateFormula(self, name='', index='', mode=''):
        try:
            item = self.item
            result = eval(self.formula.get())
            self.info.set(result)
        except Exception as e:
            self.info.set(str(e))

    ### Events handling
    def _paramChanged(self, event):
        # Handles change in palette Option
        self.evaluateFormula()

    def _selectAllText(self, event):
        """ Select all the text of the widget that triggered the event"""
        # select text
        event.widget.select_range(0, 'end')
        # move cursor to the end
        event.widget.icursor('end')
        # stop propagation
        # return 'break'

    def _insertText(self, text, widget):
        """ Insert the passed text in the entry at current cursor location"""
        # Get inserting position
        pos = widget.index(tk.INSERT)

        # Insert he text
        widget.insert(pos, text)

def insertText (target, textToInsert, position):
    """ Inserts a text into another at a position

    :param target: text to do the insertion on
    :param textToInsert: text to be inserted
    :param position: position where to insert the new text"""

    return target[:position] + textToInsert + target[position:]