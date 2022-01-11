# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import os

import pyworkflow.viewer as pwviewer
import pyworkflow.protocol.params as params
import pyworkflow.gui.plotter as plotter

import pwem.constants as emcts
import pwem.protocols as emprot
from pwem import emlib, splitRange
from pwem.viewers.viewer_chimera import mapVolsWithColorkey, chimeraScriptFileName, Chimera, ChimeraView
from pwem.wizards import ColorScaleWizardBase

from .plotter import EmPlotter, plt


class LocalResolutionViewer(pwviewer.ProtocolViewer):
    """
    Visualization tools for local resolution results.

    """
    binaryCondition = ('(colorMap == %d) ' % emcts.COLOR_OTHER)

    def __init__(self, *args, **kwargs):
        pwviewer.ProtocolViewer.__init__(self, **kwargs)

    def getImgData(self, imgFile, minMaskValue=0.1, maxMaskValue=None):
        import numpy as np
        # if image ends in .mrc or .map :mrc
        if imgFile[-4:] == ".mrc" or imgFile[-4:] == ".map":
            imgFile += ":mrc"
        img = emlib.image.ImageHandler().read(imgFile)
        imgData = img.getData()
        voldim = (img.getDimensions())[:-1]

        if maxMaskValue:
            imgData = np.ma.masked_where(imgData > maxMaskValue, imgData, copy=False)
        maxRes = np.amax(imgData)

        if minMaskValue:
            imgData = np.ma.masked_where(imgData < minMaskValue, imgData, copy=False)
        minRes = np.amin(imgData)

        return imgData, minRes, maxRes, voldim

    def getSlice(self, index, volumeData):
        return int(index * volumeData.shape[0] / 9)

    def getSliceImage(self, volumeData, sliceNumber, dataAxis):
        if dataAxis == 'y':
            imgSlice = volumeData[:, sliceNumber, :]
        elif dataAxis == 'x':
            imgSlice = volumeData[:, :, sliceNumber]
        else:
            imgSlice = volumeData[sliceNumber, :, :]

        return imgSlice

    def createChimeraScript(self, scriptFile, fnResVol,
                            fnOrigMap, sampRate, numColors=13, lowResLimit=None, highResLimit=None, showAxis=True):

        imageFile = os.path.abspath(fnResVol)

        _, minRes, maxRes, voldim = self.getImgData(imageFile)

        # Narrow the color range to the highest resolution range
        if lowResLimit is None:
            lowResLimit = min(maxRes, minRes + 5)
        if highResLimit is None:
            highResLimit = minRes

        stepColors = splitRange(highResLimit, lowResLimit,
                                splitNum=numColors)
        colorList = plotter.getHexColorList(len(stepColors), self._getColorName())

        fnVol = os.path.abspath(fnOrigMap)

        mapVolsWithColorkey(fnVol,
                            imageFile,
                            stepColors,
                            colorList,
                            voldim,
                            volOrigin=None,
                            step=-1,
                            sampling=sampRate,
                            scriptFileName=scriptFile,
                            bgColorImage='white',
                            showAxis=showAxis)

    def _getColorName(self):
        return self.colorMap.get()


class RMSDPerResidueViewer(pwviewer.ProtocolViewer):
    """ Viewer for plot the histogram of per-residue RMSD of a SetOfAtomStruct. """
    _label = 'RMSD per-residue viewer'
    _targets = [emprot.ProtRMSDAtomStructs]
    _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]

    def __init__(self, **kwargs):
      pwviewer.ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
      form.addSection(label='Visualization of per-residue RMSD')
      ColorScaleWizardBase.defineColorScaleParams(form, defaultLowest=0, defaultHighest=2, defaultIntervals=20,
                                                  defaultColorMap='seismic_r')
      group = form.addGroup('Visualization in Chimera')
      group.addParam('displayStruct', params.LabelParam,
                     label='Display output AtomStruct RMSD: ',
                     help='Display one of the AtomStruct in the set, coloured by the RMSD per residue.\n'
                          'The color palette, intervals, lowest and highest values can be chosen in the '
                          'parameters above.'
                     )

      group = form.addGroup('RMSD histogram')
      group.addParam('displayHistogram', params.LabelParam,
                     label='Display histogram of RMSD per-residue: ',
                     help='Display a histogram with the number of RMSD-per residue.\n'
                          'The intervals, lowest and highest values can be chosen in the parameters above.\n'
                          'The bars are coloured by their RMSD value: lower 25% (blue), 25-75% (grey), higher 25% (red)'
                     )

    def _getVisualizeDict(self):
      return {
        'displayStruct': self._showChimera,
        'displayHistogram': self._showHistogram,
      }

    def _showChimera(self, paramName=None):
      obj = self.protocol
      chimScript = obj._getExtraPath(chimeraScriptFileName)
      f = open(chimScript, "w")
      f.write('cd %s\n' % os.getcwd())
      f.write("cofr 0,0,0\n")  # set center of coordinates

      # building coordinate axes
      _inputStruct = obj.outputAtomStructs.getFirstItem()
      _inputVol = _inputStruct.getVolume()
      if _inputVol is not None:
        dim, sampling = _inputVol.getDim()[0], _inputVol.getSamplingRate()

        f.write("open %s\n" % _inputVol.getFileName())
        x, y, z = _inputVol.getOrigin(force=True).getShifts()
        f.write("volume #%d style surface voxelSize %f\nvolume #%d origin "
                "%0.2f,%0.2f,%0.2f\n"
                % (1, sampling, 1, x, y, z))
        strId = 3
      else:
        dim, sampling = 150., 1.
        strId = 2

      tmpFileName = os.path.abspath(obj._getExtraPath("axis_input.bild"))
      Chimera.createCoordinateAxisFile(dim,
                                       bildFileName=tmpFileName,
                                       sampling=sampling)
      f.write("open %s\n" % tmpFileName)

      # Open atomstruct and color it by the bfactor (which is actually the DAQ score)
      f.write("open %s\n" % _inputStruct.getFileName())
      defAttrFile = self.reformatRMSDFile(strId)
      f.write("defattr %s\n" % defAttrFile)

      stepColors, colorList = self.getColors()
      scolorStr = ''
      for step, color in zip(stepColors, colorList):
        scolorStr += '%s,%s:' % (step, color)
      scolorStr = scolorStr[:-1]

      f.write("color byattribute perResidueRMSD palette {} range 0,4".format(scolorStr))
      f.close()
      view = ChimeraView(chimScript)
      return [view]

    def _showHistogram(self, paramName=None):
      obj = self.protocol
      self.plotter = EmPlotter(x=1, y=1, windowTitle='RMSD', figure='active')
      rmsdPerResidueDic = obj.parsePerResidueRMSD()
      rmsdValues = list(rmsdPerResidueDic.values())

      a = self.plotter.createSubPlot("Per-residue RMSD", "RMSD", "Count")
      low, high = self.lowest.get(), self.highest.get()
      a.set_xlim([low, high])

      n = 2
      mult = 10 ** n
      stepSize = int(round(high / self.intervals.get(), n) * mult)
      bins = [i / mult for i in range(int(low * mult), int(high * mult), stepSize)]
      _, _, bars = a.hist(rmsdValues, bins=bins, linewidth=1, label="Map", rwidth=0.9)
      for bar in bars:
        if bar.get_x() < high/4:
          bar.set_facecolor('blue')
        elif bar.get_x() < 3*high/4:
          bar.set_facecolor('grey')
        else:
          bar.set_facecolor('red')
      a.grid(True)
      self.show()

    def reformatRMSDFile(self, strId):
      defattrFile = self.protocol._getExtraPath('overAllRMSD.defattr')
      with open(defattrFile, 'w') as f:
        with open(self.protocol.getPerResidueRMSDFile()) as fIn:
          for i in range(2):
            f.write(fIn.readline())
          for line in fIn:
            spec, value = line.strip().split()
            f.write('\t#{}/{}:{}\t{}\n'.format(strId, spec.split('.')[1],
                                               spec.split('.')[0][1:], value))
      return defattrFile

    def show(self, block=True):
      self.plotter.show(block=block)

    def getColors(self):
      stepColors = splitRange(self.highest.get(), self.lowest.get(),
                              splitNum=self.intervals.get())
      colorList = plotter.getHexColorList(len(stepColors), self.colorMap.get())
      return stepColors, colorList
