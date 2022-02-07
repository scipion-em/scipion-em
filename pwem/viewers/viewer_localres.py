# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
import numpy as np

import pyworkflow.viewer as pwviewer
import pyworkflow.protocol.params as params
import pyworkflow.gui.plotter as plotter
import pyworkflow.gui.dialog as dialog

import pwem.constants as emcts
import pwem.protocols as emprot
from pwem import emlib, splitRange
from pwem.objects import AtomStruct, SetOfAtomStructs
from pwem.viewers.viewer_chimera import mapVolsWithColorkey, Chimera, ChimeraView, generateColorLegend
from pwem.convert.atom_struct import *

from .plotter import EmPlotter, plt

CHIMERA_ERROR = 'Chimera program is not found were it was expected: \n\n{}\n\n' \
                'Either install ChimeraX in this path or install our ' \
                'scipion-em-chimera plugin'.format(Chimera.getProgram())


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


class ChimeraAttributeViewer(pwviewer.ProtocolViewer):
    """ Viewer for attributes of a SetOfAtomStruct or AtomStruct.
    Includes visualization in chimera and in histograms"""
    _label = 'Atomic structure attributes viewer'
    _targets = []
    _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]

    def __init__(self, **kwargs):
        pwviewer.ProtocolViewer.__init__(self, **kwargs)
        if not chimeraInstalled():
            print(CHIMERA_ERROR)

    def _defineParams(self, form):
      form.addSection(label='Visualization of attributes')
      form.addParam('attrName', params.EnumParam,
                    choices=self._getStructureAttributes(), default=self._getAttributeIndex(self.protocol._ATTRNAME),
                    label='Attribute to display: ',
                    help='Display this attribute of the structure'
                    )
      group = form.addGroup('Attribute histogram')
      group.addParam('displayHistogram', params.LabelParam,
                     label='Display histogram of attribute: ',
                     help='Display a histogram with the count of the attribute.\n'
                          'The color palette, intervals, lowest and highest values can be chosen in the '
                          'parameters below. If the highest and lowest values you input are the same, '
                          'the lowest and higher values in the date are used.'
                     )

      group = form.addGroup('Visualization in Chimera')
      group.addParam('displayStruct', params.LabelParam,
                     label='Display structure and color by attribute in Chimera: ',
                     help='Display the AtomStruct, coloured by the attribute.\n'
                          'The color palette, intervals, lowest and highest values can be chosen in the '
                          'parameters below. If the highest and lowest values you input are the same, '
                          'the lowest and highest values in the date are used.'
                     )
      group.addParam('residuesAverages', params.BooleanParam,
                     default=True, label='Average attribute on residues: ',
                     help='In case the recipients are the atoms, their values are averaged on the residues.\n'
                     )

      from pwem.wizards.wizard import ColorScaleWizardBase
      group = form.addGroup('Color settings')
      ColorScaleWizardBase.defineColorScaleParams(group, defaultLowest=0, defaultHighest=2, defaultIntervals=21,
                                                  defaultColorMap='RdBu')

    def _getOutputObject(self):
        '''Return the output object of the protocol (SetOfAtomStruct or AtomStruct)
        To change in other viewers with different output name'''
        return getattr(self.protocol, self.protocol._OUTNAME)

    def getAtomStructObject(self):
        obj = self._getOutputObject()
        if issubclass(type(obj), SetOfAtomStructs):
          _inputStruct = obj.getRepresentative()
          if _inputStruct is None:
              _inputStruct = obj.getFirstItem()
        elif issubclass(type(obj), AtomStruct):
          _inputStruct = obj
        return _inputStruct

    def _getStructureAttributes(self):
        '''Returns a list with the names of the attributes of the output object'''
        obj = self.getAtomStructObject()
        ASH = AtomicStructHandler()
        cifDic = ASH.readLowLevel(obj.getFileName())
        return list(set(cifDic[NAME]))

    def _getStructureRecipient(self):
        '''Returns a list with the names of the attributes of the output object'''
        obj = self.getAtomStructObject()
        ASH = AtomicStructHandler()
        cifDic = ASH.readLowLevel(obj.getFileName())

        recipDic = {}
        for name, recip in zip(cifDic[NAME], cifDic[RECIP]):
            recipDic[name] = recip
        return recipDic[self.getEnumText('attrName')]

    def _getAttributeIndex(self, attrName):
        '''Returns a list with the names of the attributes of the output object'''
        names = self._getStructureAttributes()
        return names.index(attrName)


    def _getVisualizeDict(self):
      return {
        'displayStruct': self._showChimera,
        'displayHistogram': self._showHistogram,
      }

    def _showChimera(self, paramName=None):
      if chimeraInstalled():
          prot = self.protocol
          chimScript = prot._getExtraPath('chimera_script.py')
          f = open(chimScript, "w")
          f.write("from chimerax.core.commands import run\n")

          f.write("run(session, 'cd %s')\n" % os.getcwd())
          f.write("run(session, 'cofr 0,0,0')\n")  # set center of coordinates

          _inputStruct = self.getAtomStructObject()
          _inputVol = _inputStruct.getVolume()
          if _inputVol is not None:
            strId = 3
            dim, sampling = _inputVol.getDim()[0], _inputVol.getSamplingRate()

            f.write("run(session, 'open %s')\n" % _inputVol.getFileName())
            x, y, z = _inputVol.getOrigin(force=True).getShifts()
            f.write("run(session, 'volume #%d style surface voxelSize %f')\n" % (1, sampling))
            f.write("run(session, 'volume #%d origin %0.2f,%0.2f,%0.2f')\n" % (1, x, y, z))
            f.write("run(session, 'hide #!%d models')\n" % 1)

          else:
            dim, sampling = 150., 1.
            strId = 2

          tmpFileName = os.path.abspath(prot._getExtraPath("axis_input.bild"))
          Chimera.createCoordinateAxisFile(dim,
                                           bildFileName=tmpFileName,
                                           sampling=sampling)
          f.write("run(session, 'open %s')\n" % tmpFileName)

          # Open atomstruct and color it by the bfactor (which is actually the DAQ score)
          cifDic = AtomicStructHandler().readLowLevel(_inputStruct.getFileName())
          if self._getStructureRecipient() in ['atoms']:
              f.write("run(session, 'open %s')\n" % self.replaceOcuppancyWithAttribute(cifDic))
              attrColorName = 'occupancy'
              
          elif self._getStructureRecipient() in ['residues'] and \
                  getCifKeyName(cifDic, 'asym') and getCifKeyName(cifDic, 'seq'):
              f.write("run(session, 'open %s')\n" % self.replaceOcuppancyWithAttribute(cifDic))
              attrColorName = 'occupancy'
          else:
              f.write("run(session, 'open %s')\n" % _inputStruct.getFileName())
              defAttrFile = self.createAttributesFile(_inputStruct, strId)
              f.write("run(session, 'open %s')\n" % defAttrFile)
              attrColorName = self.getEnumText('attrName')

          stepColors, colorList = self.getColors()
          scolorStr = ''
          for step, color in zip(stepColors, colorList):
            scolorStr += '%s,%s:' % (step, color)
          scolorStr = scolorStr[:-1]

          average = ''
          if self._getStructureRecipient() == 'atoms' and self.residuesAverages:
              average = 'average residues'

          f.write("run(session, 'color byattribute {} palette {} {}')\n".
                  format(attrColorName, scolorStr, average))
          f.write(generateColorLegend(stepColors, colorList))
          f.write("run(session, 'view')\n")

          f.close()
          view = ChimeraView(chimScript)
          return [view]

      else:
          return [self.warnMessage(CHIMERA_ERROR, 'Chimera not found')]

    def _showHistogram(self, paramName=None):
      attrname = self.getEnumText('attrName')

      cifDic = AtomicStructHandler().readLowLevel(self.getAtomStructObject().getFileName())
      #TODO: admit non float attributes
      names, values = np.array(cifDic[NAME]), np.array(cifDic[VALUE])
      recipient = self._getStructureRecipient()
      attrValues = values[names == attrname]
      attrValues = list(map(float, attrValues))

      self.plotter = EmPlotter(x=1, y=1, windowTitle=attrname)
      a = self.plotter.createSubPlot("{} counts".format(self.getEnumText('attrName')),
                                     self.getEnumText('attrName'), "{} count".format(recipient))
      low, high = self.getValuesRange()
      a.set_xlim([low, high])

      n = 2
      mult = 10 ** n
      stepSize = int(round((high-low) / self.intervals.get(), n) * mult)
      bins = [i / mult for i in range(int(low * mult), int(high * mult), stepSize)]
      _, _, bars = a.hist(attrValues, bins=bins, linewidth=1, label="Map", rwidth=0.9)

      colorSteps, colorList = self.getColors()
      colorList.reverse()
      for bar, colorS, color in zip(bars, colorSteps, colorList):
          bar.set_facecolor(color)

      a.grid(True)
      return [self.plotter]

    ###################### UTILS #########################
    def getValuesRange(self):
        if self.lowest.get() != self.highest.get():
            return self.lowest.get(), self.highest.get()
        else:
            attrname = self.getEnumText('attrName')
            cifDic = AtomicStructHandler().readLowLevel(self.getAtomStructObject().getFileName())
            names, values = np.array(cifDic[NAME]), np.array(cifDic[VALUE])
            recipient = cifDic[RECIP][0]
            attrValues = values[names == attrname]
            attrValues = list(map(float, attrValues))
            return min(attrValues), max(attrValues)

    def getColors(self):
      low, high = self.getValuesRange()
      stepColors = splitRange(high, low, splitNum=self.intervals.get())
      colorList = plotter.getHexColorList(len(stepColors), self.colorMap.get())
      return stepColors, colorList

    def createAttributesFile(self, AS, chiEleId):
        ASH = AtomicStructHandler()
        cifDic = ASH.readLowLevel(AS.getFileName())

        attrName = self.getEnumText('attrName')
        defAttrStr = 'attribute: {}\n'.format(self.getEnumText('attrName'))
        first = True
        for name, recip, spec, value in zip(cifDic[NAME], cifDic[RECIP],
                                            cifDic[SPEC], cifDic[VALUE]):
            if name == attrName:
                if first:
                    defAttrStr += 'recipient: {}\n'.format(recip)
                    first = False
                defAttrStr += '\t#{}/{}\t{}\n'.format(chiEleId, spec, value)

        defattrFile = self.protocol._getExtraPath('{}.defattr'.format(attrName))
        with open(defattrFile, 'w') as f:
            f.write(defAttrStr)
        return defattrFile

    def replaceOcuppancyWithAttribute(self, cifDic):
        '''Instead of defining the atribute in a defattr file, switch it with the occupancy and color by it.
        It notably speeds up chimera colouring'''
        attrName = self.getEnumText('attrName')
        outFile = self.protocol._getExtraPath('chimeraAttribute_{}.cif'.format(attrName))
        if not os.path.exists(outFile):
            ASH = AtomicStructHandler()

            names, values, specs = np.array(cifDic[NAME]), np.array(cifDic[VALUE]), np.array(cifDic[SPEC])
            recipient = self._getStructureRecipient()
            attrValues, attrSpecs = values[names == attrName], specs[names == attrName]
            if recipient == 'atoms':
                cifDic['_atom_site.occupancy'] = attrValues
            elif recipient == 'residues':
                atomValues = []
                resDic = self.makeResidueValuesDic(cifDic)

                chainsStr, resStr = getCifKeyName(cifDic, 'asym'), getCifKeyName(cifDic, 'seq')
                for resChain, resNumber in zip(cifDic[chainsStr], cifDic[resStr]):
                    resKey = '{}:{}'.format(resChain, resNumber)
                    if resKey in resDic:
                        #HOH atoms may be ignored and not appear in the cif
                        lastVal = resDic[resKey]
                        atomValues.append(lastVal)
                    else:
                        atomValues.append(lastVal)

                cifDic['_atom_site.occupancy'] = atomValues

            ASH._writeLowLevel(outFile, cifDic)
        return outFile

    def makeResidueValuesDic(self, cifDic):
        resDic = {}
        attrName = self.getEnumText('attrName')
        names, values, specs = np.array(cifDic[NAME]), np.array(cifDic[VALUE]), np.array(cifDic[SPEC])
        attrValues, attrSpecs = values[names == attrName], specs[names == attrName]
        for spec, val in zip(attrSpecs, attrValues):
            resDic[spec] = val
        return resDic

def chimeraInstalled():
  return Chimera.getHome() and os.path.exists(Chimera.getProgram())

def getCifKeyName(cifDic, keyBase):
    base = '_atom_site.{}' + '_{}_id'.format(keyBase)
    options = [base.format('pdbx_auth'), base.format('auth'), base.format('label')]
    for name in options:
      if name in cifDic:
        return name