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
import pwem.constants as emcts
from pwem import emlib, splitRange
from pwem.viewers.viewer_chimera import mapVolsWithColorkey
import pyworkflow.gui.plotter as plotter


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
