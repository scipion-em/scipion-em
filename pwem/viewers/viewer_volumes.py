# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
This module implements visualization program
for input volumes.
"""

import os
from shutil import which
from tkinter.messagebox import showerror

import pyworkflow.protocol.params as params
import pyworkflow.viewer as pwviewer

import pwem.objects as emobj
import pwem.protocols as emprot
from pwem.viewers import Chimera, ChimeraView, EmProtocolViewer
from pwem.emlib import image

VOLUME_SLICES = 1
VOLUME_CHIMERA = 0


class viewerProtImportVolumes(EmProtocolViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj. """

    _label = 'viewer input volume'
    _targets = [emprot.ProtImportVolumes, emprot.ProtOrigSampling]
    _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]

    def _defineParams(self, form):
        form.addSection(label='Visualization of input volumes')
        form.addParam('displayVol', params.EnumParam,
                      choices=['chimerax', 'slices'],
                      default=VOLUME_CHIMERA,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Display volume with',
                      help='*chimerax*: display volumes as surfaces with '
                           'ChimeraX.\n *slices*: display volumes as 2D slices '
                           'along z axis.\n If the number of volumes is equal to 1, '
                           'a system of coordinates is shown'
                      )

    def _getVisualizeDict(self):
        return {
            'displayVol': self._showVolumes,
        }

    def _validate(self):
        if (self.displayVol == VOLUME_CHIMERA
                and which(Chimera.getProgram()) is None):
            return ["chimera is not available. "
                    "Either install it or choose option 'slices'. "]
        return []

    # =========================================================================
    # ShowVolumes
    # =========================================================================

    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()

        elif self.displayVol == VOLUME_SLICES:
            return self._showVolumesSlices()

    def _createSetOfVolumes(self):
        try:
            setOfVolumes = self.protocol.outputVolumes
            sampling = self.protocol.outputVolumes.getSamplingRate()
        except Exception as ex:
            setOfVolumes = self.protocol._createSetOfVolumes()
            setOfVolumes.append(self.protocol.outputVolume)
            sampling = self.protocol.outputVolume.getSamplingRate()

        return sampling, setOfVolumes

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        tmpFileNameCMD = self.protocol._getExtraPath("chimera.cxc")
        f = open(tmpFileNameCMD, "w")
        sampling, _setOfVolumes = self._createSetOfVolumes()
        count = 1  # chimeraX stars counting in 1 (chimera in 0)

        if len(_setOfVolumes) == 1:
            # if we have a single volume then create axis
            # as bild file. Chimera must read the bild file first
            # otherwise system of coordinates will not
            # be in the center of the window

            dim = self.protocol.outputVolume.getDim()[0]
            tmpFileNameBILD = os.path.abspath(self.protocol._getExtraPath(
                "axis.bild"))
            Chimera.createCoordinateAxisFile(dim,
                                             bildFileName=tmpFileNameBILD,
                                             sampling=sampling)
            f.write("open %s\n" % tmpFileNameBILD)
            f.write("cofr 0,0,0\n")  # set center of coordinates
            count = 2  # skip first model because is not a 3D map

        for vol in _setOfVolumes:
            localVol = os.path.abspath(image.ImageHandler.removeFileType(
                vol.getFileName()))
            if localVol.endswith("stk"):
                errorWindow(None, "Extension .stk is not supported")
            f.write("open %s\n" % localVol)
            f.write("volume #%d style surface level 0.001 voxelSize %f\n" %
                    (count, sampling))
            count += 1

        if len(_setOfVolumes) > 1:
            f.write('tile\n')
        else:
            x, y, z = vol.getShiftsFromOrigin()
            f.write("volume #2 origin %0.2f,%0.2f,%0.2f\n" % (x, y, z))
            if vol.getHalfMaps():
                for halfMap in vol.getHalfMaps().split(','):
                    if not os.path.abspath(halfMap).endswith(".mrc"):
                        f.write("open %s\n" % (os.path.abspath(halfMap).split(".")[0] + ".mrc"))
                    else:
                        f.write("open %s\n" % os.path.abspath(halfMap))
                    f.write("volume #%d style surface level 0.001 voxelSize %f\n" %
                            (count, sampling))
                    f.write("volume #%d origin %0.2f,%0.2f,%0.2f\n" %
                            (count, x, y, z))
                    f.write("tile\n")
                    count += 1
        f.write("view\n")
        f.close()
        return [ChimeraView(tmpFileNameCMD)]

    def _showVolumesSlices(self):
        # Write an sqlite with all volumes selected for visualization.
        sampling, setOfVolumes = self._createSetOfVolumes()
        # setOfVolumes.setSamplingRate(sampling)
        if len(setOfVolumes) == 1:
            if self.protocol.outputVolume.getHalfMaps():
                setOfVolumes = self.protocol._createSetOfVolumes(suffix='tmp')
                vol = emobj.Volume()
                vol.setFileName(self.protocol.outputVolume.getFileName())
                vol.setSamplingRate(sampling)
                setOfVolumes.append(vol)
                for halfMap in self.protocol.outputVolume.getHalfMaps().split(','):
                    volHalf = emobj.Volume()
                    volHalf.setSamplingRate(sampling)
                    if not halfMap.endswith(".mrc"):
                        volHalf.setFileName(halfMap.split(".")[0] + ".mrc")
                    else:
                        volHalf.setFileName(halfMap)
                    setOfVolumes.append(volHalf)
                setOfVolumes.write()
                return [self.objectView(setOfVolumes)]
            else:
                return [self.objectView(self.protocol.outputVolume)]
        return [self.objectView(setOfVolumes)]


def errorWindow(tkParent, msg):
    try:
        showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
    except Exception as ex:
        print("Error:", msg)

