# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

import pyworkflow.protocol.params as params
from pwem import UNIT_ANGSTROM_SYMBOL
from pyworkflow import VERSION_3_0

from pwem.protocols import EMProtocol
import pwem.objects as emobj
import pyworkflow.utils as pwutils
from os.path import basename


class ProtOrigSampling(EMProtocol):
    """ Modify the origin and sampling values assigned to a 3D map
    """
    _label = 'assign Orig & Sampling'
    _program = ""
    _lastUpdateVersion = VERSION_3_0

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inVolume', params.PointerParam, label="Input Volume", important=True,
                      pointerClass='Volume',
                      help='A new object volume will be created with the assigned sampling rate and origin. '
                           'No new binary file will be created')
        form.addParam('copyFiles', params.BooleanParam,
                      label="copy (true)/link(false) volume:",
                      help="Option YES:\nA new volume file will be copied "
                      "otherwise a link to the input volume is made\n"
                      "default = false", expertLevel=params.LEVEL_ADVANCED,
                      default=False)
        form.addParam('setSampling', params.BooleanParam,
                      label="Set SamplingRate",
                      help="Option YES:\nA new volume object will be created with "
                           "the given SamplinRate. "
                           "This SamplingRate will NOT be set in the map file header.\n\n",
                      default=False)
        form.addParam('samplingRate', params.FloatParam,
                      condition='setSampling',
                      label=pwutils.Message.LABEL_SAMP_RATE)
        form.addParam('setOrigCoord', params.BooleanParam,
                      label="Set origin of coordinates",
                      help="Option YES:\nA new volume object will be created with "
                           "the given ORIGIN of coordinates. "
                           "This ORIGIN will NOT be set in the map file header.\n\n",
                      default=False)
        line = form.addLine('Offset(%s)' % UNIT_ANGSTROM_SYMBOL,
                            help="A wizard will suggest you possible "
                                 "coordinates for the ORIGIN. In MRC volume "
                                 "files, the ORIGIN coordinates will be "
                                 "obtained from the file header.\n "
                                 "In case you prefer set your own ORIGIN "
                                 "coordinates, write them here. You have to "
                                 "provide the map center coordinates in "
                                 "Angstroms (pixels x sampling).\n",
                            condition='setOrigCoord')
        line.addParam('x', params.FloatParam, condition='setOrigCoord',
                      label="x",
                      help="offset along x axis (Angstroms)")
        line.addParam('y', params.FloatParam, condition='setOrigCoord',
                      label="  y",
                      help="offset along y axis (Angstroms)")
        line.addParam('z', params.FloatParam, condition='setOrigCoord',
                      label="  z",
                      help="offset along z axis (Angstroms)")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('assignStep')

    # --------------------------- STEPS functions -----------------------------

    def assignStep(self):
        # Create a Volume object
        self.inVol = self.inVolume.get()
        self.outVol = emobj.Volume()
        self.outVol.copy(self.inVol)

        # Set sampling Rate (or copy from input)
        if self.setSampling.get():
            samplingRate = self.samplingRate
        else:
            samplingRate = self.inVol.getSamplingRate()

        self.outVol.setSamplingRate(samplingRate)

        # set Origin
        if self.setOrigCoord.get():
            origin = emobj.Transform()
            origin.setShifts(self.x, self.y, self.z)
        else:
            origin = self.inVol.getOrigin()

        self.outVol.setOrigin(origin)

        # copy or link
        copyOrLink = self.getCopyOrLink()
        fileName = self.inVol.getFileName()
        fileName = fileName.replace(':mrc', '')
        imgBase = basename(fileName)
        imgDst = self._getExtraPath(imgBase)
        copyOrLink(os.path.abspath(fileName), imgDst)
        self.outVol.setLocation(imgDst)

        # save
        self._defineOutputs(outputVolume=self.outVol)
        self._defineSourceRelation(self.inVol, self.outVol)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        message = []
        # message.append("You must set a path to export.")
        # Nothing to append
        return message

    def _summary(self):
        # volume copy or linked
        message = []
        if self.setSampling.get():
            message.append("New Sampling: %f\n"%
                                   self.samplingRate)
        if self.setOrigCoord.get():
            message.append("New Origin: %f %f %f\n"%
                               (self.x, self.y, self.z))
        return message

    def _methods(self):
        return []

# --------------------------- UTILS functions ---------------------------------

    def getCopyOrLink(self):
        # Set a function to copyFile or createLink
        # depending in the user selected option
        if self.copyFiles:
            return pwutils.copyFile
        else:
            return pwutils.createAbsLink

    def getFnPath(self, label='volume'):
        return os.path.join(self.filesPath.get(),
                            self._getFileName(label))

