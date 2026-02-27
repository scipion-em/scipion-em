# **************************************************************************
# *
# * Authors:     Oier Lauzirika Zarrabeitia (olauzirika@cnb.csic.es)
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

from enum import Enum
import numpy as np
import xmippLib

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, EnumParam, FloatParam, 
                                        StringParam, Range)
from pwem.protocols import EMProtocol
from pwem.objects import SetOfImages, Image, Transform
from pwem.convert.transformations import euler_from_matrix

def deep_getattr(obj, attr: str):
    for key in attr.split('.'):
        obj = getattr(obj, key)
    return obj



class ProtBreakSymmetryOutputs(Enum):
    Output = SetOfImages

class ProtBreakSymmetry(EMProtocol):
    """
    This protocol allows to break symmetrical angular assignments by shuffling
    the unit cell.
    """
    _label = 'break symmetry'

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input', PointerParam, pointerClass=SetOfImages, label='Input')
        form.addParam('symmetryGroup', StringParam, label='Symmetry group')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self):
        inputImages = self._getInputImages()
        outputImages = inputImages.createCopy(
            self._getPath(), 
            copyInfo=True
        )

        symmetry = xmippLib.SymList()
        symMatrices = np.array(symmetry.getSymmetryMatrices(self.symmetryGroup.get()))
        rng = np.random.default_rng()
        
        symMatrix = np.eye(4)
        for image in inputImages:
            index = rng.integers(0, len(symMatrices))
            symMatrix[:3,:3] = symMatrices[index]
            
            transform: Transform = image.getTransform()
            matrix = transform.getMatrix()
            matrix = symMatrix.T @ matrix
            transform.setMatrix(matrix)
            image.setTransform(transform)
            outputImages.append(image)
          
        self._defineOutputs(**{ProtBreakSymmetryOutputs.Output.name: outputImages})
        self._defineSourceRelation(self.input, outputImages)
            
    # --------------------------- UTILS functions -----------------------------
    def _getInputImages(self) -> SetOfImages:
        return self.input.get()
    