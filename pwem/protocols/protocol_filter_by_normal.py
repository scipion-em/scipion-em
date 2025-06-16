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
import math
import numpy as np

from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, EnumParam, FloatParam, 
                                        StringParam, Range)
from pwem.protocols import EMProtocol
from pwem.objects import EMSet, Transform, Matrix
from pwem.convert.transformations import euler_from_matrix

def deep_getattr(obj, attr: str):
    for key in attr.split('.'):
        obj = getattr(obj, key)
    return obj



class ProtSetFilterByNormalOutputs(Enum):
    Output = EMSet

class ProtSetFilterByNormal(EMProtocol):
    """
    This protocol allows to filter any set with items that have a transform
    associated such that only top views or side views remain.
    """
    _label = 'filter set by normals'

    CHOICE_SIDE_VIEW = 'Side views'
    CHOICE_TOP_VIEW = 'Top views'
    CHOICES = [
        CHOICE_SIDE_VIEW,
        CHOICE_TOP_VIEW
    ]

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input', PointerParam, pointerClass=EMSet, label='Input',
                      help='Set which items will be filtered.')
        form.addParam('selection', EnumParam, label='Selection',
                      choices=self.CHOICES, default=0,
                      help='Select which views are kept. Side views keep items '
                           'with tilt angles near 90 or 270 degrees. Top '
                           'views keep items with tilt angles near 0 or 180 '
                           'degrees.')
        form.addParam('tolerance', FloatParam, label='Tolerance (deg)',
                      default=15.0, validators=[Range(0.0, 90.0)],
                      help='Select the tolerance used for keeping elements in '
                           'degrees. This defines the deviation from the angle '
                           'mentioned in the selection. Must be in [0, 90] '
                           'range.')
        form.addParam('transformField', StringParam, label='Transform field',
                      default='_transform', expertLevel=LEVEL_ADVANCED,
                      help='Field used to retrieve the transform matrix of '
                      'the items.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self):
        inputSet = self._getInputSet()
        outputSet = inputSet.createCopy(self._getPath(), copyInfo=True)
        
        ANGLE_TARGETS = [90.0, 0.0]
        target = math.radians(ANGLE_TARGETS[self.selection.get()])
        tolerance = math.sin(math.radians(self.tolerance.get()))
        for item in inputSet:
            matrix = self._getMatrix(item)
            _, tilt, _ = euler_from_matrix(matrix, 'szyz')
            tilt -= target
            s = abs(math.sin(tilt))
            if s < tolerance:
                outputSet.append(item)
            
        self._defineOutputs(**{ProtSetFilterByNormalOutputs.Output.name: outputSet})
        self._defineSourceRelation(self.input, outputSet)
            
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        inputSet = self._getInputSet()
        item = inputSet.getFirstItem()
        try:
            self._getMatrix(item)
            
        except (AttributeError, TypeError, ValueError):
            errors.append('Input items must have a valid transform associated. '
                          'Try modifying "Matrix Field" parameter if necessary')
                    
        return errors
    
    # --------------------------- UTILS functions -----------------------------
    def _getInputSet(self) -> EMSet:
        return self.input.get()
    
    def _getMatrix(self, item) -> np.ndarray:
        transform = deep_getattr(item, self.transformField.get())
        if isinstance(transform, Transform):
            matrix = transform.getMatrix()
        elif isinstance(transform, Matrix):
            matrix = transform.getMatrix()
        else:
            raise TypeError('Unknown transform type')
        
        if matrix.shape == (4, 4):
            matrix = matrix[:3,:3]
        
        if matrix.shape != (3, 3):
            raise ValueError('Invalid matrix shape')
        
        return matrix
