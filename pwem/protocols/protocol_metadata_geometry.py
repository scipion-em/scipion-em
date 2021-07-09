# **************************************************************************
# *
# * Authors:     Roberto Marabini(roberto@cnb.csic.es)
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

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol
from pwem.constants import (SYM_CYCLIC, SYM_DIHEDRAL_X,
                            SYM_DIHEDRAL_Y, SYM_TETRAHEDRAL,
                            SYM_TETRAHEDRAL_Z3, SYM_OCTAHEDRAL,
                            SYM_OC  SYM_I222, SYM_I222r, SYM_In25,
                            SYM_In25r, SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r,
                            SCIPION_SYM_NAME)


class ProtMetadataGeometry(EMProtocol):
    """
    Protocol Metadata Geometry. Select/delete particles
    that satisfy some constrains related with geometry
    """
    _label = 'geometry metadata operator'

    def _defineParams(self, form):
        """
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """

        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam, pointerClass='EMSet',
                      label='Set to edit',
                      help='Set which items will be modified.')
        desplegable
        vector
        form.addParam('originSymmetryGroup', EnumParam,
                      choices=[LOCAL_SYM_NAME[SYM_I222] +
                               " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                               SCIPION_SYM_NAME[SYM_I222r] +
                               " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",
                               SCIPION_SYM_NAME[SYM_In25] +
                               " (" + SCIPION_SYM_NAME[SYM_In25] + ")",
                               SCIPION_SYM_NAME[SYM_In25r] +
                               " (" + SCIPION_SYM_NAME[SYM_In25r] + ")",
                               SCIPION_SYM_NAME[SYM_I2n3] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n3] + ")",
                               SCIPION_SYM_NAME[SYM_I2n3r] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n3r] + ")",
                               SCIPION_SYM_NAME[SYM_I2n5] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n5] + ")",
                               SCIPION_SYM_NAME[SYM_I2n5r] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n5r] + ")"
                               ],
                      default=SYM_I222r - SYM_I222,
                      label="Symmetry",
                      help="Select the current symmetry of your atomic structure./n"
                           "See https://github.com/I2PC/xmipp-portal/wiki/Symmetry"
                           "Symmetry for a description of the symmetry groups "
                      )

    def _insertAllSteps(self):
        self._insertFunctionStep('editItemsStep')

    def editItemsStep(self):
        """
        Goes through all items in the input set and applies the formula to each of them using exec.
        Complex python code could be run separating lines with ;  To use numpy you could do
        import numpy; item._resolution.set(numpy.random.randint(10))
        """
        inputSet = self.inputSet.get()

        modifiedSet = inputSet.createCopy(self._getExtraPath(), copyInfo=True)

        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            exec(self.formula.get())
            modifiedSet.append(item)

        outputArgs = {self.inputSet.getExtended(): modifiedSet}
        self._defineOutputs(**outputArgs)
