# **************************************************************************
# *
# * Authors:     Pablo Conesa(pconesa@cnb.csic.es)
# *              Roberto Marabini(roberto@cnb.csic.es)
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
from pwem.protocols import EMProtocol
from pwem.objects.data import SetOfParticles

from pwem.convert.transformations import (rotation_matrix,
      angle_between_vectors, vector_product)
import numpy as np

from pwem.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                                     SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r,
                                     SCIPION_SYM_NAME)

class ProtSetEditor(EMProtocol):
    """
    Protocol to edit attributes of all the items of a set using a formula.
    This could be useful for editing some values in the set. Use this
    protocol with extreme care, you can easily produce a set that is
    not consistent. Some predefined operation are offered as given
    a set of particles with projection alignment modify the alignment matrix
    so when reconstructed the volume is rotated.
    """
    _label = 'edit set'
    CHOICE_FORMULA = 0
    CHOICE_ROTATE_VECTOR = 1
    CHOICE_ROTATE_ICOSAHEDRAL = 2
    CHOICE_LABEL = {CHOICE_FORMULA: 'formula',
                    CHOICE_ROTATE_VECTOR: 'rotate to vector',
                    CHOICE_ROTATE_ICOSAHEDRAL: 'rotate icosahedral'}
    LOCAL_SYM_NAME = {SYM_I222: 'I1',
                      SYM_I222r: 'I2',
                      SYM_In25: 'I3',
                      SYM_In25r: 'I4',
                      SYM_I2n3: 'I5',
                      SYM_I2n3r: 'I6',
                      SYM_I2n5: 'I7',
                      SYM_I2n5r: 'I8'}

    def _defineParams(self, form):
        """
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """
        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam,
                      pointerClass='EMSet',
                      label='Set to edit',
                      help='Set which items will be modified.')
        form.addParam('operation', params.EnumParam,
                      choices=[self.CHOICE_LABEL[self.CHOICE_FORMULA],
                               self.CHOICE_LABEL[self.CHOICE_ROTATE_VECTOR],
                               self.CHOICE_LABEL[self.CHOICE_ROTATE_ICOSAHEDRAL]
                               ],
                      default = self.CHOICE_FORMULA,
                      label="Select operation",
                      help="Select operation to be performed in the set.\n"
                           " *rotate to vector* modifies the alignment matrix  "
                           " so a reconstruction made from the images produces a "
                           " rotated reconstruction\n"
                           " *rotated icosahedral symmetry* is identical to the previous case"
                           " but instead the provide the initial and final vectors, "
                           " the target and source symmetries are provided\n"
                           " *formula* will apply the formula to the chosen attribute"
                           " (i.e, multiply resolution by 2 )"
                      )

        form.addParam('formula', params.StringParam, label="Formula",
                      condition = "operation==%d" % self.CHOICE_FORMULA,
                      help='A python code compatible with eval, where item represents each of '
                           'the elements of the set. E.g.: item._resolution.set(item._resolution.get() +1).'
                           'You could also use modules like "import numpy;  item._resolution .... "')

        line = form.addLine('source vector',
                            condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                            help='Vector will be rotated until overlaps target vector')
        line.addParam('xs', params.FloatParam, default=0,
                      condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="x",
                      help="X Coordinate ")
        line.addParam('ys', params.FloatParam, default=0,
                      condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="y",
                      help="Y dim ")
        line.addParam('zs', params.FloatParam, default=1,
                      condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="z",
                      help="Z dim ")

        line = form.addLine('target vector',
                            condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                            help='source vector will be rotated until '
                                 'overlaps target vector')
        line.addParam('xt', params.FloatParam, default=1,
                      condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="x",
                      help="X Coordinate ")
        line.addParam('yt', params.FloatParam, default=0,
                      condition=
                      "operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="y",
                      help="Y dim ")
        line.addParam('zt', params.FloatParam, default=0,
                      condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="z",
                      help="Z dim ")

        line = form.addLine('Icosahedral symmetry',
                            condition=
                            "operation==%d" % self.CHOICE_ROTATE_ICOSAHEDRAL,
                            help='convert between icosahedral symmetries')
        line.addParam("originSymmetryGroup", params.EnumParam,
                      condition="operation==%d" % self.CHOICE_ROTATE_ICOSAHEDRAL,
                      label="origin",
                      help = "Source symmetry. "
                             "Only implemented for icosahedral symmetry",
                      choices=[self.LOCAL_SYM_NAME[SYM_I222] +
                               " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                               self.LOCAL_SYM_NAME[SYM_I222r] +
                               " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",
                               self.LOCAL_SYM_NAME[SYM_In25] +
                               " (" + SCIPION_SYM_NAME[SYM_In25] + ")",
                               self.LOCAL_SYM_NAME[SYM_In25r] +
                               " (" + SCIPION_SYM_NAME[SYM_In25r] + ")",
                               self.LOCAL_SYM_NAME[SYM_I2n3] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n3] + ")",
                               self.LOCAL_SYM_NAME[SYM_I2n3r] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n3r] + ")",
                               self.LOCAL_SYM_NAME[SYM_I2n5] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n5] + ")",
                               self.LOCAL_SYM_NAME[SYM_I2n5r] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n5r] + ")"
                               ],
                      default=SYM_I222r - SYM_I222,
                      )
        line.addParam("targetSymmetryGroup", params.EnumParam,
                      condition="operation==%d" % self.CHOICE_ROTATE_ICOSAHEDRAL,
                      label="target",
                      help = "Target symmetry. "
                             "Only implemented for icosahedral symmetry",
                      choices=[self.LOCAL_SYM_NAME[SYM_I222] +
                               " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                               self.LOCAL_SYM_NAME[SYM_I222r] +
                               " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",
                               self.LOCAL_SYM_NAME[SYM_In25] +
                               " (" + SCIPION_SYM_NAME[SYM_In25] + ")",
                               self.LOCAL_SYM_NAME[SYM_In25r] +
                               " (" + SCIPION_SYM_NAME[SYM_In25r] + ")",
                               self.LOCAL_SYM_NAME[SYM_I2n3] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n3] + ")",
                               self.LOCAL_SYM_NAME[SYM_I2n3r] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n3r] + ")",
                               self.LOCAL_SYM_NAME[SYM_I2n5] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n5] + ")",
                               self.LOCAL_SYM_NAME[SYM_I2n5r] +
                               " (" + SCIPION_SYM_NAME[SYM_I2n5r] + ")"
                               ],
                      default=SYM_I222r - SYM_I222,
                      )

    def _insertAllSteps(self):
        operation =  self.operation.get()
        if (operation == self.CHOICE_FORMULA):
            self._insertFunctionStep('formulaStep')
        elif (operation == self.CHOICE_ROTATE_VECTOR):
            self._insertFunctionStep('rotateStep')
        elif (operation == self.CHOICE_ROTATE_ICOSAHEDRAL):
            self._insertFunctionStep('rotateIcoStep')


    def rotateIcoStep(self):
        """
        Compute rotation matrix from one icosahedral
        symmetry to another and apply it
        to projection directions. Reconstructed volume
        should change symmetry
        """
        from pwem.convert.symmetry import Icosahedron
        origSym = self.originSymmetryGroup.get() + SYM_I222
        targetSym = self.targetSymmetryGroup.get() + SYM_I222
        ico = Icosahedron(orientation = SCIPION_SYM_NAME[origSym][1:])
        matrix = ico.coordinateSystemTransform(SCIPION_SYM_NAME[origSym][1:],
                                               SCIPION_SYM_NAME[targetSym][1:])
        # convert to numpy array
        # and add extra row
        matrix =  np.array(matrix)
        matrix = np.append(matrix, [[0, 0, 0, 1]], axis=0)

        inputSet = self.inputSet.get()
        modifiedSet = inputSet.createCopy(self._getExtraPath(), copyInfo=True)
        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            transformation = item.getTransform()
            transformation.composeTransform(matrix)
            modifiedSet.append(item)
        self.createOutput(self.inputSet, modifiedSet)

    def createOutput(self, inputSet, modifiedSet):
        outputArgs = {inputSet.getExtended(): modifiedSet}
        self._defineOutputs(**outputArgs)

    def rotateStep(self):
        """
        Compute rotation matrix between user provided vectors
         (x,y,z) and (x',y',z') and apply it
        to projection directions. Reconstructed volume
        should rotate
        """
        v_source = np.array([self.xs.get(), self.ys.get(), self.zs.get()])
        v_source = v_source / np.linalg.norm(v_source)
        v_target = np.array([self.xt.get(), self.yt.get(), self.zt.get()])
        v_target = v_target / np.linalg.norm(v_target   )
        matrix = rotation_matrix(angle_between_vectors(v_source, v_target),
                                 vector_product(v_source, v_target))
        inputSet = self.inputSet.get()
        modifiedSet = inputSet.createCopy(self._getExtraPath(), copyInfo=True)
        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            transformation = item.getTransform()
            transformation.composeTransform(matrix)
            modifiedSet.append(item)

        self.createOutput(self.inputSet, modifiedSet)


    def formulaStep(self):
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

        self.createOutput(self.inputSet, modifiedSet)


    def _validate(self):
        errors = []
        inputSet = self.inputSet.get()
        operation = self.operation.get()
        if operation != self.CHOICE_FORMULA:
            if not isinstance(inputSet, SetOfParticles):
                errors.append("The input data set is not a set of projections")
            elif not inputSet.hasAlignmentProj():
                errors.append("The input data set does not have alignment 3D")
        return errors
