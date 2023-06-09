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

from pwem.convert.transformations import (
    rotation_matrix, angle_between_vectors,
    vector_product)
import numpy as np

from pwem.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                            SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r,
                            SYM_DIHEDRAL_X, SYM_DIHEDRAL_Y,
                            SYM_TETRAHEDRAL_222, SYM_TETRAHEDRAL_Z3,
                            SYM_TETRAHEDRAL_Z3R, SCIPION_SYM_NAME
                            )


class ProtProjectionEditor(EMProtocol):
    """
    This protocol edits the projection directions of all the items of a set of particles using
    a formula. This could be useful for applying geometrical transformation to a set of
    particles.

    Several predefined operation are offered such as
    * apply the rotation matrix define by a origin vector and a target vector
    * rotate the projection vector around a given vector by  A degrees
    * convert between icosahedral symmetries
    * convert between dihedral symmetries
    * convert between tetrahedral symmetries
    """
    _label = 'edit projection'
    CHOICE_MOVE_VECTOR = 0
    CHOICE_ROTATE_VECTOR = 1
    CHOICE_ICOSAHEDRAL = 2
    CHOICE_DIHEDRAL = 3
    CHOICE_TETRAHEDRAL = 4
    CHOICE_LABEL = {
        CHOICE_MOVE_VECTOR: 'move source vector to target vector',
        CHOICE_ROTATE_VECTOR: 'rotate around vector',
        CHOICE_ICOSAHEDRAL: 'convert between icosahedral symmetries',
        CHOICE_DIHEDRAL: 'convert between dihedral symmetries',
        CHOICE_TETRAHEDRAL: 'convert between tetrahedral symmetries',
        }
    ICOSAHEDRAL_SYM_NAME = {
        SYM_I222: SCIPION_SYM_NAME[SYM_I222],
        SYM_I222r: SCIPION_SYM_NAME[SYM_I222r],
        SYM_In25: SCIPION_SYM_NAME[SYM_In25],
        SYM_In25r: SCIPION_SYM_NAME[SYM_In25r],
        SYM_I2n3: SCIPION_SYM_NAME[SYM_I2n3],
        SYM_I2n3r: SCIPION_SYM_NAME[SYM_I2n3r],
        SYM_I2n5: SCIPION_SYM_NAME[SYM_I2n5],
        SYM_I2n5r: SCIPION_SYM_NAME[SYM_I2n5r]
        }

    DIHEDRAL_SYM_NAME = {
        SYM_DIHEDRAL_X: SCIPION_SYM_NAME[SYM_DIHEDRAL_X],
        SYM_DIHEDRAL_Y: SCIPION_SYM_NAME[SYM_DIHEDRAL_Y]
        }

    TETRAHEDRAL_SYM_NAME = {
        SYM_TETRAHEDRAL_222: SCIPION_SYM_NAME[SYM_TETRAHEDRAL_222],
        SYM_TETRAHEDRAL_Z3: SCIPION_SYM_NAME[SYM_TETRAHEDRAL_Z3],
        SYM_TETRAHEDRAL_Z3R: SCIPION_SYM_NAME[SYM_TETRAHEDRAL_Z3R]
        }

    def _defineParams(self, form):
        """
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """
        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam,
                      pointerClass='SetOfParticles',
                      label='Set of particles to edit',
                      help='Set which items will be modified.'
                           'Must be a set of particles with transformations.')
        form.addParam(
            'operation', params.EnumParam,
            choices=[self.CHOICE_LABEL[self.CHOICE_MOVE_VECTOR],
                     self.CHOICE_LABEL[self.CHOICE_ROTATE_VECTOR],
                     self.CHOICE_LABEL[self.CHOICE_ICOSAHEDRAL],
                     self.CHOICE_LABEL[self.CHOICE_DIHEDRAL],
                     self.CHOICE_LABEL[self.CHOICE_TETRAHEDRAL],
                     ],
            default=self.CHOICE_ICOSAHEDRAL,
            label="Select operation",
            help="Select operation to be performed in the set.\n"
                 " *rotate to vector* modifies the alignment matrix "
                 " so a reconstruction made from the images produces a "
                 " rotated reconstruction\n"
                 " *convert between icosahedral symmetry* will modify "
                 " the alignment matrix so the reconstruction will have "
                 " the new symmetry provided\n"
                 " *convert between dihedral symmetry* will modify "
                 " the alignment matrix so the reconstruction will have "
                 " the new symmetry provided\n"
                 " *convert between tetrahedral symmetry* will modify "
                 " the alignment matrix so the reconstruction will have "
                 " the new symmetry provided\n"
            )
        # move vector
        line = form.addLine(
            'source vector',
            condition="operation==%d" % self.CHOICE_MOVE_VECTOR,
            help='Vector will be rotated until overlaps target vector')
        line.addParam('xs', params.FloatParam, default=0,
                      condition="operation==%d" % self.CHOICE_MOVE_VECTOR,
                      label="x",
                      help="X Coordinate ")
        line.addParam('ys', params.FloatParam, default=0,
                      condition="operation==%d" % self.CHOICE_MOVE_VECTOR,
                      label="y",
                      help="Y dim ")
        line.addParam('zs', params.FloatParam, default=1,
                      condition="operation==%d" % self.CHOICE_MOVE_VECTOR,
                      label="z",
                      help="Z dim ")

        line = form.addLine(
            'target vector',
            condition="operation==%d" % self.CHOICE_MOVE_VECTOR,
            help='source vector will be rotated until '
                 'overlaps target vector')
        line.addParam('xt', params.FloatParam, default=1,
                      condition="operation==%d" % self.CHOICE_MOVE_VECTOR,
                      label="x",
                      help="X Coordinate ")
        line.addParam(
            'yt', params.FloatParam, default=0,
            condition="operation==%d" % self.CHOICE_MOVE_VECTOR,
            label="y",
            help="Y dim ")
        line.addParam('zt', params.FloatParam, default=0,
                      condition="operation==%d" % self.CHOICE_MOVE_VECTOR,
                      label="z",
                      help="Z dim ")
        # rotate vector
        line = form.addLine(
            'rotate around vector',
            condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
            help='rotate around this vector by "angle" angles')
        line.addParam('x', params.FloatParam, default=0,
                      condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="x",
                      help="X coordinate ")
        line.addParam('y', params.FloatParam, default=0,
                      condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="y",
                      help="Y coordinate ")
        line.addParam('z', params.FloatParam, default=1,
                      condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="z",
                      help="Z coordinate")
        line.addParam('angle', params.FloatParam, default=0,
                      condition="operation==%d" % self.CHOICE_ROTATE_VECTOR,
                      label="angle (degrees)",
                      help="rotation angle ")
        # Icosahedral symmetry
        line = form.addLine(
            'Icosahedral symmetry',
            condition="operation==%d" % self.CHOICE_ICOSAHEDRAL,
            help='convert between icosahedral symmetries')
        line.addParam(
            "originSymmetryGroupI", params.EnumParam,
            condition="operation==%d" % self.CHOICE_ICOSAHEDRAL,
            label="origin",
            help="Source symmetry. ",
            choices=list(self.ICOSAHEDRAL_SYM_NAME.values()),
            default=SYM_I222r - SYM_I222,
            )
        line.addParam(
            "targetSymmetryGroupI", params.EnumParam,
            condition="operation==%d" % self.CHOICE_ICOSAHEDRAL,
            label="target",
            help="Target symmetry. ",
            choices=list(self.ICOSAHEDRAL_SYM_NAME.values()),
            default=SYM_I222 - SYM_I222
            )
        # Dihedral symmetry
        line = form.addLine(
            'Dihedral symmetry',
            condition="operation==%d" % self.CHOICE_DIHEDRAL,
            help='convert between dihedral symmetries')
        line.addParam(
            "originSymmetryGroupD", params.EnumParam,
            condition="operation==%d" % self.CHOICE_DIHEDRAL,
            label="origin",
            help="Source symmetry. ",
            choices=list(self.DIHEDRAL_SYM_NAME.values()),
            default=SYM_DIHEDRAL_X - SYM_DIHEDRAL_X,
            )
        line.addParam(
            "targetSymmetryGroupD", params.EnumParam,
            condition="operation==%d" % self.CHOICE_DIHEDRAL,
            label="target",
            help="Target symmetry. ",
            choices=list(self.DIHEDRAL_SYM_NAME.values()),
            default=SYM_DIHEDRAL_Y - SYM_DIHEDRAL_X,
            )

        # Tetrahedral symmetry
        line = form.addLine(
            'Tetrahedral symmetry',
            condition="operation==%d" % self.CHOICE_TETRAHEDRAL,
            help='convert between dihedral symmetries'
            )
        line.addParam(
            "originSymmetryGroupT", params.EnumParam,
            condition="operation==%d" % self.CHOICE_TETRAHEDRAL,
            label="origin",
            help="Source symmetry. ",
            choices=list(self.TETRAHEDRAL_SYM_NAME.values()),
            default=SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222,
            )
        line.addParam(
            "targetSymmetryGroupT", params.EnumParam,
            condition="operation==%d" % self.CHOICE_TETRAHEDRAL,
            label="target",
            help="Target symmetry. ",
            choices=list(self.TETRAHEDRAL_SYM_NAME.values()),
            default=SYM_TETRAHEDRAL_222 - SYM_TETRAHEDRAL_222,
            )

    def _insertAllSteps(self):
        operation = self.operation.get()
        if operation == self.CHOICE_MOVE_VECTOR:
            self._insertFunctionStep('rotateStep')
        elif operation == self.CHOICE_ROTATE_VECTOR:
            self._insertFunctionStep('rotateVectorStep')
        elif operation == self.CHOICE_ICOSAHEDRAL:
            self._insertFunctionStep('rotateIcosaStep')
        elif operation == self.CHOICE_DIHEDRAL:
            self._insertFunctionStep('rotateDiStep')
        elif operation == self.CHOICE_TETRAHEDRAL:
            self._insertFunctionStep('rotateTetraStep')

    def rotateTetraStep(self):
        """ Compute rotation matrix between tetrahedral symmetries
            and apply it to the projection directions
        """
        from pwem.convert.symmetry import Tetrahedral
        origSym = self.originSymmetryGroupT.get() + SYM_TETRAHEDRAL_222
        targetSym = self.targetSymmetryGroupT.get() + SYM_TETRAHEDRAL_222
        print("origSym", origSym)
        print("targetSym", targetSym)
        tetrahedral = Tetrahedral(sym=origSym)
        matrix = tetrahedral.coordinateSystemTransform(origSym, targetSym)

        # convert to numpy array
        # and add extra row
        matrix = np.array(matrix)
        matrix = np.append(matrix, [[0, 0, 0, 1]], axis=0)

        inputSet = self.inputSet.get()
        modifiedSet = inputSet.createCopy(self._getExtraPath(), copyInfo=True)
        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            transformation = item.getTransform()
            transformation.composeTransform(matrix)
            modifiedSet.append(item)
        self.createOutput(self.inputSet, modifiedSet)

    def rotateDiStep(self):
        """ Compute rotation matrix from one dihedral
        symmetry to another and apply it to the projection directions
        """
        from pwem.convert.symmetry import Dihedral
        origSym = self.originSymmetryGroupD.get() + SYM_DIHEDRAL_X
        targetSym = self.targetSymmetryGroupD.get() + SYM_DIHEDRAL_X
        # n is a required parameter not used in this case
        dihedral = Dihedral(sym=origSym, n=1)
        matrix = dihedral.coordinateSystemTransform(origSym, targetSym)

        # convert to numpy array
        # and add extra row
        matrix = np.array(matrix)
        matrix = np.append(matrix, [[0, 0, 0, 1]], axis=0)

        inputSet = self.inputSet.get()
        modifiedSet = inputSet.createCopy(self._getExtraPath(), copyInfo=True)
        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            transformation = item.getTransform()
            transformation.composeTransform(matrix)
            modifiedSet.append(item)
        self.createOutput(self.inputSet, modifiedSet)

    def rotateIcosaStep(self):
        """
        Compute rotation matrix from one icosahedral
        symmetry to another and apply it
        to projection directions. Reconstructed volume
        should change symmetry
        """
        from pwem.convert.symmetry import Icosahedron
        origSym = self.originSymmetryGroupI.get() + SYM_I222
        targetSym = self.targetSymmetryGroupI.get() + SYM_I222
        ico = Icosahedron(orientation=SCIPION_SYM_NAME[origSym][1:])
        matrix = ico.coordinateSystemTransform(SCIPION_SYM_NAME[origSym][1:],
                                               SCIPION_SYM_NAME[targetSym][1:])
        # convert to numpy array
        # and add extra row
        matrix = np.array(matrix)
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
        v_target = v_target / np.linalg.norm(v_target)
        matrix = rotation_matrix(angle_between_vectors(v_source, v_target),
                                 vector_product(v_source, v_target))
        print("rotateStep:matrix", matrix)
        inputSet = self.inputSet.get()
        modifiedSet = inputSet.createCopy(self._getExtraPath(), copyInfo=True)
        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            transformation = item.getTransform()
            transformation.composeTransform(matrix)
            modifiedSet.append(item)

        self.createOutput(self.inputSet, modifiedSet)

    def rotateVectorStep(self):
        """
        Compute rotation matrix around user provided vector
         (x,y,z) and rotate "angle" degrees the projection directions.
          Reconstructed volume should rotate around axis
        """
        v_source = np.array([self.x.get(), self.y.get(), self.z.get()])
        v_source = v_source / np.linalg.norm(v_source)
        angle = np.radians(self.angle.get())
        matrix = rotation_matrix(angle, v_source)
        print("matrix_rot_vector", matrix)
        inputSet = self.inputSet.get()
        modifiedSet = inputSet.createCopy(self._getExtraPath(), copyInfo=True)
        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            transformation = item.getTransform()
            transformation.composeTransform(matrix)
            modifiedSet.append(item)

        self.createOutput(self.inputSet, modifiedSet)

    def _validate(self):
        errors = []
        inputSet = self.inputSet.get()
        if not isinstance(inputSet, SetOfParticles):
            errors.append("The input data set is not a set of particles")
        if not inputSet.hasAlignmentProj():
            errors.append("The input data set does not have alignment 3D")
        return errors
