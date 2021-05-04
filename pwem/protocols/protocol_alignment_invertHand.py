# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

from .protocol_2d import ProtAlign2D
from pyworkflow.utils import ProgressBar


class ProtAlignmentInvertHand(ProtAlign2D):
    """ Modify the transformation matrix of a set of particles
    So that the handedness changes
    """
    _label = 'invert hand'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', params.PointerParam, pointerClass='SetOfParticles',
                      label='Input particles',
                      help='Select the particles that you want to update the new alignment.')
        form.addParallelSection(threads=0, mpi=0)

# --------------------------- INSERT steps functions --------------------------

    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        inputParticles = self.inputParticles.get()
        outputParticles = self._createSetOfParticles()
        outputParticles.copyInfo(inputParticles)
        outputParticles.setAlignment(inputParticles.getAlignment())
        total = len(inputParticles)
        width = 40
        step = total // width

        # print a progressbar since this make take for ever
        pb = ProgressBar(total=total, width=width, fmt=ProgressBar.NOBAR)
        for counter, particle in enumerate(inputParticles):
            if counter % step == 0:
                pb.update(counter)
            t = particle.getTransform()
            m = t.getMatrix()
            m[0, 2] *= -1.
            m[1, 2] *= -1.
            m[2, 1] *= -1.
            m[2, 0] *= -1.
            # particle.setTransform(t.setMatrix(m))
            outputParticles.append(particle)
        pb.finish()

        self._defineOutputs(outputParticles=outputParticles)
        self._defineSourceRelation(self.inputParticles, outputParticles)

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputParticles'):
            methods.append("Output particles not ready yet.")
        return methods

    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        # check that input set of aligned particles do have 2D alignment
        errors = []
        inputAlignmentSet = self.inputParticles.get()
        if not inputAlignmentSet.hasAlignment():
            errors.append("Input alignment set should contains some kind of alignment (2D, 3D or Projection).")
        else:
            # Just for consistency, check that the particles really contains Transform object
            first = inputAlignmentSet.getFirstItem()
            alignment = first.getTransform()
            if alignment is None:
                errors.append('Inconsistency detected in *Input alignment* !!!')
                errors.append('It has alignment: _%s_, but the alignment is missing!!!' %
                              inputAlignmentSet.getAlignment())
            
        # Add some errors if input is not valid
        return errors
