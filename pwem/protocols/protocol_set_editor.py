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
from pwem.protocols.protocol_sets import ProtSets
from pwem.objects.data import SetOfParticles

import numpy as np


class ProtSetEditor(ProtSets):
    """
    Protocol to edit attributes of all the items of a set using a formula.
    This could be useful for editing some values in the set. Use this
    protocol with extreme care, you can easily produce a set that is
    not consistent.
    """
    _label = 'edit set'

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
        # formula for main items
        form.addParam('formula', params.StringParam, label="Formula",
                      help='Any python code compatible with eval, where item represents each of '
                           'the elements of the set. E.g.: item._resolution.set(item._resolution.get() +1).'
                           'You could also use modules like "import numpy;  item._resolution .... "')

        # formula for subitems items
        form.addParam('subitemsformula', params.StringParam, label="Sub item formula",
                      help='Only valid for complex sets with subelements: Classes -> Particles, TiltSeries -> TiltImage, ..'
                           'Any python code compatible with eval, where item represents each of '
                           'the SUBELEMENTS each ITEM. E.g.: subitem.tiltAngle.set(subitem.tiltAngle.get() +1).'
                           'You could also use modules like "import numpy;  subitem._resolution .... "')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.formulaStep)

    def formulaStep(self):
        """
        Goes through all items in the input set and applies the formula to each of them using exec.
        Complex python code could be run separating lines with ;  To use numpy you could do
        import numpy; item._resolution.set(numpy.random.randint(10))
        """

        inputSet = self.inputSet.get()
        modifiedSet = inputSet.createCopy(self._getExtraPath(), copyInfo=True)

        for sourceItem in inputSet.iterItems():

            # item = sourceItem
            # exec(self.formula.get())
            updateItemCallBack = self.updateItem if self.formula.get() else None
            updateSubElemCallBack = self.updateSubElement if self.subitemsformula.get() else None

            self._append(modifiedSet, sourceItem, itemUpdateCallback=updateItemCallBack, subElemUpdateCallback=updateSubElemCallBack)

        self.createOutput(self.inputSet, modifiedSet)

    def updateItem(self, item):

        exec(self.formula.get())
    def updateSubElement(self, subitem):


        exec(self.subitemsformula.get())

    def createOutput(self, inputSet, modifiedSet):
        """ Save the output set."""
        outputArgs = {inputSet.getExtended(): modifiedSet}
        self._defineOutputs(**outputArgs)

    def _summary(self):

        return ["The this formula (%s) is/was applied to the items of the input set." % self.formula.get()]


