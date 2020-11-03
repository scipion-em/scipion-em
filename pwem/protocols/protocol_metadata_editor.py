# **************************************************************************
# *
# * Authors:     Pablo Conesa(pconesa@cnb.csic.es)
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


class ProtMetadataEditor(EMProtocol):
    """
    Protocol to edit attributes of all the items of a set using a formula.
    This could be useful for corrupting your data for testing purposes or
    editing some values in the set that make sense to do it. Use this with
    extreme care, you can easily ruin your processing.
    """
    _label = 'metadata editor'

    def _defineParams(self, form):
        """
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """

        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam, pointerClass='EMSet',
                      label='Set to edit',
                      help='Set which items will be modified.')
        form.addParam('formula', params.StringParam, label="Formula", important=True,
                      help='A python code compatible with eval, where item represents each of '
                           'the elements of the set. E.g.: item._resolution.set(item._resolution.get() +1).'
                           'You could also use modules like "import numpy;  item._resolution .... "')

    def _insertAllSteps(self):
        self._insertFunctionStep('editItemsStep')

    def editItemsStep(self):
        """
        Goes through all items in the input set and applies the formula to each of them using exec.
        Complex python code could be run separating lines with ;  To use numpy you could do
        import numpy; item._resolution.set(numpy.random.randint(10))
        """
        inputSet = self.inputSet.get()

        modifiedSet = inputSet.create(self._getExtraPath())
        modifiedSet.copyInfo(inputSet)

        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            exec(self.formula.get())
            modifiedSet.append(item)

        outputArgs = {self.inputSet.getExtended(): modifiedSet}
        self._defineOutputs(**outputArgs)
