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
from pwem.protocols import EMProtocol


class ProtSetFilter(EMProtocol):
    """
    Protocol to filter sets based on its attributes through an expresion that should evaluate to true
    or false.
    """
    _label = 'filter set'

    def _defineParams(self, form):
        """
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """

        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam, pointerClass='EMSet',
                      label='Set to filter',
                      help='Set which items will be filtered.')
        form.addParam('formula', params.StringParam, label="Passing formula", important=True,
                      help='A python code compatible with eval that should evaluate to True or False, where item represents each of '
                           'the elements of the set. E.g.: item._resolution.get() < 4).'
                           'You could also use modules like "import numpy;  item._resolution .... "')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.filterSetStep)

    def filterSetStep(self):
        """
        Goes through all items in the input set and applies the formula to each of them using exec.
        Complex python code could be run separating lines with ;  To use numpy you could do
        import numpy; item._resolution.set(numpy.random.randint(10))
        If result is True, item will be transferred to the output set
        """
        inputSet = self.inputSet.get()

        modifiedSet = inputSet.create(self._getExtraPath())
        modifiedSet.copyInfo(inputSet)

        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            exec("item.setEnabled(%s)"% self.formula.get())
            if item.isEnabled():
                modifiedSet.append(item)

        outputArgs = {self.inputSet.getExtended(): modifiedSet}
        self._defineOutputs(**outputArgs)
