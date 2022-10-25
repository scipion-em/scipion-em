# **************************************************************************
# *
# * Authors:     Daniel Marchan [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import pwem.objects as emobj
from pwem.protocols import EMProtocol
import logging
logger = logging.getLogger(__file__)


OUTPUT_CLASSES = "outputClasses"


class ProtGoodClassesExtractor(EMProtocol):
    """ Extracts items from a SetOfClasses based on the IDs of the given good averages/classes
    """

    _label = "good classes selector"
    outputsToDefine = {}

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputClasses', params.PointerParam,
                      pointerClass='SetOfClasses',
                      label='Input classes',
                      help='Set of classes to extract items from.')

        form.addParam('inputGoodClasses', params.PointerParam,
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      label='Good references',
                      help='Set of good reference to stay with from the inputClasses.')

    # -------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        """ Insert all steps
            """
        self.goodClassesIDs = []
        self._insertFunctionStep(self.selectGoodClasses)
        self._insertFunctionStep(self.createOutputStep)

    def selectGoodClasses(self):
        """ Select only the good Classes from the Averages
        """
        #print(set(self.inputGoodClasses.get().getUniqueValues('filename')))
        self.goodClassesIDs = self.inputGoodClasses.get().getIdSet()
        self.info('Good classes IDs:')
        self.info(self.goodClassesIDs)
        # Change to list of filenames
        #for clazz in self.inputGoodClasses.get().iterItems(orderBy="id", direction="ASC"):
        #    print(clazz.getRepresentative()._filename)
        #    self.info("Class filename selected: %s" % clazz.getObjName())
        #    self.goodClassesIDs.append(clazz.getObjId())

    def createOutputStep(self):
        """ Create output
                    """
        inputClasses = self.inputClasses.get()
        outputClasses = emobj.SetOfClasses2D.create(self._getExtraPath())
        outputClasses.copyInfo(inputClasses)
        outputClasses.appendFromClasses(inputClasses, filterClassFunc=self._appendClass)

        self.outputsToDefine[OUTPUT_CLASSES] = outputClasses
        self._defineOutputs(**self.outputsToDefine)

    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        errors = []
        return errors

# --------------------------- UTILS functions -----------------------------
    def _appendClass(self, item):
        return False if not item.getObjId() in self.goodClassesIDs else True