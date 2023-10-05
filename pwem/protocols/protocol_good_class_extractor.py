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
import pyworkflow.utils as pwutils
from pwem.protocols import EMProtocol
from pwem.objects import SetOfParticles
import logging
logger = logging.getLogger(__file__)


OUTPUT_PARTICLES = "outputParticles"


class ProtGoodClassesExtractor(EMProtocol):
    """ Extracts items from a SetOfClasses based on the IDs of the given good averages/classes
    """

    _label = "good classes selector"
    outputsToDefine = {}
    _possibleOutputs = {OUTPUT_PARTICLES: SetOfParticles}

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
        self._insertFunctionStep(self.extractElements)
        # self._insertFunctionStep(self.createOutputStep)

    def extractElements(self):
        # For each class (order by number of items)
        for clazz in self.inputClasses.get().iterItems(orderBy="_size", direction="DESC"):
            if clazz.getObjId() in self.goodClassesIDs:
                self._extractElementFromClass(clazz)

        output = self._getOutputSet()
        output.write()
        self._store(output)

    def _extractElementFromClass(self, clazz):
        output = self._getOutputSet()
        # Go through all items and append them
        for image in clazz:
            newImage = image.clone()
            output.append(newImage)

    def selectGoodClasses(self):
        """ Select only the good Classes from the Averages
        """
        self.goodClassesIDs = self.inputGoodClasses.get().getIdSet()
        self.info('Good classes IDs:')
        self.info(self.goodClassesIDs)
        # Change to list of filenames
        #print(set(self.inputGoodClasses.get().getUniqueValues('filename')))
        #for clazz in self.inputGoodClasses.get().iterItems(orderBy="id", direction="ASC"):
        #    print(clazz.getRepresentative()._filename)
        #    self.info("Class filename selected: %s" % clazz.getObjName())
        #    self.goodClassesIDs.append(clazz.getObjId())

    def _getOutputSet(self):
        """ Returns the output set creating it if not yet done"""
        # If output not created yet
        if not hasattr(self, OUTPUT_PARTICLES):
            outputSet = None
            self.info("Creating set from images.")
            outputSet = createSetFromImages(self.inputClasses.get(), self._getPath())
            self._defineOutputs(**{OUTPUT_PARTICLES: outputSet})

        return getattr(self, OUTPUT_PARTICLES)


    def closeOutputStep(self):
        self._closeOutputSet()

    # def createOutputStep(self):
        # """ Create output
        #             """
        # inputClasses = self.inputClasses.get()
        # outputClasses = SetOfClasses2D.create(self._getExtraPath())
        # outputClasses.copyInfo(inputClasses)
        # outputClasses.appendFromClasses(inputClasses, filterClassFunc=self._appendClass)
        #
        # outputParticles = createSetFromImages(outputClasses, self._getPath())
        #
        # self.outputsToDefine[OUTPUT_CLASSES] = outputClasses
        # self.outputsToDefine[OUTPUT_PARTICLES] = outputParticles
        #
        # self._defineOutputs(**self.outputsToDefine)

    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        errors = []
        return errors

# --------------------------- UTILS functions -----------------------------
    def _appendClass(self, item):
        return False if not item.getObjId() in self.goodClassesIDs else True

def createSetFromImages(classesSet, path):
    """ Creates a the corresponding set from the images of a set of classes"""
    images = classesSet.getImages()
    # need to instantiate the right set based on the Items.
    setClass = None
    logger.debug("Creating an image set from %s: type %s" % (classesSet.__class__, images))
    setClass = pwutils.Config.getDomain().getObjects()[images.__class__.__name__]
    set = setClass.create(outputPath=path)
    set.copyInfo(images)

    return set

