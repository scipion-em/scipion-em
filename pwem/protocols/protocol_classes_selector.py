# **************************************************************************
# *
# * Authors:     Pablo Conesa [1]
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

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params

import pwem.objects as emobj
from pwem.protocols import EMProtocol
import logging
logger = logging.getLogger(__file__)

OUTPUT_NAME = "output"


class ProtClassesSelector(EMProtocol):
    """ Extracts items from a SetOfClasses based on number of items assigned to the class
    """

    _label = "numeric classes extractor"
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputClasses', params.PointerParam,
                      pointerClass='SetOfClasses',
                      label='Input classes',
                      help='Set of classes to extract items from.')
        form.addParam('extractRepresentative', params.BooleanParam,
                      label='Extract representative?',
                      help='Set to true if you want to extract the image that represents the class otherwise all '
                           'items supporting the class will be extracted.',
                      default = False)
        form.addParam('firstNElements', params.IntParam,
                      label='Biggest N classes',
                      help='Will take the N biggest classes. Those having most element supporting them.',
                      default=5)


# -------------------------- INSERT steps functions ---------------------------

    def _insertAllSteps(self):
        """ Insert all steps
        """
        self._insertFunctionStep(self.extractElements.__name__)

    def extractElements(self):

        count= 0
        # For each class (order by number of items)
        for clazz in self.inputClasses.get().iterItems(orderBy="_size", direction="DESC"):

            # Increment the count
            count += 1

            # If haven't reached the Nelements
            if count <= self.firstNElements.get():

                print(clazz.getSize())
                self._extractElementFromClass(clazz)

        output = self._getOutputSet()
        output.write()
        self._store(output)

    def _extractElementFromClass(self, clazz):

        output = self._getOutputSet()

        # If taking representatives
        if self.extractRepresentative.get():
            rep = clazz.getRepresentative().clone()
            output.append(rep)
        else:
            # Go through all items and append them
            for image in clazz:
                newImage = image.clone()
                output.append(newImage)

    def _getOutputSet(self):
        """ Returns the output set creating it if not yet done"""

        # If output not created yet
        if not hasattr(self, OUTPUT_NAME):
            outputSet = None

            # Create a set ot items or a set of representatives
            if self.extractRepresentative.get():

                self.info("Creating set from representatives.")
                outputSet = createSetFromRepresentative(self.inputClasses.get(), self._getPath())
            else:
                self.info("Creating set from images.")
                outputSet = createSetFromImages(self.inputClasses.get(), self._getPath())

            self._defineOutputs(**{OUTPUT_NAME:outputSet})

        return getattr(self, OUTPUT_NAME)


    def _summary(self):
        summary = []
        summary.append("%s option selected." % "representative/s" if self.extractRepresentative else "items")
        summary.append("%s biggest classes to be extracted from." % self.firstNElements.get())
        return summary

    def _methods(self):
        extractedElement = "representative/s" if self.extractRepresentative else "items"
        methods = ["The %s of the %s biggest classes were selected from %s." % (extractedElement, self.firstNElements.get(), self.getObjectTag(self.inputClasses))]
        return methods

    def _validate(self):
        errors = []
        return errors

# Helpers
def createSetFromRepresentative(classesSet, path):
    """ Creates a teh corresponding set from the representative of a set of classes"""

    rep_type = classesSet.REP_TYPE

    # need to instantiate the right set based on the representative.
    # Going from the REP_TYPE to the right set is not straight forward. e.g:
    # SetOfClasses2D --> REP_TYPE = Particle --> Set for the rep type == SetOfAverages
    # So we need to make do some manual checks. TODO: have factory methods in the SetOfClasses for the representatives?
    setClass = None
    logger.debug("Creating representative set from %s: type %s" % (classesSet.__class__, rep_type))
    if rep_type.__name__ ==  emobj.Particle.__name__:
        setClass = emobj.SetOfAverages
    else:
        setClass = pwutils.Config.getDomain().getObjects()["SetOf%ss" % rep_type.__name__]
    set = setClass.create(outputPath=path)

    imgs = classesSet.getImages()
    set.copyInfo(imgs)

    return set

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



