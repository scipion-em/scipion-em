# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
This module contains protocols related to Set operations such us:
- subsets
- unions
- split
... etc
"""

import random
import sys

import pyworkflow.protocol as pwprot
from pyworkflow.object import Object,Float, Integer, String

from pwem.protocols import EMProtocol
from pwem.objects import Volume, EMSet, SetOfClasses, SetOfData
from pyworkflow.utils import ProgressBar, getListFromRangeString
from pwem.constants import ID_COLUMN, ID_ATTRIBUTE


class ProtSets(EMProtocol):
    """ Base class for all protocols related to subsets. """

    def _append(self, outputSet, item, sourceItem=None, itemUpdateCallback=None, subElemUpdateCallback=None):
        """ Add an item to the outputSet.
        If the item is a new copy of sourceItem(case of the join sets),
        then use a sourceItem since item lost the information related with the
        mapper

        :param itemUpdateCallback: callback receiving the item to apply any operation (optional)
        :param subElemUpdateCallback: callback receiving the sub-element to apply any operation (optional)
        """
        subElemList = []
        if sourceItem is None:
            sourceItem = item
        if isinstance(item, EMSet):
            for subElem in sourceItem.iterItems():
                # We need to create a clone because all items have a same _objId
                clon = subElem.clone(copyEnable=True)

                # Update the sub-element if callback is passed
                if subElemUpdateCallback:
                    subElemUpdateCallback(clon)

                subElemList.append(clon)

        # Update the main item if callback is passed
        if itemUpdateCallback:
            itemUpdateCallback(item)

        outputSet.append(item)
        if subElemList:
            for subElem in subElemList:
                item.append(subElem)
            # When adding sub-elements, some item "summary" properties may be updated: e.g. TiltSeries anglesCount.
            # Need to persist them.
            outputSet.update(item)


class ProtUnionSet(ProtSets):
    """ Protocol to join two or more sets of images.
    This protocol allows to select two or more set of images
    and will produce another set joining all elements of the 
    selected sets. It will validate that all sets are of the
    same type of elements (Micrographs, Particles or Volumes) 
    """
    _label = 'join sets'
    TYPE_CTF = 'CTFs'
    TYPE_VOLUME='Volumes'
    TYPE_VOLUME_INDEX = 3

    _unionTypes = ['Particles',
                   'Micrographs',
                   TYPE_CTF,
                   TYPE_VOLUME,
                   'Averages',
                   'All']

    def __init__(self, **kwargs):
        ProtSets.__init__(self, **kwargs)

        # We need to trace the changes of 'inputType' to
        # dynamically modify the property of pointerClass
        # of the 'inputSets' parameter
        def onChangeInputType():
            inputText = self.getEnumText('inputType')

            if inputText == 'All':
                pointerClass = 'EMSet'
            # elif inputText == 'CTFs + Micrographs':
            #     pointerClass = 'SetOfCTF'
            else:
                pointerClass = 'SetOf' + inputText
            # For relatively small set we usually want to include
            # the single element type, this will allow, for example
            # to union SetOfVolumes and Volumes in the final set
            if inputText in [self.TYPE_VOLUME]:
                pointerClass += ',%s' % inputText[:-1]  # remove last 's'
            elif inputText in [self.TYPE_CTF]:
                # remove last 's'
                pointerClass = '%s,CTFModel' % pointerClass[:-1]

            self.inputSetsParam.setPointerClass(pointerClass)

        self.inputType.trace(onChangeInputType)

    # -------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputType', pwprot.params.EnumParam,
                      choices=self._unionTypes, default=5,  # All
                      label='Input type:',
                      help='Select the type of objects that you want to union.\n'
                           'Special case All will allow you to select any type.')
        self.inputSetsParam = form.addParam('inputSets', pwprot.params.MultiPointerParam,
                                            label="Input set", important=True,
                                            pointerClass='EMSet', minNumObjects=2, maxNumObjects=0,
                                            help='Select two or more sets (of micrographs, particles,'
                                                 ' volumes, etc.) to be united. If you select 3 sets '
                                                 'with 100, 200, 200 elements, the final set will '
                                                 'contain a total of 500 elements.')
        form.addParam('ignoreDuplicates', pwprot.params.BooleanParam,
                      default=False,
                      label='Remove duplicates?',
                      help='By default, duplicated items found (same ID) '
                           'within the input sets, will cause renumbering of all the '
                           'items ids in the output set. '
                           'This is the case for example when doing several '
                           'imports (which will cause ids overlapping) '
                           'but we really want to insert as new items in the '
                           'output. \n'
                           'On the other hand, items originated in a previous common '
                           'protocol (above in the workflow) might have identical items '
                           'and you would like to remove them. '
                           'Therefore, set this option to *Yes* to remove duplicates and keep only '
                           'one copy of the item (the first occurrence).')
        form.addParam('renumber', pwprot.params.BooleanParam, default=False,
                      expertLevel=pwprot.LEVEL_ADVANCED,
                      label="Force new ids",
                      help='Perform an automatic renumbering of ids to ensure all objects have unique ids. '
                           'This will mean new objects will not be associated to the old ones.')

        # TODO: See what kind of restrictions we add,
        # like "All sets should have the same sampling rate."

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def createOutputStep(self):

        set1 = self.inputSets[0].get()  # 1st set (we use it many times)

        # Read ClassName and create the corresponding EMSet (SetOfParticles...)
        try:
            if str(set1.getClassName()) is not Volume.__name__:
                outputSetFunction = getattr(self, "_create%s" % set1.getClassName())
            else:
                outputSetFunction = self._createSetOfVolumes
            outputSet = outputSetFunction()
        except Exception:
            outputSet = set1.createCopy(self._getPath())

        # Copy info from input sets (sampling rate, etc).
        if str(set1.getClassName()) is not Volume.__name__:
            outputSet.copyInfo(set1)  # all sets must have the same info as set1!
        else:
            outputSet.setSamplingRate(set1.getSamplingRate())

        # Renumber from the beginning if either the renumber option is selected
        # or we find duplicated ids in the sets
        cleanIds = not self.ignoreDuplicates.get() and self.duplicatedIds()

        # Warn in the log in case attributes will be lost
        allSetAttributes, commonAttrs = self.commonAttributes()
        warnings = self._getHeterogeneityWarning(allSetAttributes, commonAttrs)
        if warnings:
            self.info(warnings)

        # Always ignore non-common attributes
        ignoreExtraAttributes = True

        # Get the 1st level attributes to be used for the copyAttributes
        copyAttrs = list()
        for attr in commonAttrs:
            if "." not in attr:
                copyAttrs.append(attr)

        self.info("Common attributes to all sets are: %s" % copyAttrs)

        idsList = {}
        setNum = 0
        for itemSet in self.inputSets:
            setNum += 1
            if str(itemSet.get().getClassName()) is not Volume.__name__:
                for obj in itemSet.get():
                    objId = obj.getObjId()
                    if self.ignoreDuplicates.get():
                        if objId in idsList:
                            continue
                        idsList[objId] = objId
                    # This is always TRUE, if stable we could remove the "if" and the "else".
                    if ignoreExtraAttributes:
                        newObj = itemSet.get().ITEM_TYPE()
                        newObj.copyAttributes(obj, *copyAttrs)

                        self.cleanExtraAttributes(newObj, commonAttrs)
                        if not cleanIds or setNum == 1:
                            newObj.setObjId(objId)
                    else:
                        newObj = obj

                    if (cleanIds and setNum > 1) or self.renumber.get():
                        newObj.cleanObjId()

                    self._append(outputSet, newObj, sourceItem=obj)

            else:
                obj = itemSet.get()
                objId = obj.getObjId()
                if self.ignoreDuplicates.get():
                    if objId in idsList:
                        continue
                    idsList[objId] = objId
                newObj = obj
                if (cleanIds and setNum > 1) or self.renumber.get():
                    newObj.cleanObjId()
                outputSet.append(newObj)

        self._defineOutputs(outputSet=outputSet)
        for itemSet in self.inputSets:
            self._defineSourceRelation(itemSet, outputSet)

    # Overwrite SetOfCoordinates creation
    def _createSetOfCoordinates(self, suffix=''):
        coordSet = self.inputSets[0].get()
        micSet = coordSet.getMicrographs()
        return ProtSets._createSetOfCoordinates(self, micSet, suffix)

    def cleanExtraAttributes(self, obj, verifyAttrs, prefix=""):

        for attr, value in obj.getAttributesToStore():

            prefixedAttribute = prefix + attr

            if prefixedAttribute not in verifyAttrs:
                value._objDoStore = False
                self.info("%s will be lost." % attr)

            else:
                self.cleanExtraAttributes(value, verifyAttrs,
                                          prefixedAttribute + ".")

    # def getObjDict(self, includeClass=False, includeBasic=False):
    #     return super(ProtUnionSet, self).getObjDict(
    #         includeClass=includeClass, includeBasic=includeBasic)

    def duplicatedIds(self):
        """ Check if there are duplicated ids to renumber from
        the beginning. """
        usedIds = set()  # to keep track of the object ids we have already seen
        for item_pointer in self.inputSets:
            if str(item_pointer.get().getClassName()) is not Volume.__name__:
                for objIds in item_pointer.get().getIdSet():
                    if objIds in usedIds:
                        return True
                    else:
                        usedIds.add(objIds)
            else:
                objId = item_pointer.get().getObjId()
                if objId in usedIds:
                    return True
                else:
                    usedIds.add(objId)
        return False

    def getAllSetsAttributes(self):
        allSetsAttributes = list()
        for itemSet in self.inputSets:
            if str(itemSet.get().getClassName()) is not Volume.__name__:
                item = itemSet.get().getFirstItem()
            else:
                item = itemSet.get()
            attrs = set(item.getObjDict().keys())
            allSetsAttributes.append(attrs)

        return allSetsAttributes

    def commonAttributes(self):
        """ Compute the set of common attributes to all items within
        each input set. """
        commonAttrs = None
        allSetsAttributes = self.getAllSetsAttributes()

        for attrSet in allSetsAttributes:
            if commonAttrs is None:  # first time
                commonAttrs = attrSet
            else:
                commonAttrs = commonAttrs & attrSet

        return allSetsAttributes, list(commonAttrs)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        # Are all inputSets from the same class?
        classes = {x.get().getClassName() for x in self.inputSets}
        if len(classes) > 1:
            return ["All objects should have the same type.",
                    "Types of objects found: %s" % ", ".join(classes)]
        if issubclass(type(self.inputSets[0].get()), SetOfClasses):
            return ["Is not possible to join different sets of classes.\n"
                    "If you want to join different representative, extract them "
                    "with the viewer and them run this protocol with the "
                    "resulting averages."]

        # Validate attributes like sampling rate or dimensions
        return self._checkSetsCompatibility()

    def _checkSetsCompatibility(self):
        """ Check if all input sets have a minimum compatible attributes """
        # Attributes to check -> defined by the Set subclass type that are requested to be joined
        firstInput=self.inputSets[0].get()

        # Verify inputs are sets
        if not isinstance(firstInput, EMSet):
            return []

        attrs = firstInput.getCompatibilityDict()
        errors = []
        # For each attribute
        for key, attr in attrs.items():

            # Intentional: we need a default value not None, since some
            # attributes could return None as a valid value.
            refValue = '?'

            # For pointer to a set
            for setPointer in self.inputSets:

                # Get the set:
                inputSet = setPointer.get()

                # If the set has the attribute
                if not hasattr(inputSet, attr):
                    break

                # Get the attribute and "call it" --> final ().
                setValue = getattr(inputSet, attr)()

                if refValue == '?':
                    refValue = setValue
                else:
                    if refValue != setValue:
                        errors.append("There are different %s among the input"
                                      " sets: %s and %s" % (key, refValue, setValue))
                        break

        return errors

    def _warnings(self):
        """ Warn about loosing info. """

        # Get all attributes "map"
        allSetsAttributes, commonAttributes = self.commonAttributes()

        return self._getHeterogeneityWarning(allSetsAttributes, commonAttributes)

    def _getHeterogeneityWarning(self, allSetsAttributes, commonAttributes):

        warnings = []
        # Use a set
        commonAttributes = set(commonAttributes)

        # Go through all sets attributes
        for index, setAttributes in enumerate(allSetsAttributes):
            setAttributes = set(setAttributes)
            # Get the difference
            lostAttributes = setAttributes - commonAttributes

            if len(lostAttributes) != 0:
                warnings.append("Set #%d will loose following "
                                "attributes:" % index)
                for attr in lostAttributes:
                    warnings.append(attr)

        if len(warnings):
            warnings.append("Your input sets have different attributes. "
                            "We will keep only the common ones. This may "
                            "cause the lost of important data like CTF, "
                            "alignment information,...")

        return  warnings

    def _summary(self):
        if not hasattr(self, 'outputSet'):
            return ["Protocol has not finished yet."]
        else:
            return ["We have merged the following sets:",
                    ", ".join(x.get().getNameId() for x in self.inputSets)]

    def _methods(self):
        return self._summary()


class ProtSplitSet(ProtSets):
    """ Protocol to split a set in two or more subsets.
    """
    _label = 'split sets'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSet', pwprot.params.PointerParam,
                      pointerClass='EMSet',
                      label="Input set", important=True,
                      help='Select the set of elements (images, etc) that you '
                           'want to split.')

        form.addParam('numberOfSets', pwprot.params.IntParam, default=2,
                      label="Number of subsets",
                      help='Select how many subsets do you want to create.')

        form.addParam('randomize', pwprot.params.BooleanParam, default=False,
                      label="Randomize elements",
                      help='Put the elements at random in the different '
                           'subsets.')

    # Overwrite SetOfCoordinates creation
    def _createSetOfCoordinates(self, suffix=''):
        coordSet = self.inputSet.get()
        micSet = coordSet.getMicrographs()
        return ProtSets._createSetOfCoordinates(self, micSet, suffix)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        inputSet = self.inputSet.get()
        inputClassName = str(inputSet.getClassName())
        n = self.numberOfSets.get()
        # Create as many subsets as requested by the user
        try:
            outputSetFunction = getattr(self, "_create%s" % inputClassName)
            subsets = [outputSetFunction(suffix=str(i)) for i in range(1, n + 1)]
        except Exception:
            subsets = [inputSet.createCopy(self._getPath(), suffix=str(i))
                       for i in range(1, n + 1)]

        # Iterate over the elements in the input set and assign
        # to different subsets.
        elements = self.inputSet.get()

        ns = [len(elements) // n + (1 if i < len(elements) % n else 0)
              for i in range(n)]  # number of elements in each subset
        pos, i = 0, 0  # index of current subset and index of position inside it
        orderBy = 'RANDOM()' if self.randomize else 'id'

        for elem in elements.iterItems(orderBy=orderBy, direction='ASC'):
            if i >= ns[pos]:
                pos += 1
                i = 0
            self._append(subsets[pos], elem)
            i += 1

        key = 'output' + inputClassName.replace('SetOf', '') + '%02d'
        for i in range(1, n + 1):
            subset = subsets[i - 1]
            subset.copyInfo(inputSet)
            self._defineOutputs(**{key % i: subset})
            self._defineTransformRelation(inputSet, subset)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if self.inputSet.get().getSize() < self.numberOfSets:
            errors.append("The number of subsets requested is greater than")
            errors.append("the number of elements in the input set.")
        return errors

    def _summary(self):
        if not any(x.startswith('output') for x in dir(self)):
            return ["Protocol has not finished yet."]
        else:
            return ["We have split the set %s in %d sets." %
                    (self.inputSet.getName(), self.numberOfSets.get())]


class ProtSubSet(ProtSets):
    """    
    Create a set with the elements of an original set that are also
    referenced in another set.
    
    Usually there is a bigger set with all the elements, and a smaller
    one obtained from classification, cleaning, etc. The desired result
    is a set with the elements from the original set that are also present
    somehow in the smaller set (in the smaller set they may be downsampled
    or processed in some other way).
    
    Both sets should be of the same kind (micrographs, particles, volumes)
    or related (micrographs and CTFs for example).
    """
    _label = 'subset'
    SET_INTERSECTION = 0
    SET_DIFFERENCE = 1

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        add = form.addParam  # short notation
        add('inputFullSet', pwprot.params.PointerParam, pointerClass='EMSet',
            label="Full set of items", important=True,
            help='Even if the operation can be applied to two arbitrary sets,\n'
                 'the most common use-case is to retrieve a subset of\n'
                 'elements from an original full set.\n'
                 '*Note*: the elements of the resulting set will be the same\n'
                 'ones as this input set.')
        add('chooseAtRandom', pwprot.params.BooleanParam, default=False,
            label="Make random subset",
            help='Choose elements randomly form the full set.')
        add('nElements', pwprot.params.IntParam, default=2,
            condition='chooseAtRandom',
            label="Number of elements",
            help='How many elements will be taken from the full set.')
        add('selectIds', pwprot.params.BooleanParam, default=False,
            condition='not chooseAtRandom',
            label="Make a subset from specific IDs",
            help="Choose specific elements form the full set.")
        add('range', pwprot.params.NumericRangeParam,
            label="IDs range or list",
            condition='selectIds and not chooseAtRandom',
            allowsNull=True,
            help='Select the IDs that will be the subset.\n'
                 'You have several ways to specify the IDs.\n'
                 'Example: \n'
                 '"1,3,5-8,17-20" -> [1,3, 5, 6, 7, 8, 17, 18, 19, 20]\n')
        add('inputSubSet', pwprot.params.PointerParam,
            pointerClass='EMSet', condition='not (chooseAtRandom or selectIds)',
            label="Other set",
            allowsNull=True,
            help='The elements present in this set will be used to pick \n'
                 'elements from the input full set.     \n'
                 'This means that the output set will contain elements with \n'
                 'exact the same information of input full set.\n\n'
                 'Set operation: if _intersection_ is used,\n'
                 'elements that are both in input and other set\n'
                 'will be included. If _difference_, elements that\n'
                 'are in input but not in other will picked.')
        add('setOperation', pwprot.params.EnumParam,
            condition='not (chooseAtRandom or selectIds)',
            default=self.SET_INTERSECTION,
            choices=['intersection', 'difference'],
            display=pwprot.params.EnumParam.DISPLAY_HLIST,
            label='Set operation',
            help='Set operation: if _intersection_ is used,\n'
                 'elements that are both in input and other set\n'
                 'will be included. If _difference_, elements that\n'
                 'are in input but not in other will picked.')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        inputFullSet = self.inputFullSet.get()

        inputClassName = inputFullSet.getClassName()

        try:
            outputSetFunction = getattr(self, "_create%s" % inputClassName)
            outputSet = outputSetFunction()
        except Exception:
            outputSet = inputFullSet.createCopy(self._getPath())

        outputSet.copyInfo(inputFullSet)

        if self.chooseAtRandom or self.selectIds:
            if self.chooseAtRandom:
                # Get all ids form iput set
                self.info("Creating subset from random positions from input set.")
                ids = set(random.sample(list(inputFullSet.getIdSet()), self.nElements.get()))
            else:
                self.info("Creating subset by range: %s" % self.range)
                ids = set(getListFromRangeString(self.range.get()))
        else:
            # Get the ids from both sets
            fullSetIds = inputFullSet.getIdSet()
            smallSetIds = self.inputSubSet.get().getIdSet()

            # The function to include an element or not
            # depends on the set operation
            # if it is 'intersection' we want that item is not None (found)
            # if it is 'difference' we want that item is None
            # (not found, different)
            if self.setOperation == self.SET_INTERSECTION:
                ids = fullSetIds.intersection(smallSetIds)
            else:
                ids = fullSetIds.difference(smallSetIds)

        progress = None
        nElements = len(ids)

        if nElements > 100000:  # show progressBar for large sets
            progress = ProgressBar(total=nElements, fmt=ProgressBar.NOBAR)
            progress.start()
            sys.stdout.flush()
            step = max(25000, nElements // 25000)

        i = 0

        for elem in inputFullSet.iterItems():
            if elem.getObjId() in ids:
                i += 1
                if progress and i % step == 0:
                    progress.update(i)
                self._append(outputSet, elem)

        if progress:
            progress.finish(printNewLine=True)

        if outputSet.getSize():
            key = 'output' + inputClassName.replace('SetOf', '')
            self._defineOutputs(**{key: outputSet})
            self._defineTransformRelation(inputFullSet, outputSet)
            if not (self.chooseAtRandom.get() or self.selectIds.get()):
                self._defineSourceRelation(self.inputSubSet, outputSet)
        else:
            self.summaryVar.set('Output was not generated. Resulting set '
                                'was EMPTY!!!')

    # Overwrite SetOfCoordinates creation
    def _createSetOfCoordinates(self, suffix=''):
        coordSet = self.inputFullSet.get()
        micSet = coordSet.getMicrographs()
        return ProtSets._createSetOfCoordinates(self, micSet, suffix)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        """Make sure the input data make sense."""

        # Do not allow failing sets:
        notImplentedClasses = ['SetOfClasses2D', 'SetOfClasses3D',
                               'CoordinatesTiltPair']

        errors =[]
        if not self.chooseAtRandom and not self.selectIds and not self.inputSubSet.get():
            errors.append("Subsetting without ids or random selection needs the 'Other set' parameter.")

        if not self.inputFullSet.get():
            # Since is mandatory it will not validate
            # Stop validating since following validations need this set
            return errors

        c1 = self.inputFullSet.get().getClassName()
        if c1 in notImplentedClasses:
            errors.append("%s subset is not implemented." % c1)

        # First dispatch the easy case, where we choose elements at random.
        if self.chooseAtRandom:
            if self.nElements > self.inputFullSet.get().getSize():
                errors.append("Number of elements to choose cannot be bigger than",
                        "the number of elements in the set.")


        # Now the harder case: two sets. Check for compatible classes.

        # self.inputFullSet and self.inputSubSet .get().getClassName()
        # can be SetOf...:
        #   Alignment
        #   Angles
        #   Averages
        #   Classes
        #   ClassesVol
        #   Coordinates
        #   CTF
        #   Micrographs
        #   MovieParticles
        #   Movies
        #   Particles
        #   Volumes

        if not self.inputSubSet.get():
            # Stop validating since following validations need this set
            return errors
        
        c2 = self.inputSubSet.get().getClassName()
        if c2 in notImplentedClasses:
            errors.append("%s subset is not implemented." % c2)

        if c1 == c2:
            return errors

        # Avoid combinations that make no sense.
        for classA, classesIncompatible in [
            ('SetOfParticles',
             {'SetOfMicrographs', 'SetOfMovies', 'SetOfVolumes'}),
            ('SetOfCoordinates',
             {'SetOfMicrographs', 'SetOfMovies', 'SetOfVolumes'}),
            ('SetOfVolumes',
             {'SetOfMicrographs', 'SetOfMovies', 'SetOfParticles', 'SetOfCoordinates'})]:
            if ((c1 == classA and c2 in classesIncompatible) or
                    (c2 == classA and c1 in classesIncompatible)):
                errors.append("The full set and the subset are of incompatible classes",
                        "%s and %s." % (c1, c2))
        return errors

    def _summary(self):
        if self.summaryVar.hasValue():
            return [self.summaryVar.get()]

        key = 'output' + self.inputFullSet.get().getClassName().replace('SetOf', '')

        if not hasattr(self, key):
            return ["Protocol has not finished yet."]
        else:
            if self.setOperation == self.SET_INTERSECTION:
                return ["The elements of %s that also are referenced in %s" %
                        (self.inputFullSet.getName(), self.inputSubSet.getName()),
                        "are now in %s" % getattr(self, key).getName()]
            else:
                return ["%s has elements only present in %s." %
                        (getattr(self, key).getName(),
                         self.inputFullSet.getName())
                        ]


class ProtSubSetByMic(ProtSets):
    """
    Create a subset of those particles that come from a particular set of micrographs
    """
    _label = 'particles subset by micrograph'

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        add = form.addParam  # short notation
        add('inputParticles', pwprot.params.PointerParam,
            pointerClass='SetOfParticles', label="Input particles",
            help='Set of particles from which the subset will be taken')
        add('inputMicrographs', pwprot.params.PointerParam,
            pointerClass='SetOfMicrographs', label="Input micrographs",
            help='Only the particles in this set of micrographs will be output')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep',
                                 self.inputParticles.getObjId(),
                                 self.inputMicrographs.getObjId())

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self, partsId, micsId):
        inputParticles = self.inputParticles.get()
        inputMicrographs = self.inputMicrographs.get()

        outputSet = self._createSetOfParticles()
        outputSet.copyInfo(inputParticles)

        micIds = []

        for mic in inputMicrographs:
            micIds.append(mic.getObjId())

        for particle in inputParticles:
            if particle.getMicId() in micIds:
                outputSet.append(particle)

        self._defineOutputs(outputParticles=outputSet)
        self._defineTransformRelation(inputParticles, outputSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        """Make sure the input data make sense, i.e. hasMicId.
        Thus they come from some Mic"""
        if not self.inputParticles.get().getFirstItem().hasMicId():
            return ['The _Input Particles_ must come from some Micrographs '
                    'of the workflow, i.e. particles must have micId.']

    def _summary(self):
        if not hasattr(self, 'outputParticles'):
            summary = ["Protocol has not finished yet."]
        else:
            summary = ['A subset of *%d* particles is made from a total of *%d*'
                       ' particles.' % (self.outputParticles.getSize(),
                                        self.inputParticles.get().getSize())]
        return summary


class ProtSubSetByCoord(ProtSets):
    """
    Create a subset of those particles that have a particular set of coordinates
    """
    _label = 'particles subset by coordinates'

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        add = form.addParam  # short notation
        add('inputParticles', pwprot.params.PointerParam,
            pointerClass='SetOfParticles', label="Input particles",
            help='Set of particles from which the subset will be taken')
        add('inputCoordinates', pwprot.params.PointerParam,
            pointerClass='SetOfCoordinates', label="Input coordinates",
            help='Only the particles with this set of coordinates will be output')
        add('coordTolerance', pwprot.params.FloatParam,
            label='Coordinate tolerance (px)', default=0,
            help='Two coordinates are supposed to be the same if their X and Y distance'
                 ' is smaller or equal this value')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep',
                                 self.inputParticles.getObjId(),
                                 self.inputCoordinates.getObjId(),
                                 self.coordTolerance.get())

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self, partsId, micsId, tolerance):
        inputParticles = self.inputParticles.get()
        inputCoordinates = self.inputCoordinates.get()

        outputSet = self._createSetOfParticles()
        outputSet.copyInfo(inputParticles)

        micCoordinates = {}
        for coord in inputCoordinates.iterCoordinates():
            micId = coord.getMicId()
            x, y = coord.getPosition()
            if micId not in micCoordinates:
                micCoordinates[micId] = []
            micCoordinates[micId].append((x, y))

        for particle in inputParticles:
            if particle.getMicId() in micCoordinates:
                x0, y0 = particle.getCoordinate().getPosition()
                okToAdd = False
                for x, y in micCoordinates[particle.getMicId()]:
                    if abs(x - x0) <= tolerance and abs(y - y0) <= tolerance:
                        okToAdd = True
                        break
                if okToAdd:
                    outputSet.append(particle)

        self._defineOutputs(outputParticles=outputSet)
        self._defineTransformRelation(inputParticles, outputSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        """Make sure the input data make sense, i.e. hasMicId.
        Thus they come from some Mic"""
        if not self.inputParticles.get().getFirstItem().hasCoordinate():
            return ['The _Input Particles_ must have coordinates']

    def _summary(self):
        if not hasattr(self, 'outputParticles'):
            summary = ["Protocol has not finished yet."]
        else:
            summary = ['A subset of *%d* particles is made from a total of *%d*'
                       ' particles.' % (self.outputParticles.getSize(),
                                        self.inputParticles.get().getSize())]
        return summary

class ProtCrossSubSet(ProtSets):
    """
    Create a subset of the main set based on a matching field in another set. e.g.: Use _micName field (in both fields)
    to select micrographs (main set) present in a set of coordinates (secondary set)
    """
    _label = 'Crossed subset'

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        add = form.addParam  # short notation
        add('mainSet', pwprot.params.PointerParam,
            pointerClass='EMSet', label="Main set",
            help='Set to be reduced')

        add('mainSetField', pwprot.params.StringParam,
            label='Main field', default="id",
            help='Field in the main set that contains the values in common with the secondary set. Use any of the metadata viewers to find the field name.')

        add('secSet', pwprot.params.PointerParam,
            pointerClass='EMSet', label="Secondary set",
            help='Set holding the matching field. e.g: Set of Coordinates hold the micName that can be used to filter a set of micrographs (main set)')

        add('secSetField', pwprot.params.StringParam,
            label='Secondary field', default="id",
            help='Field in the secondary set that contains the values in common with the main set. Use any of the metadata viewers to find the field name.')


    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # These arguments are mainly for skipping the step if they are the same in the resume execution.
        self._insertFunctionStep(self.createOutputStep,
                                 self.mainSet.getObjId(),
                                 self.secSet.getObjId(),
                                 self.mainSetField.get(),
                                 self.secSetField.get())

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self, mainId, secId, mainSetField, secSetField):
        mainSet = self.mainSet.get()
        secSet = self.secSet.get()

        # Instantiate and copy main properties
        outputSet = mainSet.create(self.getPath())
        outputSet.copyInfo(mainSet)

        # Get unique values of secfield in secset
        uniqueValuesinSec = secSet.getUniqueValues(secSetField)
        uniqueValuesinSec ={value:None for value in uniqueValuesinSec}

        mainSetField=self.getMainSetField(pythonName=True)
        isIdField = self._isIdField(mainSetField)
        self.info("Attribute in main set items is %s %s" % (mainSetField, "" if not isIdField else "(id field)"))

        pb = ProgressBar(mainSet.getSize(), fmt=ProgressBar.FULL)
        pb.start()

        for item in mainSet:
            valueInMain=getattr(item,mainSetField)

            valueInMain = valueInMain if isIdField else valueInMain.get()
            if  valueInMain in uniqueValuesinSec:
                self._append(outputSet,item)
            pb.increase()

        pb.finish()

        self._defineOutputs(subset=outputSet)
        self._defineTransformRelation(mainSet, outputSet)

    def getMainSetField(self, pythonName=False):
        if pythonName:
            return self._normalizeSpecialFields(self.mainSetField.get())
        else:
            return self.mainSetField.get()

    def getSecSetField(self, pythonName=False):
        if pythonName:
            return self._normalizeSpecialFields(self.secSetField.get())
        else:
            return self.secSetField.get()

    def _normalizeSpecialFields(self, field):
        if field == ID_COLUMN:
            return ID_ATTRIBUTE
        else:
            return field
    def _isIdField(self, field):
        return field in [ID_COLUMN, ID_ATTRIBUTE]
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        """Make sure the input data make sense"""
        errors=[]
        if not hasattr(self.mainSet.get().getFirstItem(), self.getMainSetField(pythonName=True)):
            errors.append('The main set does not have the field %s' % self.mainSetField.get())

        if not hasattr(self.secSet.get().getFirstItem(), self.getSecSetField(pythonName=True)):
            errors.append('The secondary set does not have the field %s' % self.secSetField.get())

        return errors
    def _summary(self):

        summary = ["Items in the main set where %s=%s of items in the secondary set were selected." % (self.mainSetField.get(), self.secSetField.get())]

        if hasattr(self, "subset"):
            summary.append('*%d* items matched the criteria' % self.subset.getSize())

        return summary


class ProtSetAggregate(EMProtocol):
    """ Aggregates any set data based on its fields"""
    _label = "data summary"
    def _defineParams(self, form):
        form.addSection(label='Input')

        add = form.addParam  # short notation
        add('inputSet', pwprot.params.PointerParam,
            pointerClass='EMSet', label="Any set",
            help='Set with the dta to be aggregated')

        add('operations', pwprot.params.StringParam,
            label='Summary operations',default="COUNT",
            help='Summary operations to apply to all fields in Fields parameter. e.g: MIN MAX AVG. Possible values are MIN, MAX, COUNT, '
                 'AVG, SUM, TOTAL, GROUP_CONCAT. For more technical information see: https://www.sqlite.org/lang_aggfunc.html')

        add('fields', pwprot.params.StringParam,
            label='Fields', default="id",
            help='Fields to apply operations on. Fields can be found in the metadata viewers.'
                  ' The header of the columns are valid names. e.g: _samplingRate id. Fields listed here should '
                 'support the operations specified: DO NOT add literal fields.',
            )

        add('groupby', pwprot.params.StringParam,
            label='Group by',
            help='Fields to make the group. An empty value will summarize the whole dataset.',
            )

    def _insertAllSteps(self):
        self._insertFunctionStep(self.aggregateSet, self.operations.get(), self.fields.get(), self.groupby.get())

    def aggregateSet(self, *args):

        mainSet = self.inputSet.get()

        # Instantiate and copy main properties
        outputSet = SetOfData.create(self.getPath())

        # Run the aggregation method
        operations = self.operations.getListFromValues(caster=str)
        self.info("Operations: %s" % operations)

        fields = self.fields.getListFromValues(caster=str)
        self.info("Fields: %s" % fields)


        if self.groupby.get():
            groupBy = self.groupby.getListFromValues(caster=str)
            self.info("Grouping by: %s" % groupBy)
        else:
            groupBy = None
            self.info("No grouping fields.")


        result = mainSet.aggregate(operations,
                                   fields, groupBy)

        pb = ProgressBar(len(result), fmt=ProgressBar.FULL)
        pb.start()

        # Dictionary to hold the scipion data type based on the key
        scipionTypes ={}

        def getScipionType(fieldName:str):

            if fieldName not in scipionTypes:

                if fieldName.startswith("COUNT"):
                    scipionType=Integer
                elif fieldName.startswith(("MIN","MAX","AVG", "SUM","TOTAL")):
                    scipionType=Float
                else:
                    scipionType=String

                self.info("Scipion type for %s is %s" %(key, scipionType.getClassName()))
                scipionTypes[key] = scipionType
            return scipionTypes[key]

        # Fill the set
        for line in result:
            newItem = Object()
            for key in line.keys():
                scipionType = getScipionType(key)
                value =line[key]
                setattr(newItem, key, scipionType(value))

            outputSet.append(newItem)
            pb.increase()

        pb.finish()

        self._defineOutputs(aggregate=outputSet)
        self._defineTransformRelation(mainSet, outputSet)



