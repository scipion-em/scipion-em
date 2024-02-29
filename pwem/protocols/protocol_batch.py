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

import os
try:
    from itertools import izip
except ImportError:
    izip = zip

from json import dumps, loads

from pyworkflow.protocol.params import (PointerParam, FileParam, StringParam)
from pyworkflow.utils import moveFile

from pwem.protocols import EMProtocol
import pwem.objects as emobj
from pwem.constants import ALIGN_NONE
from pwem.objects import (MicrographsTiltPair, ParticlesTiltPair)


class BatchProtocol(EMProtocol):
    """ Base class to all protocols that are launched
    through other GUIs (such as showj) and that
    are called "batch" protocols. They should not be
    executed from normal "form" of other protocols.
    """
    pass


class ProtUserSubSet(BatchProtocol):
    """ Create subsets from the GUI.
    This protocol will be executed mainly from the ShowJ gui or the new metadata viewer.
    For ShowJ: The enabled/disabled changes will be stored in a temporary sqlite
    For the Metadata viewer: the original sqlite is passed and a path of a txt file with the ids selected.

    Example of a ShowJ socket message:

     run protocol ProtUserSubSet inputObject=380 sqliteFile='Runs/000335_ProtImodTSNormalization/tiltseries_state.sqlite','' outputClassName=SetOfTiltSeries other='' label='create subset'

     run protocol is received by ProjectTCPRequestHandler.handle()

    """
    def __init__(self, **args):
        BatchProtocol.__init__(self, **args)
        self._selectedIds = None
        self._dbName = None # To keep the path to the database to read
        self._selectionTxt = None # To keep the file (txt) that has the selection. Extracted from sqliteFile tuple.

    def _defineParams(self, form):
        form.addSection("Input")
        form.addParam('inputObject', PointerParam, label="Input set", pointerClass='EMSet')
        form.addParam('other', StringParam, allowsNull=True, label="Other")
        form.addParam('sqliteFile', FileParam, label="Selection file")
        form.addParam('outputClassName', StringParam, label="Output type")

    def _insertAllSteps(self):
        self._insertFunctionStep('createSetStep')

    def createSetStep(self):

        sourceSet = self.inputObject.get()
        markedSet = self.createSetObject()  # Set equal to sourceSet but marked with disabled
        other = self.other.get()

        self.info("Source: %s" % sourceSet) # Scipion object input for this protocol
        self.info("Output type: %s" % self.outputClassName)
        self.info("Subset file: %s" % self.sqliteFile) # Sqlite with the selection (ShowJ case) or
        if other:
            self.info("Other: %s" % other)


        # New recommended way to create subsets: making the set responsible for his own subset process
        # Once all Sets implement appendFromSet and the if below is gone we can remove this "if"
        if getattr(markedSet, "USE_CREATE_COPY_FOR_SUBSET", False):
            markedSet.loadAllProperties()
            newSet = markedSet.createCopy(self._getPath(),
                                          copyInfo=True, copyItems=True)

            # Define outputs, may be use something more specific than "subset"
            self._defineOutputs(subset=newSet)
            self._defineSourceRelation(sourceSet, newSet)
            return

        if other and ',Volume' in other:
            volId = int(other.split(',')[0])

            if isinstance(markedSet, emobj.SetOfVolumes):
                volSet = emobj.SetOfVolumes(filename=self._dbName)
                output = volSet[volId]
            else:
                classSet = emobj.SetOfClasses3D(filename=self._dbName)
                output = classSet[volId].getRepresentative()
            self._defineOutputs(outputVolume=output)

        elif isinstance(sourceSet, emobj.SetOfImages):
            self._createSubSetFromImages(sourceSet)

        elif isinstance(sourceSet, emobj.SetOfClasses):
            self._createSubSetFromClasses(sourceSet)

        elif isinstance(sourceSet, emobj.SetOfCTF):
            outputClassName = self.outputClassName.get()
            if outputClassName.startswith('SetOfMicrographs'):
                self._createMicsSubSetFromCTF(sourceSet)

        elif isinstance(sourceSet, emobj.SetOfAtomStructs):
            self._createSubSetFromAtomStructs(sourceSet)

        elif isinstance(sourceSet, MicrographsTiltPair):
            self._createSubSetFromMicrographsTiltPair(sourceSet)

        elif isinstance(sourceSet, ParticlesTiltPair):
            self._createSubSetFromParticlesTiltPair(sourceSet)

        elif isinstance(sourceSet, EMProtocol):
            if self.other.hasValue():
                otherid = self.other.get()
                otherObj = self.getProject().mapper.selectById(int(otherid))

                if isinstance(markedSet, emobj.SetOfClasses):
                    markedSet.setImages(otherObj)
                    self._createSubSetFromClasses(markedSet)

                elif isinstance(markedSet, emobj.SetOfImages):
                    markedSet.copyInfo(otherObj)  # copy info from original images
                    self._createSubSetFromImages(markedSet)

                elif isinstance(markedSet, emobj.SetOfNormalModes):
                    self._createSimpleSubset(otherObj)
            else:
                if isinstance(markedSet, emobj.SetOfVolumes):
                    volSet = emobj.SetOfVolumes(filename=self._dbName)
                    volSet.loadAllProperties()
                    self._createSimpleSubset(volSet)

                # Go for a generic way of creating the set the the
                # input set is not registered (typically from viewers)
                else:
                    # We might want to do this before, inside the createSetObject
                    markedSet.loadAllProperties()
                    self._createSimpleSubset(markedSet)

        else:
            self._createSimpleSubset(sourceSet)

    def _createSimpleSubset(self, inputObj):
        className = inputObj.getClassName()
        modifiedSet = inputObj.getClass()(filename=self._dbName,
                                          prefix=self._dbPrefix)
        try:
            createFunc = getattr(self, '_create' + className)
            output = createFunc()
        except Exception:
            output = inputObj.createCopy(self._getPath())

        for item in modifiedSet:
            if self._itemSelected(item):
                output.append(item)

        if hasattr(modifiedSet, 'copyInfo'):
            output.copyInfo(inputObj)

        # Register outputs
        self._defineOutput(className, output)

        if inputObj.hasObjId():
            self._defineTransformRelation(inputObj, output)

        return output

    def usingShowJ(self):
        return ".sqlite" in self.sqliteFile.get()

    def _itemSelected(self, item):
        """ Returns true if the element is selected. We will have to modes: ShowJ or Metadata viewer
        ShowJ marks the items as enable=False
        Metadata viewer, passes a file with the id selected

        """
        if self.usingShowJ():
            return item.isEnabled()
        else:
            return self._itemInSelectionTxt(item)

    def _itemInSelectionTxt(self, item):

        id = item.getObjId()
        return id in self._getSelectionTxtDict()

    def _createSubSetFromImages(self, inputImages,
                                copyInfoCallback=None):
        className = inputImages.getClassName()
        setClass = inputImages.getClass()

        modifiedSet = setClass(filename=self._dbName, prefix=self._dbPrefix)
        try:
            createFunc = getattr(self, '_create' + className)
            output = createFunc()
        except Exception:
            output = inputImages.createCopy(self._getPath())

        if copyInfoCallback is None:
            # modifiedSet.loadAllProperties()
            output.copyInfo(inputImages)
        else:
            copyInfoCallback(output)

        output.appendFromImages(modifiedSet, self._itemSelected)
        # Register outputs
        self._defineOutput(className, output)

        if inputImages.hasObjId():
            self._defineTransformRelation(inputImages, output)

        # Define an informative summary of the subset operation
        sizeIn = inputImages.getSize()
        sizeOut = output.getSize()
        sizeDiff = sizeIn - sizeOut
        msg = 'A subset of _%s_ was created, ' % output.getClassName()
        msg += ('discarding *%d* items (%0.1f %%) from the input set.' %
                (sizeDiff, sizeDiff*100./sizeIn))
        self.summaryVar.set(msg)

        return output

    def _getSelectionTxtDict(self):

        if self._selectedIds is None:
            self._selectedIds = dict()
            ids = None
            self.info("Reading selection from %s" % self._dbName)
            # This file has a single line with id separated by spaces
            with open(self._selectionTxt, "r") as fh:
                line = fh.readline().strip()
                ids = line.split(" ")

            for id in ids:
                self._selectedIds[int(id)] = id

        return self._selectedIds

    def _createSubSetFromClasses(self, inputClasses):
        outputClassName = self.outputClassName.get()

        if (outputClassName.startswith('SetOfAverages') or
            outputClassName.startswith('SetOfVolumes') or
                outputClassName.startswith('SetOfParticles')):
            # We need to distinguish two cases:
            # a) when we want to create images by grouping class images
            # b) create a subset from a particular class images
            from pyworkflow.mapper.sqlite import SqliteFlatDb
            db = SqliteFlatDb(dbName=self._dbName, tablePrefix=self._dbPrefix)
            itemClassName = db.getSelfClassName()

            if itemClassName.startswith('Class'):

                if outputClassName.startswith('SetOfParticles'):
                    # Just to be sure, we check first if we want to create a subset of the items of the class
                    # If this check is not done in first place, it might be possible to enter through the second
                    # check (if REP_SET_TYPE is set), generating a SetOfVolumes like instead of a SetOfParticles like
                    # object
                    return self._createImagesFromClasses(inputClasses)

                elif hasattr(inputClasses, "REP_SET_TYPE"):
                    # Second check determines if we want to create a subset of representatives. In addition, we consider
                    # that the SetOfClasses stores the information about the type of set to be created
                    return self._createRepresentativesFromClasses(inputClasses,
                                                                  getattr(inputClasses, "REP_SET_TYPE"))
                else:
                    # Third check determines if we want to create a subset of representatives. By default,
                    # a SetOfVolumes is generated
                    return self._createRepresentativesFromClasses(inputClasses,
                                                                  outputClassName.split(',')[0])
            else:
                def callback(output):
                    self._copyInfoAndSetAlignment(inputClasses, output)

                return self._createSubSetFromImages(inputClasses.getImages(),
                                                    copyInfoCallback=callback)

        elif outputClassName.startswith('SetOfClasses'):
            return self._createClassesFromClasses(inputClasses)
        else:
            raise Exception("Unrecognized output type: '%s'" % outputClassName)

    def _createMicsSubSetFromCTF(self, inputCTFs):
        """ Create a subset of Micrographs and CTFs when analyzing the CTFs. """
        outputMics = self._createSetOfMicrographs()
        outputCtfs = self._createSetOfCTF()
        setOfMics = inputCTFs.getMicrographs()
        if setOfMics is None:
            raise Exception('Could not create SetOfMicrographs subset from '
                            'this SetOfCTF, the micrographs were not set.')
        outputMics.copyInfo(setOfMics)

        modifiedSet = emobj.SetOfCTF(filename=self._dbName,
                                     prefix=self._dbPrefix)

        count = 0
        for ctf in modifiedSet:
            if self._itemSelected(ctf):
                mic = ctf.getMicrograph()
                outputMics.append(mic)
                outputCtfs.append(ctf)
                count += 1

        # Register outputs
        outputCtfs.setMicrographs(outputMics)
        # NOTE: I've split the define output in 2 steps.
        # It seems with python3 outputCTF was processed first and needs mics to be saved first.
        self._defineOutputs(outputMicrographs=outputMics)
        self._defineOutputs(outputCTF=outputCtfs)
        self._defineTransformRelation(setOfMics, outputMics)
        self._defineCtfRelation(outputMics, outputCtfs)
        msg = 'From input %s of size %s created output ' % (inputCTFs.getClassName(),
                                                            inputCTFs.getSize())
        msg += 'SetOfMicrographs and SetOfCTF of size %d' % count
        self.summaryVar.set(msg)

        return outputMics, outputCtfs

    def _createRepresentativesFromClasses(self, inputClasses, outputClassName):
        """ Create a new set of images joining all images
        assigned to each class.
        """
        inputImages = inputClasses.getImages()

        if isinstance(outputClassName, str):
            createFunc = getattr(self, '_create' + outputClassName)
            output = createFunc()
        else:
            output = outputClassName.create(self.getPath())

        modifiedSet = inputClasses.getClass()(filename=self._dbName, prefix=self._dbPrefix)
        self.info("Creating REPRESENTATIVES of images from classes, "
                  "sqlite file: %s" % self._dbName)

        count = 0

        output.copyInfo(inputImages)
        output.setSamplingRate(None)
        # For now this is to avoid having a wrong alignment.
        # THis is because is getting the alignment info from the input images and this does not have to match.
        # This created an error when scaling averages #903
        output.setAlignment(ALIGN_NONE)
        for cls in modifiedSet:
            if self._itemSelected(cls):
                img = cls.getRepresentative()
                if not output.getSamplingRate():
                    output.setSamplingRate(cls.getSamplingRate()
                                           if img.getSamplingRate() is None
                                           else img.getSamplingRate())
                img.copyObjId(cls)
                output.append(img)
                count += 1
        # Register outputs
        self._defineOutput('Representatives', output)
        if inputClasses.hasObjId():
            self._defineSourceRelation(inputClasses, output)
        else:
            self._defineSourceRelation(inputImages, output)

        selectmsg = 'we selected %s items' % count if count > 1 else 'was selected 1 item'
        msg = 'From input %s of size %s %s to create output %s' % (inputClasses.getClassName(),
                                                                   inputClasses.getSize(),
                                                                   selectmsg,
                                                                   output.getClassName())
        self.summaryVar.set(msg)
        return output

    @staticmethod
    def _copyInfoAndSetAlignment(inputClasses, output):
        """ This method is used when creating subset of images from classes.
        We need to copy the information from the original classes images
        and also set the proper alignment contained in the classes.
        """
        inputImages = inputClasses.getImages()
        # Copy all info form the original 'classified' images
        output.copyInfo(inputImages)
        # Take the alignment of the first class
        cls = inputClasses.getFirstItem()
        output.setAlignment(cls.getAlignment())

    def _createImagesFromClasses(self, inputClasses):
        """ Create a new set of images joining all images
        assigned to each class.
        """
        inputImages = inputClasses.getImages()
        className = inputImages.getClassName()
        try:
            createFunc = getattr(self, '_create' + className)
            output = createFunc()
        except Exception as e:
            output = inputImages.createCopy(self._getPath())

        modifiedSet = inputClasses.getClass()(filename=self._dbName, prefix=self._dbPrefix)
        self.info("Creating subset of images from classes, sqlite file: %s" % self._dbName)
        self._copyInfoAndSetAlignment(inputClasses, output)
        output.appendFromClasses(modifiedSet, filterClassFunc=self._itemSelected)
        # Register outputs
        self._defineOutput(className, output)
        if inputClasses.hasObjId():
            self._defineSourceRelation(inputClasses, output)
        self._defineTransformRelation(inputImages, output)
        count = len([cls for cls in modifiedSet if cls.isEnabled()])
        selectmsg = 'we selected %s items' % count if count > 1 else 'was selected 1 item'
        msg = 'From input %s of size %s %s to create output %s of size %s' % (inputClasses.getClassName(),
                                                                              inputClasses.getSize(),
                                                                              selectmsg,
                                                                              output.getClassName(),
                                                                              output.getSize())
        self.summaryVar.set(msg)
        return output

    def _createClassesFromClasses(self, inputClasses):
        """ Create a new set of images joining all images
        assigned to each class.
        """
        # inputImages = inputClasses.getImages()
        className = inputClasses.getClassName()
        # createFunc = getattr(self, '_create' + className)

        try:
            createFunc = getattr(self, '_create' + className)
            output = createFunc()
        except Exception as e:
            output = inputClasses.getClass().create(self._getPath())

        modifiedSet = inputClasses.getClass()(filename=self._dbName, prefix=self._dbPrefix)
        self.info("Creating subset of classes from classes, sqlite file: %s" % self._dbName)
        output.copyInfo(inputClasses)
        output.appendFromClasses(modifiedSet,filterClassFunc=self._itemSelected)
        # Register outputs
        self._defineOutput(className, output)
        if inputClasses.hasObjId():
            self._defineTransformRelation(inputClasses, output)
        else:
            self._defineSourceRelation(inputClasses.getImages(), output)
        count = len([cls for cls in modifiedSet if cls.isEnabled()])
        selectmsg = 'we selected %s items' % count if count > 1 else 'was selected 1 item'
        msg = 'From input %s of size %s %s to create output %s' % (inputClasses.getClassName(),
                                                                   inputClasses.getSize(),
                                                                   selectmsg,
                                                                   output.getClassName())
        self.summaryVar.set(msg)
        return output

    def _createSubSetFromMicrographsTiltPair(self, micrographsTiltPair):
        """ Create a subset of Micrographs Tilt Pair. """
        output = MicrographsTiltPair(filename=self._getPath('micrographs_pairs.sqlite'))
        modifiedSet = MicrographsTiltPair(filename=self._dbName,
                                          prefix=self._dbPrefix)
        inputU = micrographsTiltPair.getUntilted()
        inputT = micrographsTiltPair.getTilted()
        outputU = emobj.SetOfMicrographs(filename=self._getPath('mics_untilted.sqlite'))
        outputT = emobj.SetOfMicrographs(filename=self._getPath('mics_tilted.sqlite'))
        outputU.copyInfo(inputU)
        outputT.copyInfo(inputT)

        # noinspection DuplicatedCode
        for micPair, u, t in izip(modifiedSet, inputU, inputT):
            if micPair.isEnabled():
                output.append(micPair)
                outputU.append(u)
                outputT.append(t)
        output.setUntilted(outputU)
        output.setTilted(outputT)
        # Register outputs
        outputDict = {'outputMicrographsTiltPair': output}
        self._defineOutputs(**outputDict)
        self._defineTransformRelation(micrographsTiltPair, output)
        return output

    def _createSubSetFromAtomStructs(self, setOfPDBs):
        """ Create a subset of SetOfAtomStruct. """
        output = emobj.SetOfAtomStructs(filename=self._getPath('atomstructs.sqlite'))
        modifiedSet = emobj.SetOfAtomStructs(filename=self._dbName, prefix=self._dbPrefix)

        for pdb in modifiedSet:
            if pdb.isEnabled():
                output.append(pdb)

        # Register outputs
        outputDict = {'outputAtomStructs': output}
        self._defineOutputs(**outputDict)
        self._defineTransformRelation(setOfPDBs, output)
        return output

    def _createSubSetFromParticlesTiltPair(self, particlesTiltPair):
        """ Create a subset of Particles Tilt Pair. """
        output = ParticlesTiltPair(filename=self._getPath('particles_pairs.sqlite'))

        inputU = particlesTiltPair.getUntilted()
        inputT = particlesTiltPair.getTilted()
        outputU = emobj.SetOfParticles(filename=self._getPath('particles_untilted.sqlite'))
        outputT = emobj.SetOfParticles(filename=self._getPath('particles_tilted.sqlite'))
        outputU.copyInfo(inputU)
        outputT.copyInfo(inputT)

        modifiedSet = ParticlesTiltPair(filename=self._dbName,
                                        prefix=self._dbPrefix)

        # noinspection DuplicatedCode
        for pair, u, t in izip(modifiedSet, inputU, inputT):
            if pair.isEnabled():
                output.append(pair)
                outputU.append(u)
                outputT.append(t)
        # Register outputs
        output.setUntilted(outputU)
        output.setTilted(outputT)
        # Link output to the same coordinates pairs than input
        output.setCoordsPair(particlesTiltPair.getCoordsPair())

        outputDict = {'outputParticlesTiltPair': output}
        self._defineOutputs(**outputDict)
        self._defineTransformRelation(particlesTiltPair, output)
        return output

    # noinspection PyAttributeOutsideInit
    def createSetObject(self):
        """ Moves the sqlite with the enable/disable status to its own
        path to keep it and names it subset.sqlite"""

        _dbName, self._dbPrefix = self.sqliteFile.get().split(',')

        if self.usingShowJ():

            self._dbName = self._getPath('subset.sqlite')
            os.rename(_dbName, self._dbName)
        else:
            self._dbName = self.inputObject.get().getFileName()
            self._selectionTxt = _dbName

        # Prefix: used to create subsets from Particles from a class (specific table in a set of classes)
        if self._dbPrefix.endswith('_'):
            self._dbPrefix = self._dbPrefix[:-1]

        from pwem.utils import loadSetFromDb

        # Ignoring self._dbPrefix here, since we want to load
        # the top-level set in the sqlite file
        setObj = loadSetFromDb(self._dbName)
        return setObj

    def _summary(self):
        summary = []
        msg = self.summaryVar.get()
        if msg is None:
            msg = self.getDefaultSummary()
        summary.append(msg)
        return summary

    def getDefaultSummary(self):
        inputStr = ''
        inputObj = self.inputObject.get()
        if inputObj is not None:
            inputStr += inputObj.getClassName()
            if isinstance(inputObj, emobj.EMSet):
                inputStr += ' of size %s' % inputObj.getSize()
        output = ''
        for _, attr in self.iterOutputAttributes():
            output += attr.getClassName()
            if isinstance(attr, emobj.EMSet):
                output += ' of size %s' % attr.getSize()

        msg = 'From input %s created output %s ' % (inputStr, output)

        return msg

    def _methods(self):
        return self._summary()

    def _defineOutput(self, className, output):
        outputDict = {'output' + className.replace('SetOf', ''): output}
        self._defineOutputs(**outputDict)


class ProtCreateMask(BatchProtocol):

    def _defineParams(self, form):
        form.addHidden('inputObj', PointerParam, pointerClass='EMObject')
        form.addHidden('maskFile', StringParam)

    def _insertAllSteps(self):
        self._insertFunctionStep('createMaskStep')

    def createMaskStep(self):
        inputObj = self.inputObj.get()
        maskSrc = self.maskFile.get()
        basename = os.path.basename(maskSrc)
        maskDst = self._getPath(basename)
        moveFile(maskSrc, maskDst)
        samplingRate = None
        if hasattr(inputObj, "getSamplingRate"):
            samplingRate = inputObj.getSamplingRate()
        else:
            for key, attr in inputObj.iterInputAttributes():
                if hasattr(attr.get(), "getSamplingRate"):
                    samplingRate = attr.get().getSamplingRate()
        if not samplingRate:
            raise Exception("sampling rate required")

        mask = emobj.Mask()
        mask.setFileName(maskDst)
        mask.setSamplingRate(samplingRate)
        self._defineOutputs(outputMask=mask)
        self._defineSourceRelation(self.inputObj, self.outputMask)

    def _summary(self):
        summary = list()
        summary.append('From input %s created mask %s'
                       % (self.getObjectTag("inputObj"),
                          self.getObjectTag("outputMask")))
        return summary

    def _methods(self):
        return self._summary()


class ProtCreateFSC(BatchProtocol):

    def _defineParams(self, form):
        form.addHidden('inputObj', PointerParam,
                       pointerClass='EMObject')
        form.addHidden('fscValues', StringParam,
                       help='String representation of the list with FSC values')
        form.addHidden('fscLabels', StringParam,
                       help='String with fsc labels')

    def setInputObj(self, obj):
        self.inputObj.set(obj)

    def setInputFscList(self, fscList):
        fscStr = ''
        fscLabel = ""
        numberFsc = 0
        for fsc in fscList:
            if numberFsc == 0:
                numberFsc += 1
            else:
                fscStr += "|"
                fscLabel += "|"
            fscStr += dumps(fsc.getData())
            fscLabel += dumps(fsc.getObjLabel())
        self.fscValues.set(fscStr)
        self.fscLabels.set(fscLabel)

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        fscSet = self._createSetOfFSCs()
        fscSet.setObjLabel("setOfFSCs")
        dataStringList = self.fscValues.get().split("|")
        labelStringList = self.fscLabels.get().split("|")
        for fsc, label in zip(dataStringList, labelStringList):
            _fsc = emobj.FSC(objLabel=loads(label))
            freq, value = loads(fsc)
            _fsc.setData(freq, value)
            fscSet.append(_fsc)
        self._defineOutputs(outputFSCs=fscSet)

    def _summary(self):
        summary = list()
        summary.append('From input %s created FSC %s'
                       % (self.getObjectTag("inputObj"),
                          self.getObjectTag("outputFSCs")))
        return summary

    def _methods(self):
        return self._summary()

