# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *              Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *              Daniel Marchan (da.marchan@cnb.csic.es) - Refactor streaming
# *
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
import enum
from os.path import exists, getmtime
from datetime import datetime

import pyworkflow.object as pwobj
import pyworkflow.protocol.constants as pwcts
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils

import pwem.objects as emobj

from .protocol import EMProtocol


class ProtMicrographs(EMProtocol):
    pass


class ProtCTFMicOutputs(enum.Enum):
    outputCTF = emobj.SetOfCTF


class ProtCTFMicrographs(ProtMicrographs):
    """ Base class for all protocols that estimates the CTF"""

    _possibleOutputs = ProtCTFMicOutputs
    PARALLEL_BATCH_SIZE = 8  # By default in Scipion

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = pwcts.STEPS_PARALLEL
        self.isFirstTime = pwobj.Boolean(False)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=pwutils.Message.LABEL_CTF_ESTI)
        form.addParam('recalculate', params.BooleanParam, default=False,
                      condition='recalculate',
                      label="Do recalculate ctf?")

        form.addParam('continueRun', params.PointerParam, allowsNull=True,
                      condition='recalculate', label="Input previous run",
                      pointerClass=self.getClassName())
        form.addHidden('sqliteFile', params.FileParam, condition='recalculate',
                       allowsNull=True)

        form.addParam('inputMicrographs', params.PointerParam, important=True,
                      condition='not recalculate',
                      label=pwutils.Message.LABEL_INPUT_MIC,
                      pointerClass='SetOfMicrographs')

        form.addParam('AutoDownsampling', params.BooleanParam, default=False,
                      label='Automatic Downsampling Factor',
                      help='Recommended value to downsample')

        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
                      label='Manual CTF Downsampling factor',
                      condition='not AutoDownsampling',  # 'not recalculate',
                      help='Set to 1 for no downsampling. Non-integer downsample '
                           'factors are possible. This downsampling is only used '
                           'for estimating the CTF and it does not affect any '
                           'further calculation. Ideally the estimation of the '
                           'CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) '
                           'and not occupying the whole power spectrum (since '
                           'this downsampling might entail aliasing).')

        self._defineProcessParams(form)

        line = form.addLine('Resolution', condition='not recalculate',
                            help='Give a value in digital frequency '
                                 '(i.e. between 0.0 and 0.5). These cut-offs '
                                 'prevent the typical peak at the center of the'
                                 ' PSD and high-resolution terms where only '
                                 'noise exists, to interfere with CTF '
                                 'estimation. The default lowest value is 0.05 '
                                 'but for micrographs with a very fine sampling '
                                 'this may be lowered towards 0. The default '
                                 'highest value is 0.35, but it should be '
                                 'increased for micrographs with signals '
                                 'extending beyond this value. However, if '
                                 'your micrographs extend further than 0.35, '
                                 'you should consider sampling them at a finer '
                                 'rate.')
        line.addParam('lowRes', params.FloatParam, default=0.05, label='Lowest')
        line.addParam('highRes', params.FloatParam, default=0.35, label='Highest')
        line = form.addLine('Defocus search range (microns)',
                            condition='not recalculate',
                            expertLevel=pwcts.LEVEL_ADVANCED,
                            help='Select _minimum_ and _maximum_ values for '
                                 'defocus search range (in microns). Underfocus'
                                 ' is represented by a positive number.')
        line.addParam('minDefocus', params.FloatParam, default=0.25,
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=4.,
                      label='Max')

        form.addParam('windowSize', params.IntParam, default=512,
                      expertLevel=pwcts.LEVEL_ADVANCED,
                      label='Window size', condition='not recalculate',
                      help='The PSD is estimated from small patches of this '
                           'size. Bigger patches allow identifying more '
                           'details. However, since there are fewer windows, '
                           'estimations are noisier.')

        form.addParallelSection(threads=2, mpi=1)

    def _defineProcessParams(self, form):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific operation."""
        pass

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """ Insert the steps to perform CTF estimation, or re-estimation,
        on a set of micrographs.
        """
        self.initializeStep()
        self._defineCtfParamsDict()
        self.initialIds = self._insertInitialSteps()
        #     self._insertFinalSteps(ctfIds)
        waitCondition = self._getFirstJoinStepName() == 'createOutputStep'
        ctfIds = []

        if self.recalculate:
             if self.isFirstTime:
                 # Insert previous estimation or re-estimation an so on...
                 self._insertPreviousSteps()
                 self.isFirstTime.set(False)

             ctfIds = self._insertRecalculateSteps()
             # For now the streaming is not allowed for recalculate CTF
             waitCondition = False

        self._insertFunctionStep(self.createOutputStep, prerequisites=ctfIds,
                                 wait=waitCondition)

    def initializeStep(self):
        self.insertedDict = {}
        self.micsFn =  self.getInputMicrographs().getFileName()
        # Important to have both:
        self.insertedIds = []  # Contains images that have been inserted in a Step (checkNewInput).
        self.processedIds = [] # Contains images that have been processed in a Step (checkNewOutput).
        self.streamClosed = self.getInputMicrographs().isStreamClosed()
        self.convertCIStep = []
        self.batchSize = self._getStreamingBatchSize()

    def defineBatchSize(self, batchSize, newIds):
        if batchSize == 1:
            batchSize = self.PARALLEL_BATCH_SIZE # Call the program one by one but group images in one step
        else:
            batchSize = len(newIds) # Call the program with all the available images

        return batchSize

    def _insertInitialSteps(self):
        """ Override this function to insert some steps before the
        estimate ctfs steps.
        Should return a list of ids of the initial steps. """

        return []

    def _insertNewMicsSteps(self, newIds, batchSize):
        """ Insert steps to process new mics (from streaming)
        Params:
            inputMics: input mics set to be check
        """
        deps = []
        # Loop through the image IDs in batches
        for i in range(0, len(newIds), batchSize):
            batch_ids = newIds[i:i + batchSize]
            stepId = self._insertFunctionStep(self._insertNewMics, batch_ids,
                                                *self._getCtfArgs(),
                                              prerequisites=self.initialIds)
            for micId in batch_ids:
                self.insertedIds.append(micId)

            deps.append(stepId)

        return deps

    def _insertNewMics(self, micListIds, *args):
        """ Insert steps of new micrographs taking into account a list of micIds
        This function can be used from several base protocols that support
        streaming and batch:
        - ProtCTFMicrographs
        Params:
            micListIds: the input micrographs ids to be inserted into steps
            *args: argument list to be passed to step functions
        Returns:
        """
        micList = []
        inputMicSet = self._loadInputSet(self.micsFn)

        for micId in micListIds:
            mic = inputMicSet.getItem("id", micId).clone()
            micList.append(mic)

        # Now handle the steps depending on the streaming batch size
        batchSize = self.batchSize
        if batchSize == 1:  # This is one by one, as before the batch size
            for mic in micList:
                self.estimateCtfStep(mic, *args)
        else: # Greedy, take all available ones
            self.estimateCtfListStep(micList, *args)

    def _insertRecalculateSteps(self):
        recalDeps = []
        # For each psd insert the steps to process it
        self.recalculateSet = emobj.SetOfCTF(filename=self.sqliteFile.get(),
                                             objDoStore=False)
        inputMics = self.getInputMicrographs()
        for ctf in self.recalculateSet:
            line = ctf.getObjComment()
            if ctf.isEnabled() and line:
                # CTF Re-estimation
                # Make estimation steps independent between them
                objId = ctf.getObjId()
                stepId = self._insertFunctionStep('reEstimateCtfStep',
                                                  objId, prerequisites=[])
                recalDeps.append(stepId)
                self.micDict[objId] = inputMics[objId].clone()
        return recalDeps

    def _insertFinalSteps(self, deps):
        """ This should be implemented in subclasses"""
        return deps

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overwritten in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    # -------------------------- STEPS functions ------------------------------
    def _insertCtfStep(self, mic, prerequisites, *args):
        """ Basic method to insert an estimation step for a given micrograph. """
        micStepId = self._insertFunctionStep('estimateCtfStep',
                                             mic, *args,
                                             prerequisites=prerequisites)
        return micStepId

    def estimateCtfStep(self, mic, *args):
        """ Step function that will be common for all CTF protocols.
        The function _estimateCTF will be called, that should be implemented by each
        CTF estimation protocol. Then, it will take care of registering the micrographs id as processed.
        """
        self.info("Estimating CTF of micrograph: %s " % mic.getObjId())
        self._estimateCTF(mic, *args)
        self.processedIds.append(mic.getObjId())

    def _estimateCTF(self, mic, *args):
        """ Do the CTF estimation with the specific program
        and the parameters required.
        Params:
         micFn: micrograph filename
         micDir: micrograph directory
        """
        raise Exception(pwutils.Message.ERROR_NO_EST_CTF)

    def reEstimateCtfStep(self, micId):
        """ CTF - re-estimation that is common for all programs.
        The _restimateCTF function will be called with proper parameters.
        """
        ctf = self.recalculateSet[micId]
        mic = self.micDict[micId]
        self.info("Estimating CTF of micrograph: %s " % mic.getObjId())
        self._reEstimateCTF(mic, ctf)

    def _reEstimateCTF(self, mic, ctf):
        """ Do the re-estimation of this mic (original one)
        and the parameters that comes in the comment field of the ctf object.
        Params:
         micFn: micrograph filename
         micDir: micrograph directory
        """
        raise Exception(pwutils.Message.ERROR_NO_EST_CTF)

    # Group of functions to estimate several micrographs if the batch size is
    # defined. In some programs it might be more efficient to estimate many
    # at once and not one by one

    def _insertCtfListStep(self, micList, prerequisites, *args):
        """ Basic method to insert an estimation step for a given micrograph. """
        micStepId = self._insertFunctionStep('estimateCtfListStep',
                                             micList, *args,
                                             prerequisites=prerequisites)
        return micStepId

    def estimateCtfListStep(self, micList, *args):
        micIds = [mic.getObjId() for mic in micList]
        self.info("Estimating CTF for micrographs: %s"
                  % micIds)
        self._estimateCtfList(micList, *args)
        self.processedIds.extend(micIds)

    def _estimateCtfList(self, micList, *args):
        """ This function can be implemented by subclasses if it is a more
        efficient way to estimate many micrographs at once.
         Default implementation will just call the _estimateCTF. """
        for mic in micList:
            self._estimateCTF(mic, *args)

    def _createCtfModel(self, mic, updateSampling=False):
        """ This should be implemented in subclasses
        in order to create a CTF model from program results.
        """
        pass

    def createOutputStep(self):
        """ This function is shared by Xmipp and CTFfind
        estimation, or recalculate, protocols.
        if is recalculate, it will iterated for each CTF model, see
        if was recalculated and update with new defocus values.
        Else, the function that should be implemented in each subclass.
        """
        if self.recalculate:
            ctfSet = self._createSetOfCTF("_recalculated")
            prot = self.continueRun.get() or self
            if hasattr(prot, ProtCTFMicOutputs.outputCTF.name):
                micSet = prot.outputCTF.getMicrographs()
                # We suppose this is reading the ctf selection
                # (with enabled/disabled) to only consider the enabled ones
                # in the final SetOfCTF
                # with the recalculate parameters
                newCount = 0
                for ctfModel in self.recalculateSet:
                    if ctfModel.isEnabled() and ctfModel.getObjComment():
                        mic = ctfModel.getMicrograph()
                        # Update the CTF models that where recalculated and append
                        # later to the set, we don't want to copy the id here since
                        # it is already correct
                        newCtf = self._createCtfModel(mic, updateSampling=False)
                        ctfModel.copy(newCtf, copyId=False)
                        ctfModel.setEnabled(True)
                        newCount += 1
                    ctfSet.append(ctfModel)
                ctfSet.setMicrographs(micSet)
                self._defineOutputs(**{ProtCTFMicOutputs.outputCTF.name:ctfSet})
                self._defineCtfRelation(micSet, ctfSet)
                self._computeDefocusRange(ctfSet)
                self.summaryVar.set("CTF Re-estimation of %d micrographs"
                                    % newCount)
            else:
                raise Exception(
                    pwutils.redStr("The outputCTF do not exist, all CTFs failed."))
        else:
            self._createOutputStep()
            if self.outputCTF.getSize() == 0:
                raise Exception(pwutils.redStr("outputCTF has size zero, all CTFs failed."
                                               "Please review processing steps above."))

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []

        if self.recalculate:
            if self.isFinished():
                if self.summaryVar.hasValue():
                    summary.append(self.summaryVar.get())
            else:
                summary.append(pwutils.Message.TEXT_NO_CTF_READY)
        else:
            if not hasattr(self, ProtCTFMicOutputs.outputCTF.name):
                summary.append(pwutils.Message.TEXT_NO_CTF_READY)
            else:
                summary.append("CTF estimation of %d micrographs."
                               % self.inputMicrographs.get().getSize())

        return summary

    def _methods(self):
        methods = []

        if hasattr(self, ProtCTFMicOutputs.outputCTF.name) and self.isFinished():
            methods.append(self.methodsVar.get())
        else:
            methods.append(pwutils.Message.TEXT_NO_CTF_READY)

        return methods

    # -------------------------- UTILS functions ------------------------------
    def _defineCtfParamsDict(self):
        """ This function define a dictionary with parameters used
        for CTF estimation that are common for all micrographs. """
        # Get pointer to input micrographs
        inputMics = self.getInputMicrographs()
        acq = inputMics.getAcquisition()
        sampling = inputMics.getSamplingRate()
        downFactor = self.getAttributeValue('ctfDownFactor', 1.0)
        if downFactor != 1.0:
            sampling *= downFactor
        self._params = {'voltage': acq.getVoltage(),
                        'sphericalAberration': acq.getSphericalAberration(),
                        'magnification': acq.getMagnification(),
                        'ampContrast': acq.getAmplitudeContrast(),
                        'samplingRate': sampling,
                        'scannedPixelSize': inputMics.getScannedPixelSize(),
                        'windowSize': self.windowSize.get(),
                        'lowRes': self.lowRes.get(),
                        'highRes': self.highRes.get(),
                        # Convert from microns to Angstroms
                        'minDefocus': self.minDefocus.get() * 1e+4,
                        'maxDefocus': self.maxDefocus.get() * 1e+4
                        }

    def getCtfParamsDict(self):
        """ Return a copy of the global params dict,
        to avoid overwriting values. """
        return self._params

    def getRecalCtfParamsDict(self, ctfModel):
        """ This function get the acquisition info of the micrographs"""
        mic = ctfModel.getMicrograph()

        acq = mic.getAcquisition()
        mag = acq.getMagnification()
        scannedPixelSize = mic.getSamplingRate() * mag / 10000
        return {'voltage': acq.getVoltage(),
                'sphericalAberration': acq.getSphericalAberration(),
                'magnification': mag,
                'ampContrast': acq.getAmplitudeContrast(),
                'scannedPixelSize': scannedPixelSize,
                'samplingRate': mic.getSamplingRate()
                }

    def _ctfCounter(self, values):
        """ This function return the number of CTFs that was recalculated.
        """
        numberOfCTF = len(values) / 2
        msg = "CTF Re-estimation of %d micrographs" % numberOfCTF
        self.summaryVar.set(msg)

    def _getInputCtf(self):
        if self.continueRecal:
            sqliteFile = self._getPath()
        #             return self.outputCTF.get()
        else:
            return self.inputCtf.get()

    def _iterMicrographs(self, inputMics=None):
        """ Iterate over micrographs and yield
        micrograph name. """
        if inputMics is None:
            inputMics = self.getInputMicrographs()

        for mic in inputMics:
            micFn = mic.getFileName()
            yield micFn, mic

    def _computeDefocusRange(self, ctfSet):
        """ Compute the minimum and maximu defocus in a set of CTFs.
        The protocol methodsVar will be updated with new values.

        Params:
            ctfSet: the set of CTFs to compute min and max
        """
        defocusList = []

        for ctf in ctfSet:
            defocusList.append(ctf.getDefocusU())
            defocusList.append(ctf.getDefocusV())

        minD = min(defocusList) / 10000.
        maxD = max(defocusList) / 10000.

        self.methodsVar.set("Estimated  defocus range defocus was"
                            " %0.3f - %0.3f microns. " % (minD, maxD))

        self._store(self.methodsVar)

    def _defocusMaxMin(self, defocusList):
        """ This function return the minimum and maximum of the defocus
        of a SetOfMicrographs.
        """
        raise Exception("DEPRECATED")

    def getInputMicrographsPointer(self):
        return self.inputMicrographs

    def getInputMicrographs(self):
        return self.getInputMicrographsPointer().get()

    def _getCtfArgs(self):
        """ Should be implemented in sub-classes to define the argument
        list that should be passed to the estimation step function.
        """
        return []

    # ------ Methods for Streaming CTF --------------
    def _stepsCheck(self):
        # To allow streaming ctf estimation we need to detect:
        #   1) new micrographs ready to be estimated
        #   2) new output ctfs that have been produced and add then
        #      to the output set.
        # For now the streaming is not allowed for recalculate CTF
        if self.recalculate:
            return
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # Check if there are new micrographs to process from the input set
        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = datetime.fromtimestamp(getmtime(self.micsFn))
        self.debug('Last check: %s, modification: %s'
                   % (pwutils.prettyTime(self.lastCheck),
                      pwutils.prettyTime(mTime)))

        if self.lastCheck > mTime and self.insertedIds: # If this is empty it is dut to a static "continue" action or it is the first round
            return None

        # Open input movies.sqlite and close it as soon as possible
        micSet = self._loadInputSet(self.micsFn)
        micSetIds = micSet.getIdSet()
        newIds = [idMic for idMic in micSetIds if idMic not in self.insertedIds]

        self.streamClosed = micSet.isStreamClosed()
        self.lastCheck = datetime.now()
        micSet.close()

        outputStep = self._getFirstJoinStep()

        if self.isContinued() and not self.insertedIds:  # For "Continue" action and the first round
            doneIds, size_done_ids = self._getAllDoneIds()
            skipIds = list(set(newIds).intersection(set(doneIds)))
            newIds = list(set(newIds).difference(set(doneIds)))
            self.info("Skipping Mics with ID: %s, seems to be done" % skipIds)
            self.insertedIds = doneIds  # During the first round of "Continue" action it has to be filled

        # Now handle the steps depending on the streaming batch size
        batchSize = self.batchSize
        if batchSize > 1:
            if len(newIds) < batchSize and not self.streamClosed:
                return # No register any step if the batch size is not reach unless is the lass iter
        else:
            batchSize = self.defineBatchSize(batchSize, newIds)

        if newIds:
            fDeps = self._insertNewMicsSteps(newIds, batchSize)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    def _checkNewOutput(self):
        # Load previously done items (from text file)
        doneList, size_done = self._getAllDoneIds()
        processedIds = self.processedIds
        newDoneIds = [micId for micId in processedIds if micId not in doneList]
        allDone = len(doneList) + len(newDoneIds)
        maxMicSize = self._loadInputSet(self.micsFn).getSize()

        # Update the file with the newly done mics
        # or exit from the function if no new done mics
        self.debug('_checkNewOutput: ')
        self.debug('   listOfMics: %s, doneList: %s, newDone: %s'
                   % (maxMicSize, len(doneList), len(newDoneIds)))

        # We have finished when there is not more input mics (stream closed)
        # and the number of processed mics is equal to the number of inputs
        self.finished = self.streamClosed and allDone == maxMicSize
        streamMode = pwobj.Set.STREAM_CLOSED if self.finished else pwobj.Set.STREAM_OPEN

        self.debug('   streamMode: %s newDone: %s' % (streamMode,
                                                      not (newDoneIds == [])))

        if not self.finished and not newDoneIds:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        newDone = []
        inputMicSet = self._loadInputSet(self.micsFn)
        for micId in newDoneIds:
            mic = inputMicSet.getItem("id", micId).clone()
            newDone.append(mic)

        self._updateOutputCTFSet(newDone, streamMode)
        self.debug('   finished: %s ' % self.finished)
        self.debug('        self.streamClosed (%s) AND' % self.streamClosed)
        self.debug('        allDone (%s) == len(self.listOfMics (%s)'
                   % (allDone, maxMicSize))

        if self.finished:  # Unlock createOutputStep if finished all jobs
            self._updateStreamState(streamMode)
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(pwcts.STATUS_NEW)

        self._store()

    def _loadInputSet(self, micsFile):
        """ Load the input set of movies and create a list. """
        self.debug("Loading input db: %s" % micsFile)
        micSet = emobj.SetOfMicrographs(filename=micsFile)
        micSet.loadAllProperties()
        self.streamClosed = micSet.isStreamClosed()
        micSet.close()
        self.debug("Closed db.")
        return micSet

    def _getAllDoneIds(self):
        done_ids = []
        size_output = 0

        if hasattr(self, ProtCTFMicOutputs.outputCTF.name):
            size_output += self.outputCTF.getSize()
            done_ids.extend(list(self.outputCTF.getIdSet()))

        return done_ids, size_output

    def _updateOutputCTFSet(self, micList, streamMode):
        doneFailed = []
        micDoneList = [mic for mic in micList]
        # Do no proceed if there is not micrograph ready
        if not micDoneList:
            return []

        outputName = ProtCTFMicOutputs.outputCTF.name
        outputCtf = getattr(self, outputName, None)

        # If there is not outputCTF yet, it means that is the first
        # time we are updating output CTFs, so we need to first create
        # the output set
        firstTime = outputCtf is None

        if firstTime:
            outputCtf = self._createSetOfCTF()
            outputCtf.setMicrographs(self.getInputMicrographsPointer())
        else:
            outputCtf.enableAppend()

        for micFn, mic in self._iterMicrographs(micList):
            try:
                ctf = self._createCtfModel(mic)
                outputCtf.append(ctf)
            except Exception as ex:
                print(pwutils.yellowStr("Missing CTF?: Couldn't update CTF set with mic: %s" % micFn))
                doneFailed.append(mic)

        self.debug(" _updateOutputCTFSet Stream Mode: %s " % streamMode)
        self._updateOutputSet(outputName, outputCtf, streamMode)
        if doneFailed:
            self._writeFailedList(doneFailed)

        if firstTime:  # define relation just once
            # Using a pointer to define the relations is more robust to
            # scheduling and id changes between the protocol run.db and
            # the main project database.get
            self._defineCtfRelation(self.getInputMicrographsPointer(),
                                    outputCtf)

        return micDoneList

    def _updateStreamState(self, streamMode):
        outputName = ProtCTFMicOutputs.outputCTF.name
        outputCtf = getattr(self, outputName, None)

        # If there are not outputCTFs yet, it means that is the first
        # time we are updating output CTF, so we need to first create
        # the output set
        firstTime = outputCtf is None

        if firstTime:
            outputCtf = self._createSetOfCTF()
        else:
            outputCtf.enableAppend()

        self.debug(" _updateStreamState Stream Mode: %s " % streamMode)
        self._updateOutputSet(outputName, outputCtf, streamMode)

    def _getAllFailed(self):
        return self._getExtraPath('FAILED_all.TXT')

    def _getMicrographDir(self, mic):
        """ Return an unique dir name for results of the micrograph. """
        return self._getTmpPath('mic_%04d' % mic.getObjId())

    def _writeFailedList(self, micList):
        """ Write to a text file the items that have failed. """
        with open(self._getAllFailed(), 'a') as f:
            for mic in micList:
                f.write('%d\n' % mic.getObjId())

    def _readFailedList(self):
        """ Read from a text file the id's of the items that have failed. """
        failedFile = self._getAllFailed()
        failedList = []
        if exists(failedFile):
            with open(failedFile) as f:
                failedList += [int(line.strip()) for line in f]

        return failedList


class ProtPreprocessMicrographs(ProtMicrographs):
    pass
