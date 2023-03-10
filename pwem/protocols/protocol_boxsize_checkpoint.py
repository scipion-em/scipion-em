# **************************************************************************
# *
# * Authors:     Daniel Marchán (da.marchan@cnb.csic.es)
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
from pyworkflow.object import Integer
from pyworkflow.protocol.params import  StringParam
import numpy as np
import re
import time
import os

BOXSIZE= 'boxSize'

class ProtBoxSizeCheckpoint(EMProtocol):
    """
    Protocol to make a validation operations on particle picking boxsize.
    For sanity check all the generated outputs are even numbers.
    """

    _label = 'box size checkpoint'
    _possibleOutputs = {BOXSIZE: Integer}
    outputsToDefine = {}


    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _defineParams(self, form):
        """
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """
        form.addSection(label='Input')

        form.addParam('boxSize1', params.IntParam,
                      label='Particle box size reference (px)',
                      allowsPointers=True,
                      important=True,
                      help='This is the reference size of the boxed particles (in pixels).\n'
                           'For sanity check if it is not, it will be transform to an even number.')

        form.addParam('boxSize2', params.IntParam,
                      label='Particle box size secondary reference (px)',
                      allowsPointers=True,
                      important=True,
                      help='This is a secondary size of the boxed particles (in pixels).\n'
                           'For sanity check if it is not, it will be transform to an even number.')

        form.addParam('disagreeFactor', params.FloatParam, default=0.2,
                         label='Maximum proportional difference',
                         help='This proportion is calulated with the following formula:\n'
                              'Proportional diff = 1 - min(boxSize1, boxSize2)/max(boxSize1, boxSize2)')

        form.addParam('boolBoxSize1', params.BooleanParam, default=True,
                      label='If disagreed stayed with the primary reference?',
                      help='Select yes if you want to use the primary reference as output.')

        form.addParam('boolBoxSizeAvg', params.BooleanParam, default=False,
                      label='Average the particle boxsize?',
                      help='Select yes if you want to use an average of the box sizes as output.')

        form.addParam('boolTimer', params.BooleanParam, default=False,
                      label='Use timer?',
                      help='Select yes if you want to use a timer to take the decision.')

        form.addParam('timeout', StringParam, default="1h",
                      condition="boolTimer==%d" % True,
                      label='Time to wait:',
                      help='Time in seconds that the protocol will remain '
                           'running. A correct format is an integer number in '
                           'seconds or the following syntax: {days}d {hours}h '
                           '{minutes}m {seconds}s separated by spaces '
                           'e.g: 1d 2h 20m 15s,  10m 3s, 1h, 20s or 25.')


    def _insertAllSteps(self):
        self.initParams()
        self._checkNewInput()

    def initParams(self):
        self.outputDone = False

    def createOutput(self, modifiedSet):
        pass

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        if self.outputDone:
            return None

        boxSize1 = self.boxSize1.get()
        boxSize2 = self.boxSize2.get()
        fDeps = self._insertNewOperationsStep(boxSize1, boxSize2)
        self.updateSteps()

    def _insertNewOperationsStep(self, boxSize1, boxSize2):
        deps = []
        stepId = self._insertFunctionStep('applyComparisonStep', boxSize1, boxSize2, prerequisites=[])
        deps.append(stepId)
        return deps

    def _checkNewOutput(self):
        if self.outputDone:
            self.createResultsOutput()

    def applyComparisonStep(self, boxSize1, boxSize2):
        """
        Applies the formula to make the comparison of the two box size. It has four cases:
            -  Stays within the limits of the maximum proportional difference and outputs the average box size (default).
            -  Stays within the limits of the maximum proportional difference and outputs the primary reference.
            -  Goes out of the maximum proportional difference and outputs the primary reference (default).
            -  Goes out of the maximum proportional difference and outputs the average box size.
        """

        diff_p = 1 - min(boxSize1, boxSize2)/max(boxSize1, boxSize2)
        self.info("Calculating the difference between reference box size %d and secondary box size %d"
                  %(boxSize1, boxSize2))

        if diff_p <= self.disagreeFactor.get():
            self.info("The proportional difference was within the threshold %f<%f"
                      %(diff_p, self.disagreeFactor.get()))
            if self.boolBoxSizeAvg.get():
                boxSize = np.mean([boxSize1, boxSize2])
                self.info("Using average box size %d" %boxSize)
            else:
                boxSize = boxSize1
                self.info("Using reference box size %d" %boxSize)
        else:
            self.info("The proportional difference was out of threshold %f>%f"
                      % (diff_p, self.disagreeFactor.get()))
            if self.boolBoxSize1.get():
                boxSize = boxSize1
                self.info("Using reference box size %d" % boxSize)
            else:
                boxSize = np.mean([boxSize1, boxSize2])
                self.info("Using average box size %d" % boxSize)

        if self.boolTimer.get():
            self._insertFunctionStep(self.waitingStep)

        self.registerEvenBoxSize(boxSize)
        self.outputDone = True

    def waitingStep(self):
        lastTimeCheck = 0
        timeout = self.getTimeOutInSeconds(self.timeout.get())
        finished = False
        while not finished:
            if lastTimeCheck % 300 == 0:
                self.info("Remaining time: %s seconds. Next check in 5 minutes" % str(timeout - lastTimeCheck))
            time.sleep(10)
            lastTimeCheck += 10
            finished = (lastTimeCheck > timeout or
                        self.waitingHasFinished())

    def getTimeOutInSeconds(self, timeOut):
        timeOutFormatRegexList = {'\d+s':1, '\d+m':60, '\d+h':3600,
                                  '\d+d':72000}
        try:
            return int(timeOut)
        except Exception:
            seconds = 0
        for regex, secondsUnit in timeOutFormatRegexList.items():
            matchingTimes = re.findall(regex, timeOut)
            for matchTime in matchingTimes:
                seconds += int(matchTime[:-1]) * secondsUnit

        return seconds

    def _getStopWaitingFilename(self):
        return self._getExtraPath("APPROVE_CHECKPOINT.TXT")

    def getActions(self):
        if self.isRunning():
            return [('APPROVE CHECKPOINT', self.stopWait)]
        else:
            return []

    def stopWait(self):
        f = open(self._getStopWaitingFilename(), 'w')
        f.close()

    def waitingHasFinished(self):
        return os.path.exists(self._getStopWaitingFilename())

    def registerOutput(self, outputName, value):
        self.outputsToDefine[outputName] = value

    def registerEvenBoxSize(self, boxSize):
        boxSize = transform2EvenNumber(boxSize)
        self.registerOutput(BOXSIZE, Integer(boxSize))

    def createResultsOutput(self):
        """ The output is an even Integer. Other protocols can use it in those
            Params if it has set allowsPointer=True
        """
        self._defineOutputs(**self.outputsToDefine)

    def _summary(self):
        summary = []

        return summary

    def _validate(self):
        errors = []

        if not self.boolBoxSize1.get() and not self.boolBoxSizeAvg.get():
            errors.append('One of the two options should be selected as True.')

        if self.boolTimer.get():
            timeInSec = self.getTimeOutInSeconds(self.timeout.get())
            if timeInSec == 0:
                message.append('Time format is wrong. A correct format is an '
                               'integer number in seconds or the following syntax: '
                               '{days}d {hours}h {minutes}m {seconds}s separated by'
                               ' spaces e.g: 1d 2h 20m 15s,  10m 3s, 1h, 20s or 25.')

        return errors


def transform2EvenNumber(var):
    if var % 2 != 0:
        var = round(var / 2) * 2

    return var