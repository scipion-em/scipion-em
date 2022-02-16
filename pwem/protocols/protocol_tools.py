# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
import re
import time

from pyworkflow.protocol.params import  StringParam

from pwem.protocols import EMProtocol


class ProtManualCheckpoint(EMProtocol):
    """
    This protocol is kept running for a time determined by a parameter or until
    the user determines it is convenient. This protocol is useful if we want a
    given protocol to be launched at the time the user sees appropriate.
    """
    _label = 'manual check point'

    def _defineParams(self, form):
        form.addSection(label='Criteria')
        form.addParam('timeout', StringParam, default="1h",
                      label='Time to wait:',
                      help='Time in seconds that the protocol will remain '
                           'running. A correct format is an integer number in '
                           'seconds or the following syntax: {days}d {hours}h '
                           '{minutes}m {seconds}s separated by spaces '
                           'e.g: 1d 2h 20m 15s,  10m 3s, 1h, 20s or 25.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.waitingStep)

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

    def _validate(self):
        message = []
        timeInSec = self.getTimeOutInSeconds(self.timeout.get())
        if timeInSec == 0:
            message.append('Time format is wrong. A correct format is an '
                           'integer number in seconds or the following syntax: '
                           '{days}d {hours}h {minutes}m {seconds}s separated by'
                           ' spaces e.g: 1d 2h 20m 15s,  10m 3s, 1h, 20s or 25.')
        return message

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
