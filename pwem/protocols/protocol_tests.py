# **************************************************************************
# *
# * Authors:     R. Marabini (roberto@cnb.csic.es)
# *              Pablo Conesa (pconesa@cnb.csic.es)
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
""" This module is to add protocols that need to be there for testing purposes
but nos available for regular users. This will hide them and prevent the their
 appearance in any GUI search or tree."""
import time

import pyworkflow.protocol as pwprot
import pyworkflow.protocol.params as params


STRESS_NG = 'stress-ng'


class ProtTests(pwprot.Protocol):
    @classmethod
    def isDisabled(cls):
        """ Return True for all test protocols.
        Disabled protocols will not be offered in the available protocols."""
        return True


class ProtStress(ProtTests):
    """ stress  will  stress  test  a  computer system in various selectable
       ways. Several options require the program stress-ng.
    """
    _label = 'stress'
    _program = STRESS_NG

    def __init__(self, **kwargs):
        ProtTests.__init__(self, **kwargs)
        self.xmippMic = {}

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('noCpu', params.IntParam, default=0,
                      label="No. CPU stressors",
                      help="start N workers spinning on sqrt(rand())")
        form.addParam('noMem', params.IntParam, default=0,
                      label="No. memory stressors",
                      help="start N workers continuously calling mmap(2)/munmap"
                           "(2) and writing to the allocated memory. It will "
                           "use 100% per worker ")
        form.addParam('amountMem', params.IntParam, default=256,
                      label="Memory per stressors (Mb)",
                      help="allocate N bytes per vm worker. each stressors will"
                           " use 100% of a CPU")
        form.addParam('noIO', params.IntParam, default=0,
                      label="No IO stressors",
                      help="start N workers spinning on sync() (Disk io)")
        form.addParam('timeout', params.IntParam, default=60,
                      label="TimeOut (sec)",
                      help="timeout after N seconds. Total execution time is "
                           "timeout plus delay")
        form.addParam('delay', params.IntParam, default=0,
                      label="delay (sec)",
                      help="wait this seconds before stressing the system "
                           "seconds")
        form.addParam('extraParams', params.StringParam,
                      label='Additional parameters', default=" ",
                      help='Additional parameters for stress-ng '
                           '(http://kernel.ubuntu.com/~cking/stress-ng/)')

    # ---------------------- INSERT steps functions ----------------------------

    def _insertAllSteps(self):
        if self.delay != 0:
            self._insertFunctionStep('delayStep')
        self._insertFunctionStep('stressStep')

    # ------------------------ STEPS functions ---------------------------------
    def delayStep(self):
        time.sleep(self.delay)

    def stressStep(self):
        if self.noCpu == 0 and self.noIO == 0 and self.noMem == 0:
            time.sleep(self.timeout)
        else:
            args = " "
            if self.noCpu != 0:
                args += " --cpu %d" % self.noCpu
            if self.noIO != 0:
                args += " --io %d" % self.noIO
            if self.noMem != 0:
                args += " --mmap %d" % self.noMem
            if self.amountMem != 0:
                args += " --mmap-bytes %dM" % self.amountMem
            args += " --timeout %ds" % self.timeout
            args += " --metrics-brief"
            if self.extraParams is not None:
                args += " %s " % self.extraParams
            self.runJob(self._program, args)

    # ------------------------ INFO functions ----------------------------------
    def _validate(self):
        message = []
        from shutil import which
        if self.noCpu == 0 and self.noIO == 0 and self.noMem == 0:
            pass
        else:
            if which(self._program) is None:
                message = ["Cannot find executable %s" % self._program]
        return message

    def _summary(self):
        message = "%d CPU, %d Mem, %d IO stressor" % (self.cpu, self.mem,
                                                      self.io)
        return [message]

    def _methods(self):
        return []


# class ProtOutputTest(ProtTests):
#     """ Protocol to test scalar output and input linking"""
#     _label = 'test output'
#
#     def __init__(self, **args):
#         Protocol.__init__(self, **args)
#         self.name = params.String(args.get('name', None))
#
#     def _defineParams(self, form):
#
#         section = form.addSection("Input")
#         section.addParam('iBoxSize', params.IntParam, allowsPointers=True,
#                          default=10,
#                          label='Input box size as Integer',
#                          validators=[params.Positive])
#
#         section.addParam('nullableInteger', params.IntParam, allowsPointers=True,
#                          label='Nullable Integer', allowsNull=True)
#
#     def _createOutputStep(self):
#         # New Output would be an Integer
#         boxSize = Integer(10)
#
#         if self.iBoxSize.hasValue():
#             boxSize.set(2*int(self.iBoxSize.get()))
#
#         self._defineOutputs(oBoxSize=boxSize)
#
#     def _insertAllSteps(self):
#         self._insertFunctionStep('_createOutputStep')
#
#
# class ProtMultiPointerTest(ProtTests):
#     """ Class to test how multipointer params are exported to json"""
#     def _defineParams(self, form):
#
#         # This should cover Multipointer params that points to attributes...
#         # therefore extended attribute of pointers should be used
#         form.addParam('mpToAttr', params.MultiPointerParam,
#                       label="Multipointer to attribute",
#                       pointerClass='String',
#                       help="Should point to String inside another protocol")
#
#         # This should cover Multipointer params that points to protocols...
#         # therefore extended attribute of pointers should NOT be used
#         form.addParam('mpToProts', params.MultiPointerParam,
#                       label="Multipointer to sets",
#                       pointerClass='Protocol',
#                       help="Should point to another protocol")

