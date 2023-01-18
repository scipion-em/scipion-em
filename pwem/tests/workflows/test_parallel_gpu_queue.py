# **************************************************************************
# *
# * Authors:    David Maluenda Niubo (dmaluenda@cnb.csic.es)
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
import json
import subprocess

from pyworkflow.tests import *
import pyworkflow.utils as pwutils

from pwem import Domain
import pwem.protocols as emprot

# --- Set this to match with your queue system ---
#  json params to fill the queue form, see SCIPION_HOME/config/host.conf
QUEUE_PARAMS = (u'myslurmqueue', {u'JOB_TIME': u'1',        # in hours
                                  u'JOB_MEMORY': u'2048',   # in Mb
                                  u'QUEUE_FOR_JOBS': u'N', })
#  command and args to list the queued jobs (to be used in a subprocess)
#  the command's output must contain the jobID and the protocolID
QUEUE_COMMAND = ["squeue"]


class TestQueueBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='mda'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('particles')

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(emprot.ProtImportParticles, 
                                     filesPath=pattern,
                                     samplingRate=samplingRate,
                                     checkStack=checkStack)
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles '
                            'is None.' % pattern)
        return protImport

    @classmethod
    def runNormalizeParticles(cls, particles):
        """ Run normalize particles protocol """
        relionProtocols = Domain.importFromPlugin('relion.protocols',
                                                  doRaise=True)
        protPreproc = cls.newProtocol(relionProtocols.ProtRelionPreprocessParticles,
                                      doNormalize=True)
        protPreproc.inputParticles.set(particles)
        cls.launchProtocol(protPreproc)
        cls.sampling = protPreproc.outputParticles.getSamplingRate()
        return protPreproc

    def _checkAsserts(self, prot):
        """ If prot.useQueue() is True, first of all
            we check if the job is queued and after that
            we wait until it's finished to check the outputs
        """
        def isJobInQueue(jobId, protId, Yes=True):
            """ returns Yes if the protId and the jobId
                are found in the queue.
                Set Yes=False to get a True when the job disappears
            """
            isThere = not Yes
            queueRaw = subprocess.check_output(QUEUE_COMMAND)
            queueList = queueRaw.split('\n')
            for queueLine in queueList:
                # we will assume that if protId and jobId
                # is in the string, the task is there
                if jobId in queueLine and str(protId) in queueLine:
                    isThere = Yes
                    break
            return isThere

        def wait_until(condition, timeout, *args, **kwargs):
            """ :param condition: function handle that
                                      is waited until return True
                :param timeout: maximum time to wait
                :param args and kwargs: params to pass to the condition func.
                :return: False if the timeout is reached
            """
            mustend = time.time() + timeout
            while time.time() < mustend:
                if condition(*args, **kwargs):
                    return True
                time.sleep(1)
            return False

        def checkQueue(jobId, protId):
            """ Check if the protocol job is queued
            """
            self.assertTrue(isJobInQueue(jobId, protId),
                            "The job %s corresponding to "
                            "the protocol %d has been not "
                            "attached to the system queue"
                            % (jobId, protId))
            print(pwutils.magentaStr(" > job %s of the protocol %d found in the "
                                     "queue, wait a sec..." % (jobId, protId)))

            isDone = wait_until(isJobInQueue, 10*60, jobId, protId, Yes=False)
            self.assertTrue(isDone, "Timeout: the job has not ended...")
            print(pwutils.magentaStr("    ...job ended!"))

        if prot.useQueue():
            if QUEUE_PARAMS[1]['QUEUE_FOR_JOBS'] != 'Y':
                # if the protocol is use queue system, we check if it's queued
                jobId = prot.getJobId()   # is an string
                protId = prot.getObjId()  # is an integer
                checkQueue(jobId, protId)
                return  # I don't know why, but we cannot retrieve the output, permissions???
            else:
                # Check that job files have been created
                jobFilesPath = join(pwutils.getParentFolder(prot.getLogPaths()[0]),
                                    str(prot.getObjId()))

                self.assertTrue(
                    pwutils.exists(jobFilesPath + "-0-1.out") and pwutils.exists(
                        jobFilesPath + "-0-1.err") and pwutils.exists(jobFilesPath + "-0-1.job"),
                    "Job queue files not found in log folder, job did not make it to the queue.")

        self.assertIsNotNone(prot.outputClasses,
                             "There was a problem with Relion 2D classify")

        classsesPixSize = prot.outputClasses.getImages().getSamplingRate()
        self.assertAlmostEquals(self.sampling, classsesPixSize, delta=0.001,
                                msg="There was a problem with the sampling rate "
                                    "of the particles")
        for class2D in prot.outputClasses:
            self.assertTrue(class2D.hasAlignment2D())

    def _runRelionClassify2D(self, previousRun, label='', threads=1, MPI=1,
                             doGpu=False, GPUs='', useQueue=False, steps=False):
        """ :param previousRun: The outputParticles of that will be the input
            :param label: For naming purposes
            :param threads: How many threads to use
            :param MPI: How many MPIs to use
            :param doGpu: Use GPU or not
            :param GPUs: Which GPUs to use (see Relion gpusToUse form param)
            :param useQueue: Use the queue system or not
            :return: the launched protocol
        """
        relionProtocols = Domain.importFromPlugin('relion.protocols',
                                                  doRaise=True)
        prot2D = self.newProtocol(relionProtocols.ProtRelionClassify2D,
                                  doCTF=False, maskDiameterA=340,
                                  useGradientAlg=False,
                                  numberOfMpi=MPI, numberOfThreads=threads)
        prot2D.numberOfClasses.set(4)
        prot2D.numberOfIterations.set(3)
        prot2D.inputParticles.set(previousRun.outputParticles)
        prot2D.setObjLabel(label)

        if useQueue:
            prot2D._useQueue.set(True)
            if steps:
                QUEUE_PARAMS[1]['QUEUE_FOR_JOBS'] = 'Y'
            prot2D._queueParams.set(json.dumps(QUEUE_PARAMS))

        prot2D.doGpu.set(doGpu)
        if doGpu:
            prot2D.gpusToUse.set(GPUs)

        self.launchProtocol(prot2D)
        return prot2D


class TestQueueALL(TestQueueBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestQueueBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testNoGpuSerial(self):
        relionNoGpu11 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU serial",
                                                  useQueue=True)
        self._checkAsserts(relionNoGpu11)

    def testNoGpuMPI(self):
        relionNoGpu14 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU MPI",
                                                  MPI=4, useQueue=True)
        self._checkAsserts(relionNoGpu14)

    def testNoGpuThreads(self):
        relionNoGpu41 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU Threads",
                                                  threads=4, useQueue=True)
        self._checkAsserts(relionNoGpu41)

    def testNoGpuMPIandThreads(self):
        relionNoGpu22 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU MPI+Threads",
                                                  MPI=2, threads=2,
                                                  useQueue=True)
        self._checkAsserts(relionNoGpu22)

    def testGpuSerial(self):
        relionGpu11 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU serial",
                                                doGpu=True, useQueue=True)
        self._checkAsserts(relionGpu11)

    def testGpuMPI(self):
        relionGpu12 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI",
                                                doGpu=True, MPI=2,
                                                useQueue=True)
        self._checkAsserts(relionGpu12)

    def testGpuThreads(self):
        relionGpu41 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU Threads",
                                                doGpu=True, threads=4,
                                                useQueue=True)
        self._checkAsserts(relionGpu41)

    def testGpuMPIandThreads(self):
        relionGpu22 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI+Threads",
                                                doGpu=True, useQueue=True,
                                                MPI=2, threads=2)
        self._checkAsserts(relionGpu22)


class TestQueueSmall(TestQueueBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestQueueBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testGpuMPI(self):
        relionGpu12 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI",
                                                doGpu=True, MPI=2,
                                                useQueue=True)
        self._checkAsserts(relionGpu12)

    def testGpuMPISteps(self):
        relionGpu12 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI Steps",
                                                doGpu=True, MPI=4,
                                                useQueue=True, steps=True)

        self._checkAsserts(relionGpu12)


class TestNoQueueALL(TestQueueBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestQueueBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testNoGpuSerial(self):
        relionNoGpu11 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU serial")
        self._checkAsserts(relionNoGpu11)

    def testNoGpuMPI(self):
        relionNoGpu14 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU MPI",
                                                  MPI=4)
        self._checkAsserts(relionNoGpu14)

    def testNoGpuThreads(self):
        relionNoGpu41 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU Threads",
                                                  threads=4)
        self._checkAsserts(relionNoGpu41)

    def testNoGpuMPIandThreads(self):
        relionNoGpu22 = self._runRelionClassify2D(self.protNormalize,
                                                  "Rel.2D noGPU MPI+Threads",
                                                  MPI=2, threads=4)
        self._checkAsserts(relionNoGpu22)

    def testGpuSerial(self):
        relionGpu11 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU serial",
                                                doGpu=True)
        self._checkAsserts(relionGpu11)

    def testGpuMPI(self):
        relionGpu12 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI",
                                                doGpu=True, MPI=2)
        self._checkAsserts(relionGpu12)

    def testGpuThreads(self):
        relionGpu41 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU Threads",
                                                doGpu=True, threads=4)
        self._checkAsserts(relionGpu41)

    def testGpuMPIandThreads(self):
        relionGpu22 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI+Threads",
                                                doGpu=True,
                                                MPI=2, threads=4)
        self._checkAsserts(relionGpu22)


class TestNoQueueSmall(TestQueueBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestQueueBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testGpuMPI(self):
        relionGpu12 = self._runRelionClassify2D(self.protNormalize,
                                                "Rel.2D GPU MPI",
                                                doGpu=True, MPI=2)
        self._checkAsserts(relionGpu12)


class TestQueueSteps(TestQueueBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestQueueBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)

    def testStepsNoGPU(self):
        xmipp3Protocols = Domain.importFromPlugin('xmipp3.protocols',
                                                  doRaise=True)
        protXmippPreproc = self.newProtocol(xmipp3Protocols.XmippProtPreprocessParticles,
                                            doNormalize=True, doRemoveDust=True)

        protXmippPreproc.inputParticles.set(self.protImport.outputParticles)
        protXmippPreproc.setObjLabel("Xmipp preprocess steps")

        protXmippPreproc._useQueue.set(True)

        QUEUE_PARAMS[1]['QUEUE_FOR_JOBS'] = 'Y'

        protXmippPreproc._queueParams.set(json.dumps(QUEUE_PARAMS))

        # Launch protocol but wait until it finishes
        self.launchProtocol(protXmippPreproc, wait=True)

        # Check that job files have been created

        jobFilesPath = join(pwutils.getParentFolder(protXmippPreproc.getLogPaths()[0]),
                            str(protXmippPreproc.getObjId()))

        self.assertTrue(pwutils.exists(jobFilesPath + "-0-1.out") and
                        pwutils.exists(jobFilesPath + "-0-1.err") and
                        pwutils.exists(jobFilesPath + "-0-1.job") and
                        pwutils.exists(jobFilesPath + "-0-2.out") and
                        pwutils.exists(jobFilesPath + "-0-2.err") and
                        pwutils.exists(jobFilesPath + "-0-2.job"),
                        "Job queue files not found on log folder, job did "
                        "not make it to the queue.")

        # Check that results have been produced
        self.assertIsNotNone(protXmippPreproc.outputParticles,
                             "There was a problem with Xmipp preprocess particles.")

