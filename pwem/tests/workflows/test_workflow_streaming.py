# ***************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] Science for Life Laboratory, Stockholm University
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
# ***************************************************************************/

import time
import os
from glob import glob
import threading

import pyworkflow.utils as pwutils
import pyworkflow.tests as pwtests

from pwem import Domain, emlib
import pwem.protocols as emprot
import relion


# Load the number of movies for the simulation, by default equal 5, but
# can be modified in the environement

def _getVar(varSuffix, varType, default=None):
    return varType(os.environ.get('SCIPION_TEST_STREAM_%s' % varSuffix, default))


MOVS = _getVar('MOVS', int, 10)
PATTERN = _getVar('PATTERN', str, '')
DELAY = _getVar('DELAY', int, 10)  # in seconds
# Change the timeout to stop waiting for new files
TIMEOUT = _getVar('TIMEOUT', int, 60)


class TestStreamingWorkflow(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('movies')
        cls.importThread = threading.Thread(name="createInputLinks",
                                            target=cls._createInputLinks)
        cls.importThread.start()
        # Wait until the first link is created
        time.sleep(5)

    @classmethod
    def _createInputLinks(cls):
        # Create a test folder path
        pattern = PATTERN if PATTERN else cls.ds.getFile('ribo/Falcon*mrcs')
        files = glob(pattern)
        nFiles = len(files)
        nMovies = MOVS

        for i in range(nMovies):
            # Loop over the number of input movies if we want more for testing
            f = files[i % nFiles]
            _, cls.ext = os.path.splitext(f)
            moviePath = cls.proj.getTmpPath('movie%06d%s' % (i + 1, cls.ext))
            pwutils.createAbsLink(f, moviePath)
            time.sleep(DELAY)

    def finalChecks(self, importMovies, protocols):
        """ Checks to do:
             - All movies are imported
             - Last processing protocol have processed all data
             - Summary Monitor is finished after a while
             - All protocols are finished at the end
        """
        lastProt = protocols[-2]  # last processing protocol (CTF or whatever)
        monitorProt = protocols[-1]  # last protocol usually monitor summary

        def _loadProt(prot):
            # Load the last version of the protocol from its own database
            prot2 = pwtests.getProtocolFromDb(prot.getProject().path,
                                              prot.getDbPath(),
                                              prot.getObjId())
            # Close DB connections
            prot2.getProject().closeMapper()
            prot2.closeMappers()
            return prot2

        def importFinished():
            importMovies.load()
            importSize = importMovies.getSize()
            return importSize == MOVS

        def protocolsFinished():
            lastOutput = ([getattr(lastProt, k) for k, v
                           in lastProt.iterOutputAttributes()])[-1]
            lastOutput.load()
            return lastOutput.getSize() == MOVS

        t0 = time.time()
        precessingTimeout = 15  # minutes (lastProt must process all data, yet)
        while not importFinished() or not protocolsFinished():
            if time.time() - t0 > precessingTimeout * 60:
                self.assertTrue(importFinished(),
                                "Some error occurred with the acquisition. "
                                "Not all Movies have been imported after %dmin."
                                % precessingTimeout)
                self.assertTrue(protocolsFinished(),
                                "Some error occurred with some protocol. "
                                "The outputSize of the last protocol is different "
                                "than the imported movies after %dmin."
                                % precessingTimeout)
            time.sleep(5)

        counter = 1
        protMon2 = _loadProt(monitorProt)
        waitingTimeout = 2  # minutes (summary monitor is refreshed every 60s)
        while protMon2.isActive() and counter < waitingTimeout * 10:
            time.sleep(6)
            counter += 1
            protMon2 = _loadProt(monitorProt)

        for prot1 in protocols:
            prot2 = _loadProt(prot1)
            self.assertTrue(_loadProt(prot1).isFinished(),
                            "%s has not finished after the processing ends "
                            "(i.e. %dmin after %s has processed all input files)."
                            % (prot2.getObjLabel(), waitingTimeout,
                               lastProt.getObjLabel()))

    def test_pattern(self):
        protocols = []  # append here all protocols to make final checks
        # ----------- IMPORT MOVIES -------------------
        protImport = self.newProtocol(emprot.ProtImportMovies,
                                      objLabel='import movies',
                                      importFrom=emprot.ProtImportMovies.IMPORT_FROM_FILES,
                                      filesPath=os.path.abspath(self.proj.getTmpPath()),
                                      filesPattern="movie*%s" % self.ext,
                                      amplitudConstrast=0.1,
                                      sphericalAberration=2.,
                                      voltage=300,
                                      samplingRate=3.54,
                                      dataStreaming=True,
                                      timeout=TIMEOUT)
        self.proj.launchProtocol(protImport, wait=False)
        self._waitOutput(protImport, 'outputMovies')
        protocols.append(protImport)

        # ----------- ALIGNMENT --------------------------
        XmippProtMovieCorr = Domain.importFromPlugin('xmipp3.protocols',
                                                     'XmippProtMovieCorr',
                                                     doRaise=True)
        protOF = self.newProtocol(XmippProtMovieCorr,
                                  objLabel='Movie alignment',
                                  doSaveMovie=False,
                                  doComputePSD=False,
                                  alignFrame0=3,
                                  alignFrameN=10,
                                  sumFrame0=3,
                                  sumFrameN=10,
                                  doApplyDoseFilter=False)
        protOF.inputMovies.set(protImport)
        protOF.inputMovies.setExtended('outputMovies')
        self.proj.launchProtocol(protOF, wait=False)
        self._waitOutput(protOF, 'outputMicrographs')
        protocols.append(protOF)

        # --------- CTF ESTIMATION ---------------------------
        ProtCTFFind = Domain.importFromPlugin('cistem.protocols',
                                              'CistemProtCTFFind', doRaise=True)
        protCTF = self.newProtocol(ProtCTFFind,
                                   objLabel='ctffind4')
        protCTF.inputMicrographs.set(protOF)
        protCTF.inputMicrographs.setExtended('outputMicrographs')
        self.proj.launchProtocol(protCTF, wait=False)
        self._waitOutput(protCTF, 'outputCTF')
        protocols.append(protCTF)

        # --------- SUMMARY MONITOR --------------------------
        protMonitor = self.newProtocol(emprot.ProtMonitorSummary,
                                       objLabel='summary')
        protMonitor.inputProtocols.append(protImport)
        protMonitor.inputProtocols.append(protOF)
        protMonitor.inputProtocols.append(protCTF)
        self.proj.launchProtocol(protMonitor, wait=False)
        protocols.append(protMonitor)

        # ---------- STREAMING CHECKS ------------------------
        self.finalChecks(protImport.outputMovies, protocols)


class TestBaseRelionStreaming(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion_tutorial')
        cls.importThread = threading.Thread(name="createInputLinksR",
                                            target=cls._createInputLinks)
        cls.importThread.start()
        # Wait until the first link is created
        time.sleep(5)

    @classmethod
    def _createInputLinks(cls):
        # Create a test folder path
        pattern = cls.ds.getFile('allMics')
        files = glob(pattern)
        nFiles = len(files)

        for i in range(nFiles):
            # Loop over the number of input movies if we want more for testing
            f = files[i % nFiles]
            _, cls.ext = os.path.splitext(f)
            moviePath = cls.proj.getTmpPath('movie%06d%s' % (i + 1, cls.ext))
            pwutils.createAbsLink(f, moviePath)
            time.sleep(10)

    def _waitUntilMinSize(self, emSet, size=5):
        counter = 0
        while not emSet.getSize() >= size:
            time.sleep(5)
            emSet.load()
            if counter > 1000:
                break
            counter += 1


class TestRelionExtractStreaming(TestBaseRelionStreaming):
    def testRisosome(self):
        # First, import a set of micrographs
        print("Importing a set of micrographs...")
        protImport = self.newProtocol(emprot.ProtImportMicrographs,
                                      filesPath=os.path.abspath(self.proj.getTmpPath()),
                                      filesPattern="*%s" % self.ext,
                                      samplingRateMode=1,
                                      magnification=79096,
                                      scannedPixelSize=56, voltage=300,
                                      sphericalAberration=2.0,
                                      dataStreaming=True,
                                      fileTimeout=3,
                                      timeout=60)
        protImport.setObjLabel('import 20 mics (streaming)')
        self.proj.launchProtocol(protImport, wait=False)
        self._waitOutput(protImport, 'outputMicrographs')

        # Now estimate CTF on the micrographs with ctffind
        print("Performing CTFfind...")
        ProtCTFFind = Domain.importFromPlugin('cistem.protocols',
                                              'CistemProtCTFFind', doRaise=True)
        protCTF = self.newProtocol(ProtCTFFind,
                                   minDefocus=12000, maxDefocus=30000,
                                   runMode=1,
                                   numberOfMpi=1, numberOfThreads=1)
        protCTF.inputMicrographs.set(protImport.outputMicrographs)
        protCTF.setObjLabel('CTF ctffind')
        self.proj.launchProtocol(protCTF, wait=False)
        self._waitOutput(protCTF, 'outputCTF')

        # Now pick particles on the micrographs with Relion
        print("Performing Relion Autopicking (LoG)...")

        ProtRelionAutopickLoG = Domain.importFromPlugin(
            'relion.protocols', 'ProtRelionAutopickLoG')

        protPick = self.newProtocol(
            ProtRelionAutopickLoG,
            objLabel='Streaming relion - autopic (LoG)',
            inputMicrographs=protImport.outputMicrographs,
            boxSize=60,
            minDiameter=260,
            maxDiameter=380)
        self.proj.launchProtocol(protPick, wait=False)
        self._waitOutput(protPick, 'outputCoordinates')

        ProtRelionExtractParticles = Domain.importFromPlugin(
            'relion.protocols',
            'ProtRelionExtractParticles',
            doRaise=True)
        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       objLabel='extract box=64',
                                       boxSize=64,
                                       doInvert=True
                                       )
        protExtract.inputCoordinates.set(protPick.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        self.launchProtocol(protExtract)

        x, y, _ = protExtract.outputParticles.getDim()
        self.assertEqual(protExtract.outputParticles.getDim(), (64, 64, 1),
                         "Dimension of the particles should be 64 x 64 and it "
                         "is %d x %d" % (x, y))


class TestRelionPickStreaming(TestBaseRelionStreaming):
    def testRisosome(self):
        print("Importing 2D averages (subset of 4)")
        ih = emlib.image.ImageHandler()
        classesFn = self.ds.getFile('import/classify2d/extra/'
                                    'relion_it015_classes.mrcs')

        outputName = 'input_averages.mrcs'
        inputTmp = os.path.abspath(self.proj.getTmpPath())
        outputFn = self.proj.getTmpPath(outputName)

        for i, index in enumerate([5, 16, 17, 18, 24]):
            ih.convert((index, classesFn), (i + 1, outputFn))

        protAvgs = self.newProtocol(emprot.ProtImportAverages,
                                    objLabel='avgs - 5',
                                    filesPath=inputTmp,
                                    filesPattern=outputName,
                                    samplingRate=7.08
                                    )
        self.launchProtocol(protAvgs)

        # First, import a set of micrographs
        print("Importing a set of micrographs...")
        protImport = self.newProtocol(emprot.ProtImportMicrographs,
                                      objLabel='import 20 mics (streaming)',
                                      filesPath=os.path.abspath(self.proj.getTmpPath()),
                                      filesPattern="*%s" % self.ext,
                                      samplingRateMode=1,
                                      magnification=79096,
                                      scannedPixelSize=56, voltage=300,
                                      sphericalAberration=2.0,
                                      dataStreaming=True,
                                      fileTimeout=10,
                                      timeout=60)
        self.proj.launchProtocol(protImport, wait=False)
        self._waitOutput(protImport, 'outputMicrographs')

        # Now estimate CTF on the micrographs with ctffind
        print("Performing CTFfind...")
        ProtCTFFind = Domain.importFromPlugin('cistem.protocols',
                                              'CistemProtCTFFind', doRaise=True)
        protCtf = self.newProtocol(
            ProtCTFFind,
            objLabel='ctffind',
            inputMicrographs=protImport.outputMicrographs,
            minDefocus=12000, maxDefocus=30000,
            slowSearch=False
        )
        self.proj.launchProtocol(protCtf, wait=False)
        self._waitOutput(protCtf, 'outputCTF')

        self._waitUntilMinSize(protCtf.outputCTF)

        # Select some good averages from the iterations mrcs a

        ProtRelion2Autopick = Domain.importFromPlugin('relion.protocols',
                                                      'ProtRelion2Autopick')
        protPick = self.newProtocol(
            ProtRelion2Autopick, objLabel='autopick refs',
            inputMicrographs=protImport.outputMicrographs,
            ctfRelations=protCtf.outputCTF,
            runType=relion.RUN_COMPUTE,
            inputReferences=protAvgs.outputAverages
        )
        self.launchProtocol(protPick)


class TestFrameStacking(pwtests.BaseTest):
    """ Test the cases where the input movies are input as individual frames.
    """

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('movies')

    @classmethod
    def _createFrames(cls, delay=0):
        # Create a test folder path
        pattern = cls.ds.getFile('ribo/Falcon*mrcs')
        files = glob(pattern)

        nFiles = len(files)
        nMovies = MOVS
        ih = emlib.image.ImageHandler()

        for i in range(nMovies):
            # Loop over the number of input movies if we want more for testing
            f = files[i % nFiles]
            _, _, _, nFrames = ih.getDimensions(f)

            for j in range(1, nFrames + 1):
                outputFramePath = cls.proj.getTmpPath('movie%06d_%03d.mrc'
                                                      % (i + 1, j))
                ih.convert((j, f), outputFramePath)
                time.sleep(delay)

    def test_noStream(self):
        """ Test that we can input individual frames even if not
        processing in streaming.
        """
        self._createFrames()

        # ----------- IMPORT MOVIES -------------------
        protImport = self.newProtocol(emprot.ProtImportMovies,
                                      objLabel='import stack no stream',
                                      importFrom=emprot.ProtImportMovies.IMPORT_FROM_FILES,
                                      filesPath=os.path.abspath(self.proj.getTmpPath()),
                                      filesPattern="movie*.mrc",
                                      amplitudConstrast=0.1,
                                      sphericalAberration=2.,
                                      voltage=300,
                                      samplingRate=3.54,
                                      dataStreaming=False,
                                      inputIndividualFrames=True,
                                      numberOfIndividualFrames=16,
                                      stackFrames=True,
                                      writeMoviesInProject=True,
                                      deleteFrames=True)
        self.launchProtocol(protImport)
        self.assertSetSize(protImport.outputMovies, MOVS, msg="Wrong output set size!!")

    def test_Stream(self):
        # Create a separated thread to simulate real streaming with
        # individual frames
        thread = threading.Thread(name="createFrames",
                                  target=lambda: self._createFrames(delay=1))
        thread.start()
        time.sleep(5)

        # ----------- IMPORT MOVIES -------------------
        protImport = self.newProtocol(emprot.ProtImportMovies,
                                      objLabel='import stack streaming',
                                      importFrom=emprot.ProtImportMovies.IMPORT_FROM_FILES,
                                      filesPath=os.path.abspath(self.proj.getTmpPath()),
                                      filesPattern="movie*.mrc",
                                      amplitudConstrast=0.1,
                                      sphericalAberration=2.,
                                      voltage=300,
                                      samplingRate=3.54,
                                      dataStreaming=True,
                                      timeout=60,
                                      fileTimeout=5,
                                      inputIndividualFrames=True,
                                      numberOfIndividualFrames=16,
                                      stackFrames=True,
                                      writeMoviesInProject=True,
                                      deleteFrames=False)
        self.launchProtocol(protImport)
        self.assertSetSize(protImport.outputMovies, MOVS, msg="Wrong output set size!!")
        thread.join()
