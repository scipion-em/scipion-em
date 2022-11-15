# **************************************************************************
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
import  logging
logger = logging.getLogger(__name__)

import pyworkflow.tests as pwtests

from pwem import Domain
import pwem.protocols as emprot

from .test_workflow import TestWorkflow


class TestMixedBPV(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.dataset = pwtests.DataSet.getDataSet('xmipp_tutorial')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')

    def test_workflow(self):
        try:
            from itertools import izip
        except ImportError:
            izip = zip
        # First, import a set of micrographs
        protImport = self.newProtocol(emprot.ProtImportMicrographs,
                                      filesPath=self.micsFn,
                                      samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")
        #         self.validateFiles('protImport', protImport)

        # Import a set of volumes
        logger.info("Import Volume")
        protImportVol = self.newProtocol(emprot.ProtImportVolumes, filesPath=self.vol1,
                                         samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(),
                             "There was a problem with the import")
        #        self.validateFiles('protImportVol', protImportVol)

        # Perform a downsampling on the micrographs
        logger.info("Downsampling...")
        XmippProtPreprocessMicrographs = Domain.importFromPlugin('xmipp3.protocols',
                                                                 'XmippProtPreprocessMicrographs',
                                                                 doRaise=True)
        protDownsampling = self.newProtocol(XmippProtPreprocessMicrographs,
                                            doDownsample=True, downFactor=5,
                                            doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertIsNotNone(protDownsampling.outputMicrographs,
                             "There was a problem with the downsampling")
        #         self.validateFiles('protDownsampling', protDownsampling)

        # Estimate CTF on the downsampled micrographs
        logger.info("Performing CTFfind...")
        ProtCTFFind = Domain.importFromPlugin('cistem.protocols',
                                              'CistemProtCTFFind', doRaise=True)
        protCTF = self.newProtocol(ProtCTFFind, numberOfThreads=4,
                                   minDefocus=22000, maxDefocus=25000)
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protCTF)
        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem with the CTF estimation")

        valuesList = [[24000, 24000], [22548, 22518], [23058, 23029]]
        for ctfModel, values in izip(protCTF.outputCTF, valuesList):
            self.assertAlmostEquals(ctfModel.getDefocusU(), values[0], delta=1000)
            self.assertAlmostEquals(ctfModel.getDefocusV(), values[1], delta=1000)
            self.assertAlmostEquals(ctfModel.getMicrograph().getSamplingRate(),
                                    6.185, delta=0.001)

        #         self.validateFiles('protCTF', protCTF)

        logger.info("Running Eman import coordinates...")
        protPP = self.newProtocol(emprot.ProtImportCoordinates,
                                  importFrom=emprot.ProtImportCoordinates.IMPORT_FROM_EMAN,
                                  filesPath=self.crdsDir,
                                  filesPattern='*_info.json', boxSize=110)
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protPP)
        self.assertIsNotNone(protPP.outputCoordinates,
                             "There was a problem with the Eman import coordinates")

        # Extract the SetOfParticles.
        logger.info("Run extract particles with other downsampling factor")
        XmippProtExtractParticles = Domain.importFromPlugin(
            'xmipp3.protocols',
            'XmippProtExtractParticles')
        OTHER = Domain.importFromPlugin('xmipp3.constants', 'OTHER')
        protExtract = self.newProtocol(XmippProtExtractParticles,
                                       boxSize=64,
                                       downsampleType=OTHER,
                                       downFactor=8.0,
                                       doFlip=False, runMode=1, doInvert=False)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles,
                             "There was a problem with the extract particles")
        #         self.validateFiles('protExtract', protExtract)


class TestMixedBPV2(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.dataset = pwtests.DataSet.getDataSet('xmipp_tutorial')
        cls.crdsDir = cls.dataset.getFile('boxingDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')

    def test_workflow(self):
        # First, import a set of micrographs
        protImport = self.newProtocol(emprot.ProtImportMicrographs,
                                      filesPath=self.micsFn,
                                      samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")
        #         self.validateFiles('protImport', protImport)

        # Import a set of volumes
        logger.info("Import Volume")
        protImportVol = self.newProtocol(emprot.ProtImportVolumes, filesPath=self.vol1,
                                         samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(),
                             "There was a problem with the import")
        #        self.validateFiles('protImportVol', protImportVol)

        # Perform a downsampling on the micrographs
        logger.info("Downsampling...")
        XmippProtPreprocessMicrographs = Domain.importFromPlugin(
            'xmipp3.protocols',
            'XmippProtPreprocessMicrographs',
            doRaise=True)
        protDownsampling = self.newProtocol(XmippProtPreprocessMicrographs,
                                            doDownsample=True, downFactor=5,
                                            doCrop=False, runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertIsNotNone(protDownsampling.outputMicrographs,
                             "There was a problem with the downsampling")
        #         self.validateFiles('protDownsampling', protDownsampling)

        # Estimate CTF on the downsampled micrographs
        logger.info("Performing CTFfind...")
        ProtCTFFind = Domain.importFromPlugin('cistem.protocols',
                                              'CistemProtCTFFind', doRaise=True)
        protCTF = self.newProtocol(ProtCTFFind, numberOfThreads=4,
                                   minDefocus=22000, maxDefocus=25000)
        protCTF.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protCTF)
        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem with the CTF estimation")
        #         self.validateFiles('protCTF', protCTF)

        logger.info("Running Eman import coordinates...")
        protPP = self.newProtocol(emprot.ProtImportCoordinates,
                                  importFrom=emprot.ProtImportCoordinates.IMPORT_FROM_EMAN,
                                  filesPath=self.crdsDir,
                                  filesPattern='*_info.json', boxSize=110)
        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protPP)
        self.assertIsNotNone(protPP.outputCoordinates,
                             "There was a problem with the Eman import coordinates")

        logger.info("<Run extract particles with Same as picking>")
        XmippProtExtractParticles = Domain.importFromPlugin(
            'xmipp3.protocols',
            'XmippProtExtractParticles')

        SAME_AS_PICKING = Domain.importFromPlugin('xmipp3.constants',
                                                  'SAME_AS_PICKING')
        protExtract = self.newProtocol(XmippProtExtractParticles, boxSize=110,
                                       downsampleType=SAME_AS_PICKING,
                                       doFlip=True, doInvert=True, runMode=1)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles,
                             "There was a problem with the extract particles")
        # self.validateFiles('protExtract', protExtract)

        logger.info("Run Preprocess particles")
        XmippProtCropResizeParticles = Domain.importFromPlugin(
            'xmipp3.protocols',
            'XmippProtCropResizeParticles')
        protCropResize = self.newProtocol(XmippProtCropResizeParticles,
                                          doResize=True, resizeOption=1,
                                          resizeDim=110)
        protCropResize.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protCropResize)

        self.assertIsNotNone(protCropResize.outputParticles,
                             "There was a problem with resize/crop the particles")

        logger.info("Run ML2D")
        XmippProtML2D = Domain.importFromPlugin('xmipp3.protocols',
                                                'XmippProtML2D')
        protML2D = self.newProtocol(XmippProtML2D, numberOfClasses=8, maxIters=2,
                                    numberOfMpi=2, numberOfThreads=2)
        protML2D.inputParticles.set(protCropResize.outputParticles)
        self.launchProtocol(protML2D)
        self.assertIsNotNone(protML2D.outputClasses,
                             "There was a problem with ML2D")
        # self.validateFiles('protML2D', protML2D)

        #         #FIXME: Check the initial model with EMAn and restore the next step
        #         return

        logger.info("Run Initial Model")

        EmanProtInitModel = Domain.importFromPlugin('eman2.protocols',
                                                    'EmanProtInitModel',
                                                    doRaise=True)
        protIniModel = self.newProtocol(EmanProtInitModel, numberOfIterations=1,
                                        numberOfModels=2,
                                        shrink=5, symmetry='icos',
                                        numberOfThreads=3)
        protIniModel.inputSet.set(protML2D.outputClasses)
        self.launchProtocol(protIniModel)
        self.assertIsNotNone(protIniModel.outputVolumes,
                             "There was a problem with Initial Model")
