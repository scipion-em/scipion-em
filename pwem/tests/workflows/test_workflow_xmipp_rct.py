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
import pyworkflow.tests as pwtests

import pwem.protocols as emprot

from .test_workflow import TestWorkflow


# update this test when RCT workflow are implemented
class TestXmippRCTWorkflow(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.dataset = pwtests.DataSet.getDataSet('rct')
        cls.allCrdsDir = cls.dataset.getFile('positions')
        cls.micsUFn = cls.dataset.getFile('untilted')
        cls.micsTFn = cls.dataset.getFile('tilted')
        cls.classesSqlite = cls.dataset.getFile('classes')

    def testXmippRCTWorkflowBasic(self):
        # First, import a set of micrographs
        protImport = self.newProtocol(emprot.ProtImportMicrographsTiltPairs,
                                      patternUntilted=self.micsUFn,
                                      patternTilted=self.micsTFn,
                                      samplingRate=2.28, voltage=100,
                                      sphericalAberration=2.9)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographsTiltPair,
                             "There was a problem with the import")
        # self.validateFiles('protImportRCT', protImport)

        # Then simulate a particle picking
        print("Running fake particle picking...")
        from pwem import Domain
        XmippProtParticlePickingPairs = Domain.importFromPlugin(
            'xmipp3.protocols',
            'XmippProtParticlePickingPairs',
            doRaise=True)
        protPicking = self.newProtocol(XmippProtParticlePickingPairs,
                                       importFolder=self.allCrdsDir)

        protPicking.inputMicrographsTiltedPair.set(protImport.outputMicrographsTiltPair)

        self.proj.launchProtocol(protPicking, wait=True)

        self.assertIsNotNone(protPicking.outputCoordinatesTiltPair,
                             "There was a problem with the faked picking")
        # self.validateFiles('protPicking', protPicking)

        # Extract particles
        print("Run extract particles with Same as picking")
        XmippProtExtractParticlesPairs = Domain.importFromPlugin(
            'xmipp3.protocols',
            'XmippProtExtractParticlesPairs')
        SAME_AS_PICKING = Domain.importFromPlugin('xmipp3.constants',
                                                  'SAME_AS_PICKING')
        protExtract = self.newProtocol(XmippProtExtractParticlesPairs,
                                       downFactor=2,
                                       boxSize=60,
                                       doInvert=False,
                                       downsampleType=SAME_AS_PICKING)
        protExtract.inputCoordinatesTiltedPairs.set(protPicking.outputCoordinatesTiltPair)

        self.proj.launchProtocol(protExtract, wait=True)

        # self.validateFiles('protExtract', protExtract)
        self.assertIsNotNone(protExtract.outputParticlesTiltPair,
                             "There was a problem with the extract particles")

        # Classify using Xmipp CL2D
        print("Run CL2D")
        XmippProtCL2D = Domain.importFromPlugin('xmipp3.protocols',
                                                'XmippProtCL2D')
        protCL2D = self.newProtocol(XmippProtCL2D,
                                    numberOfClasses=10,
                                    numberOfInitialClasses=1,
                                    numberOfIterations=3, numberOfMpi=2)
        protCL2D.inputParticles.set(protExtract.outputParticlesTiltPair.getUntilted())
        self.launchProtocol(protCL2D)
        self.assertIsNotNone(protCL2D.outputClasses,
                             "There was a problem with CL2D")
        # self.validateFiles('protCL2D', protCL2D)

        # Random Conical Tilt
        print("Run Random Conical Tilt")
        XmippProtRCT = Domain.importFromPlugin('xmipp3.protocols',
                                               'XmippProtRCT')
        protRCT = self.newProtocol(XmippProtRCT)
        protRCT.inputParticlesTiltPair.set(protExtract.outputParticlesTiltPair)
        protRCT.inputParticles.set(protCL2D.outputClasses)

        self.proj.launchProtocol(protRCT, wait=True)

        # self.validateFiles('protExtract', protExtract)
        self.assertIsNotNone(protRCT.outputVolumes, "There was a problem with the RCT")
