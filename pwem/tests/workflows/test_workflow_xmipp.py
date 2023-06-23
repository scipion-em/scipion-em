# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
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

import numpy as np

import pyworkflow.tests as pwtests
from pwem import Domain
import pwem.protocols as emprot

from .test_workflow import TestWorkflow


class TestXmippWorkflow(TestWorkflow):

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.dataset = pwtests.DataSet.getDataSet('xmipp_tutorial')
        cls.allCrdsDir = cls.dataset.getFile('posAllDir')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.vol1 = cls.dataset.getFile('vol1')

    def testXmippWorkflow(self):
        # First, import a set of micrographs
        protImport = self.newProtocol(emprot.ProtImportMicrographs,
                                      filesPath=self.micsFn,
                                      samplingRate=1.237, voltage=300)
        self.launchProtocol(protImport)
        self.assertSetSize(protImport.outputMicrographs, 3)
        self.assertIsNotNone(protImport.outputMicrographs.getFileName(),
                             "There was a problem with the import")
        self.validateFiles('protImport', protImport)

        # Import a set of volumes
        print("Import Volume")
        protImportVol = self.newProtocol(emprot.ProtImportVolumes,
                                         filesPath=self.vol1,
                                         samplingRate=9.896)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(),
                             "There was a problem with the import")
        #        self.validateFiles('protImportVol', protImportVol)

        # Perform a downsampling on the micrographs
        print("Downsampling...")
        xmippProtcols = Domain.importFromPlugin('xmipp3.protocols',
                                                doRaise=True)
        protDownsampling = self.newProtocol(xmippProtcols.XmippProtPreprocessMicrographs,
                                            doDownsample=True,
                                            downFactor=5,
                                            doCrop=False,
                                            runMode=1)
        protDownsampling.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protDownsampling)
        self.assertSetSize(protDownsampling.outputMicrographs, 3)
        self.assertIsNotNone(protDownsampling.outputMicrographs,
                             "There was a problem with the downsampling")
        self.validateFiles('protDownsampling', protDownsampling)

        print("Importing coordinates")
        protPP = self.newProtocol(emprot.ProtImportCoordinates,
                                  importFrom=emprot.ProtImportCoordinates.IMPORT_FROM_XMIPP,
                                  filesPath=self.allCrdsDir,
                                  filesPattern='*.pos', boxSize=110)

        protPP.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protPP)
        self.protDict['protPicking'] = protPP
        self.assertSetSize(protPP.outputCoordinates, 143 + 138 + 122 + 110 + 138 + 122,
                           "There was a problem with the import of coordinates")

        # Now estimate CTF on the downsampled micrographs
        print("Performing CTF...")
        protCTF = self.newProtocol(xmippProtcols.XmippProtCTFMicrographs,
                                   numberOfThreads=4,
                                   minDefocus=2.2,
                                   maxDefocus=2.5)  # Defocus is in microns
        protCTF.inputMicrographs.set(protDownsampling.outputMicrographs)
        self.launchProtocol(protCTF)
        self.assertSetSize(protCTF.outputCTF, 3)
        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem with the CTF estimation")
        # After CTF estimation, the output micrograph should have CTF info
        self.validateFiles('protCTF', protCTF)

        print("Run extract particles with other downsampling factor")
        protExtract = self.newProtocol(xmippProtcols.XmippProtExtractParticles,
                                       boxSize=64, downsampleType=emprot.OTHER,
                                       doFlip=True, downFactor=8,
                                       runMode=1, doInvert=True)
        protExtract.inputCoordinates.set(protPP.outputCoordinates)
        protExtract.ctfRelations.set(protCTF.outputCTF)
        protExtract.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protExtract)
        self.assertIsNotNone(protExtract.outputParticles,
                             "There was a problem with the extract particles")
        self.validateFiles('protExtract', protExtract)

        print("Run Extract Coordinates without applying shifts")
        protExtractCoords = self.newProtocol(emprot.ProtExtractCoords)
        protExtractCoords.inputParticles.set(protExtract.outputParticles)
        protExtractCoords.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protExtractCoords)
        # The size of the set of coordinates must be the same as the input set of particles
        self.assertSetSize(protExtractCoords.outputCoordinates,
                           size=protExtract.outputParticles.getSize(),
                           msg="There was a problem with the coordinates extraction")
        # Check if the scaling factor is being calculated and applied correctly
        scale = protExtract.outputParticles.getSamplingRate() / protImport.outputMicrographs.getSamplingRate()
        inParticleCoord = protExtract.outputParticles.getFirstItem().getCoordinate()
        x, y = inParticleCoord.getPosition()
        np.testing.assert_allclose(protExtractCoords.outputCoordinates.getFirstItem().getPosition(),
                                   (int(x * scale), int(y * scale)), rtol=0.5)

        print("Run Screen Particles")
        protScreen = self.newProtocol(xmippProtcols.XmippProtScreenParticles,
                                      autoParRejection=xmippProtcols.XmippProtScreenParticles.REJ_MAXZSCORE,
                                      maxZscore=3.0)
        protScreen.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protScreen)
        self.assertSetSize(protScreen.outputParticles,
                           msg="There was a problem with Screen Particles")

        print("Run ML2D")
        protML2D = self.newProtocol(xmippProtcols.XmippProtML2D,
                                    numberOfClasses=1,
                                    maxIters=4,
                                    doMlf=True,
                                    numberOfMpi=2,
                                    numberOfThreads=1)
        protML2D.inputParticles.set(protScreen.outputParticles)
        self.launchProtocol(protML2D)
        self.assertIsNotNone(protML2D.outputClasses,
                             "There was a problem with ML2D")
        self.validateFiles('protML2D', protML2D)

        print("Run CL2D")
        protCL2D = self.newProtocol(xmippProtcols.XmippProtCL2D,
                                    numberOfClasses=2,
                                    numberOfInitialClasses=1,
                                    numberOfIterations=4, numberOfMpi=2)
        protCL2D.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protCL2D)
        self.assertIsNotNone(protCL2D.outputClasses, "There was a problem with CL2D")
        self.validateFiles('protCL2D', protCL2D)

        print("Run Only Align2d")
        protOnlyAlign = self.newProtocol(xmippProtcols.XmippProtCL2DAlign,
                                         maximumShift=5, numberOfIterations=2,
                                         numberOfMpi=2, numberOfThreads=1, useReferenceImage=False)
        protOnlyAlign.inputParticles.set(protExtract.outputParticles)
        self.launchProtocol(protOnlyAlign)
        self.assertIsNotNone(protOnlyAlign.outputParticles, "There was a problem with Only align2d")
        self.validateFiles('protOnlyAlign', protOnlyAlign)

        print("Run Extract Coordinates applying shifts")
        protExtractCoordsShifts = self.newProtocol(emprot.ProtExtractCoords, applyShifts=True)
        protExtractCoordsShifts.setObjLabel('Extract Coordinates applying shifts')
        inputParticles = protOnlyAlign.outputParticles
        inputMics = protDownsampling.outputMicrographs
        protExtractCoordsShifts.inputParticles.set(inputParticles)
        protExtractCoordsShifts.inputMicrographs.set(inputMics)
        self.launchProtocol(protExtractCoordsShifts)
        # The size of the set of coordinates must be the same as the input set of particles
        outputParticles = protExtractCoordsShifts.outputCoordinates
        self.assertSetSize(outputParticles,
                           size=protExtract.outputParticles.getSize(),
                           msg="There was a problem with the coordinates extraction")
        # Check if the scaling factor is being calculated and applied correctly
        scale = inputParticles.getSamplingRate() / protDownsampling.outputMicrographs.getSamplingRate()
        inParticleCoord = inputParticles.getFirstItem().getCoordinate()
        shifts = protExtractCoordsShifts.getShifts(inputParticles.getFirstItem().getTransform(),
                                                   inputParticles.getAlignment())
        x, y = inParticleCoord.getPosition()
        xCoor, yCoor = x - int(shifts[0]), y - int(shifts[1])
        np.testing.assert_allclose(protExtractCoordsShifts.outputCoordinates.getFirstItem().getPosition(),
                                   (int(xCoor * scale), int(yCoor * scale)), rtol=0.5)

        print("Run kerdensom")
        ProtKerdensom = self.newProtocol(xmippProtcols.XmippProtKerdensom,
                                         useMask=False, SomXdim=2, SomYdim=2,
                                         SomReg0=800, SomReg1=400, SomSteps=2)
        ProtKerdensom.inputParticles.set(protOnlyAlign.outputParticles)
        self.launchProtocol(ProtKerdensom)
        self.assertIsNotNone(ProtKerdensom.outputClasses, "There was a problem with kerdensom")
        # self.validateFiles('ProtKerdensom', ProtKerdensom)



if __name__ == "__main__":
    pwtests.unittest.main()
