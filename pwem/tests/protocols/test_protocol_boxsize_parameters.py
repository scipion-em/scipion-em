# ***************************************************************************
# * Authors:    Daniel Marchán (da.marchan@cnb.csic.es)
# *
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
# ***************************************************************************/
from pyworkflow.object import Pointer
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols.protocol_import import ProtImportMicrographs
from pwem.protocols.protocol_boxsize_parameters import ProtBoxSizeParameters


class TestBoxSizeParameters(BaseTest):
    """ Test box size parameters protocol"""

    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('relion_tutorial')

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        # Run needed protocols
        cls.runImportMicrograph()
        cls.runBoxSizeExtraction()

    @classmethod
    def runImportMicrograph(cls):

        """ Run an Import micrograph protocol. """
        protImport = cls.newProtocol(
            ProtImportMicrographs,
            samplingRateMode=0,
            filesPath=TestBoxSizeParameters.ds.getFile('micrographs/*.mrc'),
            samplingRate=3.54,
            magnification=59000,
            voltage=300,
            sphericalAberration=2)

        cls.launchProtocol(protImport)
        cls.protImport = protImport

    @classmethod
    def runBoxSizeExtraction(cls):
        """ Run box size related parameters on extraction boxSize. """
        protBoxSizeExtraction = cls.newProtocol(ProtBoxSizeParameters,
                                                inputMicrographs=cls.protImport.outputMicrographs,
                                                boxSize=46,
                                                boolExtractPartBx=True)
        protBoxSizeExtraction.inputMicrographs.set(cls.protImport.outputMicrographs)
        cls.launchProtocol(protBoxSizeExtraction)
        cls.protBoxSizeExtraction = protBoxSizeExtraction


    def testBoxSizeParametersExtractGautomatch(self):
        # No training mode picking, box size not provided by user
        prot = self._runBoxSizeParamsTestExtrGaut(label='Picking box size related parameters'
                                                                           ' extraction and gautomatch')

        self.assertEqual(prot.boxSizeExtraction, 102, "Estimated box size extraction does not match.")
        self.assertEqual(prot.radiusGautomatch, 180, "Estimated radius for "
                                                                        "gautomatch does not match.")
        self.assertEqual(prot.minIntPartDistanceGautomatch, 216, "Estimated min inter "
                                                                                    "particle distance for gautomatch "
                                                                                    "does not match.")
        self.assertEqual(prot.sigmaDiameterGautomatch, 288, "Estimated sigma diameter for "
                                                                               "gautomatch does not match.")
        self.assertEqual(prot.averageDiameterGautomatch, 362, "Estimated average diameter "
                                                                                 "for gautomatch does not match.")

    def testBoxSizeParametersRelionTopazConsensus(self):
        # No training mode picking, box size not provided by user
        prot = self._runBoxSizeParamsTestRelionTopazConsensus(label='Picking box size related parameters '
                                                                    'relion, topaz and radius consensus')

        self.assertEqual(prot.minLoGFilterRelion, 228, "Estimated minLoGFilterRelion for relion does not match.")
        self.assertEqual(prot.maxLoGFilterRelion, 252, "Estimated maxLoGFilterRelion for relion does not match.")
        self.assertEqual(prot.radiusTopaz, 30, "Estimated radius for topaz does not match.")
        self.assertEqual(prot.numPartPerImgTopaz, 300, "Estimated average diameter numPartPerImgTopaz "
                                                       "for topaz does not match.")
        self.assertEqual(prot.radiusConsensus, 62, "Estimated radius for picking consensus does not match.")

    def _runBoxSizeParamsTestExtrGaut(cls, label):
        protBoxSizeParams = cls.newProtocol(ProtBoxSizeParameters,
                                            inputMicrographs=cls.protImport.outputMicrographs,
                                            boolExtractPartBx=True,
                                            boolGautomatchParams=True)

        protBoxSizeParams.boxSize.setPointer(Pointer(cls.protBoxSizeExtraction, extended="boxSizeExtraction"))
        protBoxSizeParams.setObjLabel(label)
        cls.launchProtocol(protBoxSizeParams)

        return protBoxSizeParams

    def _runBoxSizeParamsTestRelionTopazConsensus(cls, label):
        protBoxSizeParams = cls.newProtocol(ProtBoxSizeParameters,
                                            inputMicrographs=cls.protImport.outputMicrographs,
                                            boolRelionParams=True,
                                            boolTopazParams=True,
                                            boolConsensusParams=True)

        protBoxSizeParams.boxSize.setPointer(Pointer(cls.protBoxSizeExtraction, extended="boxSizeExtraction"))
        protBoxSizeParams.setObjLabel(label)
        cls.launchProtocol(protBoxSizeParams)

        return protBoxSizeParams

