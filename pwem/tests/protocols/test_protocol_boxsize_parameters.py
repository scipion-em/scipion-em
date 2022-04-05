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



from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.plugin import Domain
from pwem.protocols.protocol_import import ProtImportMicrographs
from pwem.protocols.protocol_boxsize_parameters import ProtBoxSizeParameters

SphireProtCRYOLOPicking = Domain.importFromPlugin('sphire.protocols', 'SphireProtCRYOLOPicking', doRaise=True)
INPUT_MODEL_GENERAL = Domain.importFromPlugin('sphire.constants', 'INPUT_MODEL_GENERAL', doRaise=True)


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
        cls.runCryoloPicking()

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
    def runCryoloPicking(cls):
        """ Run cryolo picking protocol. """
        protCryolo = cls.newProtocol(SphireProtCRYOLOPicking,
                                     # inputMicrographs=cls.protImport.outputMicrographs,
                                     boxSize=0,
                                     input_size=750,
                                     boxSizeFactor=1.05,
                                     streamingBatchSize=10)
        protCryolo.inputMicrographs.set(cls.protImport.outputMicrographs)
        cls.launchProtocol(protCryolo)
        cls.protCryolo = protCryolo

    def testBoxSizeParametersExtractGautomatch(self):
        # No training mode picking, box size not provided by user
        prot = self._runBoxSizeParamsTestExtrGaut(label='Picking box size related parameters extraction and gautomatch')

        self.assertEqual(prot.boxSizeExtraction, 69, "Estimated box size extraction does not match.")
        self.assertEqual(prot.radiusGautomatch, 122, "Estimated radius for gautomatch does not match.")
        self.assertEqual(prot.minIntPartDistanceGautomatch, 146, "Estimated min inter particle distance"
                                                                 " for gautomatch does not match.")
        self.assertEqual(prot.sigmaDiameterGautomatch, 195, "Estimated sigma diameter for gautomatch does not match.")
        self.assertEqual(prot.averageDiameterGautomatch, 244, "Estimated average diameter "
                                                              "for gautomatch does not match.")

    def testBoxSizeParametersRelionTopazConsensus(self):
        # No training mode picking, box size not provided by user
        prot = self._runBoxSizeParamsTestRelionTopazConsensus(label='Picking box size related parameters relion and topaz')

        self.assertEqual(prot.minLoGFilterRelion, 154, "Estimated minLoGFilterRelion for relion does not match.")
        self.assertEqual(prot.maxLoGFilterRelion, 170, "Estimated maxLoGFilterRelion for relion does not match.")
        self.assertEqual(prot.radiusTopaz, 20, "Estimated radius for topaz does not match.")
        self.assertEqual(prot.numPartPerImgTopaz, 300, "Estimated average diameter numPartPerImgTopaz "
                                                       "for topaz does not match.")
        self.assertEqual(prot.radiusConsensus, 41, "Estimated radius for picking consensus does not match.")

    def _runBoxSizeParamsTestExtrGaut(cls, label):
        protBoxSizeParams = cls.newProtocol(ProtBoxSizeParameters,
                                            inputMicrographs=cls.protImport.outputMicrographs,
                                            boxSize=cls.protCryolo.boxsize,
                                            boolExtractPartBx=True,
                                            boolGautomatchParams=True)

        protBoxSizeParams.boxSize.set(cls.protCryolo.boxsize)
        protBoxSizeParams.setObjLabel(label)
        cls.launchProtocol(protBoxSizeParams)

        return protBoxSizeParams

    def _runBoxSizeParamsTestRelionTopazConsensus(cls, label):
        protBoxSizeParams = cls.newProtocol(ProtBoxSizeParameters,
                                            inputMicrographs=cls.protImport.outputMicrographs,
                                            boxSize=cls.protCryolo.boxsize,
                                            boolRelionParams=True,
                                            boolTopazParams=True,
                                            boolConsensusParams=True)

        protBoxSizeParams.boxSize.set(cls.protCryolo.boxsize)
        protBoxSizeParams.setObjLabel(label)
        cls.launchProtocol(protBoxSizeParams)

        return protBoxSizeParams

