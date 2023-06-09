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
from pwem.protocols.protocol_boxsize_checkpoint import ProtBoxSizeCheckpoint


class TestBoxSizeCheckpoint(BaseTest):
    """ Test box size checkpoint protocol"""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)


    def testBoxSizeCheckpoint(self):
        # Calculate a manual operation
        protNoInput = self._runBoxSizeCheckpointNoInputs(label='No inputs BoxSize checkpoint')
        self.assertEqual(protNoInput.boxSize, 300, "Computation without input does not match.")

        # Feed the protocol with boxSize1
        # Calculate an average and surpassing thresholds
        protWithInput = self.newProtocol(ProtBoxSizeCheckpoint,
                                        boxSize2= 400,
                                         boolBoxSize1=False,
                                        boolBoxSizeAvg= True)

        pntToOutput = Pointer(protNoInput, extended="boxSize")
        protWithInput.boxSize1.setPointer(pntToOutput)
        protWithInput.setObjLabel("Surpass threshold and average")
        self.launchProtocol(protWithInput)
        self.assertEqual(protWithInput.boxSize.get(), 350, "Computation without input does not match.")

        # Feed the protocol with boxSize2 input as well
        # Calculate with timer operation
        protWithBothInput = self.newProtocol(ProtBoxSizeCheckpoint,
                                             boolTimer=True,
                                             timeout="2m"
                                             )
        pntToOutput = Pointer(protNoInput, extended="boxSize")
        protWithBothInput.boxSize1.setPointer(pntToOutput)
        pntToOutput = Pointer(protWithInput, extended="boxSize")
        protWithBothInput.boxSize2.setPointer(pntToOutput)
        protWithBothInput.setObjLabel("2 inputs within threshoold and timer")
        self.launchProtocol(protWithBothInput)
        self.assertEqual(protWithBothInput.boxSize.get(), 300 ,"Computation without input does not match.")


    def _runBoxSizeCheckpointNoInputs(cls, label):
        protBoxSizeCheckpoint = cls.newProtocol(ProtBoxSizeCheckpoint,
                                            boxSize1= 300,
                                            boxSize2= 310)

        protBoxSizeCheckpoint.setObjLabel("Within threshold and stay with reference")
        cls.launchProtocol(protBoxSizeCheckpoint)

        return protBoxSizeCheckpoint
