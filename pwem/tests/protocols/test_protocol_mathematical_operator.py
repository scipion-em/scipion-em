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
from pwem.protocols.protocol_mathematical_operator import ProtMathematicalOperator



class TestMathematicalOperator(BaseTest):
    """ Test mathematical operator protocol"""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def testMethods(self):

        prot = self.newProtocol(ProtMathematicalOperator,
                                      boolMain=True)
        self.assertTrue(prot.formulaNeedsParam(ProtMathematicalOperator.PARAM1), "formulaNeedsParam false negative for PARAM1")
        self.assertTrue(prot.formulaNeedsParam(ProtMathematicalOperator.PARAM2), "formulaNeedsParam false negative for PARAM2")

        # validation
        self.assertEqual(0, len(prot.validate()), "validation for default value do not work")

        prot.expression.set("SDSKJSL")
        self.assertFalse(prot.formulaNeedsParam(ProtMathematicalOperator.PARAM1), "formulaNeedsParam false positive for PARAM1")
        self.assertFalse(prot.formulaNeedsParam(ProtMathematicalOperator.PARAM2), "formulaNeedsParam false positive for PARAM2")

        prot.expression.set(ProtMathematicalOperator.PARAM1 + "LOLO")
        self.assertTrue(prot.formulaNeedsParam(ProtMathematicalOperator.PARAM1), "formulaNeedsParam false negative for PARAM1")
        self.assertFalse(prot.formulaNeedsParam(ProtMathematicalOperator.PARAM2), "formulaNeedsParam false positive for PARAM2")


    def testMathematicalOperatorManual(self):

        # Calculate a manual operation
        protNoInput = self.newProtocol(ProtMathematicalOperator)

        protNoInput.expression.set('2 + (3 * 2)')
        protNoInput.setObjLabel("No inputs")
        self.launchProtocol(protNoInput)

        self.assertEqual(protNoInput.result, 8, "Computation without input does not match.")

        # Feed the protocol with inputs
        # Calculate a manual operation
        protWithInput = self.newProtocol(ProtMathematicalOperator)
        pntToOutput = Pointer(protNoInput, extended="result")
        protWithInput.input1.setPointer(pntToOutput)
        protWithInput.expression.set('X1 + 1')
        protWithInput.setObjLabel("X1 as input")
        self.launchProtocol(protWithInput)

        self.assertEqual(protWithInput.result.get(), 9, "Computation without input does not match.")

        # Feed the protocol with X2 input as well
        # Calculate a manual operation
        protWithBothInput = self.newProtocol(ProtMathematicalOperator)
        pntToOutput = Pointer(protNoInput, extended="result")
        protWithBothInput.input1.setPointer(pntToOutput)

        pntToOutput = Pointer(protWithInput, extended="result")
        protWithBothInput.input2.setPointer(pntToOutput)

        protWithBothInput.expression.set('X1 - (X1*X2)')
        protWithBothInput.setObjLabel("X1 and X2 as input")
        self.launchProtocol(protWithBothInput)

        self.assertEqual(protWithBothInput.result.get(), -64, "Computation without input does not match.")

