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
from pwem.protocols.protocol_import import ProtImportMicrographs
from pwem.protocols.protocol_mathematical_operator import ProtMathematicalOperator



class TestMathematicalOperator(BaseTest):
    """ Test mathematical operator protocol"""

    @classmethod
    def setData(cls):
        cls.ds = DataSet.getDataSet('relion_tutorial')

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.setData()
        # Run needed protocols
        cls.runImportMicrograph()
        cls.runManualOperator()

    @classmethod
    def runImportMicrograph(cls):

        """ Run an Import micrograph protocol. """
        protImport = cls.newProtocol(
            ProtImportMicrographs,
            samplingRateMode=0,
            filesPath=TestMathematicalOperator.ds.getFile('micrographs/*.mrc'),
            samplingRate=3.54,
            magnification=59000,
            voltage=300,
            sphericalAberration=2)

        cls.launchProtocol(protImport)
        cls.protImport = protImport

    @classmethod
    def runManualOperator(cls):
        protBoxSize = cls.newProtocol(ProtMathematicalOperator,
                                      boolMain=True)

        protBoxSize.expression.set('100 + (3 * 10)')
        protBoxSize.setObjLabel('box size')
        cls.launchProtocol(protBoxSize)
        cls.protBoxSize = protBoxSize

    def testMathematicalOperatorBoxSize(self):
        # Transform box size (px) to box size (A), X1 is a set attribute and X2 an input pointer
        prot = self._runMathematicalOperatorBoxSize(label='Transform box size (px) to box size (A)')
        self.assertEqual(prot.result, 460, "Box size in Angstroms does not match.")

    def testMathematicalOperatorManual(self):
        # Calculate a manual operation
        prot = self._runMathematicalOperatorManual(label='Manual operation')
        self.assertEqual(prot.result, 8, "Manual operation does not match.")


    def _runMathematicalOperatorBoxSize(cls, label):
        protBoxSizeParams = cls.newProtocol(ProtMathematicalOperator,
                                            inputSet1=cls.protImport.outputMicrographs,
                                            input2=cls.protBoxSize.result)

        protBoxSizeParams.attribute1.set('item._samplingRate')
        protBoxSizeParams.input2.set(cls.protBoxSize.result)
        protBoxSizeParams.expression.set('X1 * X2')
        protBoxSizeParams.setObjLabel(label)
        cls.launchProtocol(protBoxSizeParams)

        return protBoxSizeParams

    def _runMathematicalOperatorManual(cls, label):
        protBoxSizeParams = cls.newProtocol(ProtMathematicalOperator,
                                            boolMain=True)

        protBoxSizeParams.expression.set('2 + (3 * 2)')
        protBoxSizeParams.setObjLabel(label)
        cls.launchProtocol(protBoxSizeParams)

        return protBoxSizeParams