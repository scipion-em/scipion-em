# **************************************************************************
# *
# * Authors:    Daniel Marchán (da.marchan@cnb.csic.es)
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
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow.object import Integer, Float, String

class ProtMathematicalOperator(EMProtocol):
    """
    Protocol to make mathematical operations on different inputs
    """

    _label = 'mathematical operator'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _defineParams(self, form):
        """
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """
        form.addSection(label='Input')
        # Manual
        form.addParam('boolMain', params.BooleanParam, default=False,
                      label='Manual Operation',
                      help='Select yes if you want to make a manual operation, e.g: 2 + 2 ')

        # Attribute X1
        group1 = form.addGroup('Input X1', condition="boolMain==%d" % False)
        group1.addParam('bool1', params.BooleanParam, default=True,
                        label='Is X1 a Set attribute?',
                        condition="boolMain==%d" % False,
                        help='Select yes if X1 is a set attribute, e.g: sampling rate')
        group1.addParam('inputSet1', params.PointerParam, pointerClass='EMSet, SetOfImages',
                        condition="bool1==%d and boolMain==%d" % (True, False),
                        label='Set 1',
                        help='Please select the Set from where you want to extract the attribute.')
        group1.addParam('attribute1', params.StringParam, label="Select attribute X1:",
                        condition="bool1==%d and boolMain==%d" % (True, False),
                        help='Use the wizard to list all the sets attributes and select one.')
        group1.addParam('input1', params.IntParam,
                        label='Input X1:',
                        condition="bool1==%d and boolMain==%d" % (False, False),
                        default=0,
                        allowsPointers=True,
                        help='Write the value you want to consider as X1 or'
                             ' if it is a input pointer select from the list of available inputs.')
        # Attribute X2
        group2 = form.addGroup('Input X2', condition="boolMain==%d" % False)
        group2.addParam('bool2', params.BooleanParam, default=False,
                        label='Is X2 a Set attribute?',
                        condition="boolMain==%d" % False,
                        help='Select yes if X2 is a set attribute, e.g: sampling rate')
        group2.addParam('inputSet2', params.PointerParam, pointerClass='EMSet, SetOfImages',
                        condition="bool2==%d and boolMain==%d" % (True, False),
                        label='Set 2',
                        help='Please select the Set from where you want to extract the attribute.')
        group2.addParam('attribute2', params.StringParam, label="Select attribute X2:",
                        condition="bool2==%d and boolMain==%d" % (True, False),
                        help='Use the wizard to list all the sets attributes and select one.')
        group2.addParam('input2', params.IntParam,
                        label='Input X2:',
                        condition="bool2==%d and boolMain==%d" % (False, False),
                        default=0,
                        allowsPointers=True,
                        help='Write the value you want to consider as X2 or'
                             ' if it is a input pointer select from the list of available inputs.')
        group3 = form.addGroup('Result')
        group3.addParam('typeResult', params.EnumParam, default=0,
                        choices=['Int', 'Float', 'String'],
                        label='Result type:',
                        help='Choose the variable type you want your result to be, e.g: Integer')
        group3.addParam('expression', params.StringParam, label="Result =",
                        important=True,
                        help='Write the mathematical expression you want to calculate, e.g: (X1 + X2) * X1 * 0.9 \n'
                             'Please note that X1 and X2 would be replaced by the value you selected for each case.'
                             ' It is important to use upper cases when writing the expression. ')

        form.addParallelSection(threads=1, mpi=1)

    def _insertAllSteps(self):
        self._checkNewInput()

    def createOutput(self, modifiedSet):
        pass

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        if hasattr(self, 'attr1') and hasattr(self, 'attr2'):
            return None

        if self.bool1:
            str1 = self.attribute1.get()
            str1 = str1[str1.find('_'):]
            self.attr1 = getattr(self.inputSet1.get(), str1)
        else:
            self.attr1 = self.input1.get()

        if self.bool2:
            str2 = self.attribute2.get()
            str2 = str2[str2.find('_'):]
            self.attr2 = getattr(self.inputSet2.get(), str2)
        else:
            self.attr2 = self.input2.get()

        fDeps = self._insertNewOperationStep(self.attr1, self.attr2)
        self.updateSteps()

    def _insertNewOperationStep(self, attr1, attr2):
        deps = []
        stepId = self._insertFunctionStep('formulaStep', attr1, attr2, prerequisites=[])
        deps.append(stepId)
        return deps

    def _checkNewOutput(self):
        if hasattr(self, 'result'):
            self.createResultOutput(self.result)

    def formulaStep(self, attr1, attr2):
        """
        Applies the formula to each of the attributes or the numeric formula provided.
        Complex python code could be run separating lines with ;
        To use numpy you could do:
        import numpy;
        item._resolution.set(numpy.random.randint(10))
        """
        str1 = str(attr1)
        str2 = str(attr2)
        strInit = self.expression.get()
        strInit = strInit.replace('X1', str1)
        self.formula = strInit.replace('X2', str2)
        self.result = eval(self.formula)

        print('To calculate: %s' % self.formula)
        print('Result = %s' % self.result)

    def createResultOutput(self, result):
        """ The output can be an Integer, Float or String. Other protocols can use it in those
            Params if it has set allowsPointer=True
        """
        if hasattr(self, "result"):
            tpR = self.typeResult.get()
            if tpR == 0:
                resultNum = Integer(result)
            elif tpR == 1:
                resultNum = Float(result)
            elif tpR == 2:
                resultNum = String(result)

            self._defineOutputs(result=resultNum)

    def _summary(self):
        summary = []

        summary.append('To calculate: %s' % self.expression.get())
        if hasattr(self, "formula"):
            summary.append('To calculate: %s' % self.formula)
        if hasattr(self, "result"):
            summary.append('Result = %s' % self.result)

        return summary

    def _validate(self):
        errors = []
        if not self.boolMain.get():
            if self.bool1.get():
                if self.attribute1.get() == None:
                    errors.append('Input X1 is None please select another value.')
            else:
                if self.input1.get() == None:
                    errors.append('Input X1 is None please select another value.')

            if self.bool2.get():
                if self.attribute2.get() == None:
                    errors.append('Input X2 is None please select another value.')
            else:
                if self.input2.get() == None:
                    errors.append('Input X2 is None please select another value.')

        return errors
