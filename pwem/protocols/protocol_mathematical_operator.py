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
    PARAM1 = 'X1'
    PARAM2 = 'X2'
    _label = 'mathematical operator'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _defineParams(self, form):
        """attr2
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """
        form.addSection(label='Input')

        # Attribute X1
        form.addParam('input1', params.IntParam,
                      label='Input %s:' % self.PARAM1,
                      default=3,
                      allowsPointers=True,
                      help='Write the value you want to consider as %s or'
                             ' if it is a input pointer select from the list of available inputs.' % self.PARAM1)

        # Attribute X2
        form.addParam('input2', params.IntParam,
                      label='Input %s:' % self.PARAM2,
                      default=2,
                      allowsPointers=True,
                      help='Write the value you want to consider as %s or'
                             ' if it is a input pointer select from the list of available inputs.' % self.PARAM2)

        # Result
        form.addParam('expression', params.StringParam, label="Formula: ",
                      important=True,
                      default = "%s + (%s*3)" % (self.PARAM1, self.PARAM2),
                      help='Write the mathematical formula you want to calculate, e.g: (%s + %s) * %s * 0.9 \n'
                             'Please note that %s and %s would be replaced by the value you selected for each case.'
                             ' It is important to use upper cases when writing the expression.' % (self.PARAM1, self.PARAM2, self.PARAM1, self.PARAM1, self.PARAM2))

        form.addParam('typeResult', params.EnumParam, default=1,
                      choices=['Int', 'Float', 'String'],
                      label='Result type:',
                      help='Choose the variable type you want your result to be, e.g: Integer',
                      expertLevel=params.LEVEL_ADVANCED)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.computeStep, self.getParam1(), self.getParam2() , prerequisites=[])

    def computeStep(self, attr1, attr2):
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
        strInit = strInit.replace(self.PARAM1, str1)
        self.formula = strInit.replace(self.PARAM2, str2)
        self.info("FORMULA after replacement: %s" % self.formula)
        result = eval(self.formula)
        self.info("RESULT: %s" % result )
        self.createResultOutput(result)

        self.info('To calculate: %s' % self.formula)
        self.info('Result = %s' % result)

    def createResultOutput(self, result):
        """ The output can be an Integer, Float or String. Other protocols can use it in those
            Params if it has set allowsPointer=True
        """

        tpR = self.typeResult.get()
        if tpR == 0:
            resultNum = Integer(result)
        elif tpR == 1:
            resultNum = Float(result)
        elif tpR == 2:
            resultNum = String(result)

        self._defineOutputs(result=resultNum)

    def formulaNeedsParam(self, param):
        """ return True if the formula uses the param passed"""

        return param in self.expression.get()

    def formulaNeedsInput(self):
        """ Returns true if there is need for any input"""
        return self.formulaNeedsParam(self.PARAM1) or self.formulaNeedsParam(self.PARAM2)

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
        if self.formulaNeedsParam(self.PARAM1):
            if self.getParam1() is None:
                errors.append('Input %s should be provided.' % self.PARAM1)

        if self.formulaNeedsParam(self.PARAM2):
            if self.getParam2() is None:
                errors.append('Input %s should be provided.' % self.PARAM2)

        return errors


    def getParam2(self):

        return self.input2.get()

    def getParam1(self):

        return self.input1.get()