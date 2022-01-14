# **************************************************************************
# *
# * Authors:     Pablo Conesa(pconesa@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
from pwem.objects.data import SetOfCoordinates
import numpy as np

class ProtSetFilter(EMProtocol):
    """
    Protocol to filter sets based on its attributes through an expression that
    should evaluate to true or false. Some predefined expresions are stored (i.e.
    distance to center, distance between coordinates)
    """
    _label = 'filter set'
    CHOICE_FORMULA = 0
    CHOICE_DISTANCE_CENTER = 1
    CHOICE_DISTANCE_BETWEEN_COORDS = 2
    CHOICE_RANKED = 3
    CHOICE_LABEL = {CHOICE_FORMULA: 'formula',
                    CHOICE_DISTANCE_CENTER: 'distance to center',
                    CHOICE_DISTANCE_BETWEEN_COORDS:
                        'distance between particles',
                    CHOICE_RANKED: 'ranking'}

    def _defineParams(self, form):
        """
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """
        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam, pointerClass='EMSet',
                      label='Set to filter',
                      help='Set which items will be filtered.')
        form.addParam('operation', params.EnumParam,
                      choices=[self.CHOICE_LABEL[self.CHOICE_FORMULA],
                               self.CHOICE_LABEL[self.CHOICE_DISTANCE_CENTER],
                               self.CHOICE_LABEL[self.CHOICE_DISTANCE_BETWEEN_COORDS],
                               self.CHOICE_LABEL[self.CHOICE_RANKED]
                               ],
                      default = self.CHOICE_FORMULA,
                      label="Select operation",
                      help="Select operation to be performed in the set.\n"
                           " *distance to center* keep coordinates that are farther from"
                           " the center than a given value.\n"
                           " *distance-between-coordinates* if two coordinates \n"
                           " are closer than a given value then keep only one."
                           " *formula* will apply the formula to the chosen attribute"
                           " (i.e, resolution less than 2 )"
                      )

        form.addParam('formula', params.StringParam, label="Passing formula",
                      condition = "operation==%d" % self.CHOICE_FORMULA,
                      help='A python code compatible with eval that should evaluate to True or False, where item represents each of '
                           'the elements of the set. E.g.: item._resolution.get() < 4).'
                           'You could also use modules like "import numpy;  item._resolution .... "')
        form.addParam('distance', params.FloatParam, default=0,
                      condition="operation==%d or operation==%d" %
                                (self.CHOICE_DISTANCE_CENTER,
                                 self.CHOICE_DISTANCE_BETWEEN_COORDS),
                      label="distance (A)",
                      help="distance from coordinates to center or " \
                           "distance between coordinates ")
        form.addParam('keepFirst', params.BooleanParam, default=True,
                      condition="operation==%d" %
                                self.CHOICE_DISTANCE_BETWEEN_COORDS,
                      label="keep first coordinate",
                      help="If 2 or more coordinates are closer than distance"
                           "keep the first one or delete all"
                      )

        form.addParam('threshold', params.StringParam, label="Threshold: ",
                      condition="operation==%d" % (self.CHOICE_RANKED),
                      help='Number/proportion of items to keep:\n\tNumber: n>=1 \n\tProportion: 0<n<1\n\tPercentage: n%\n\n'
                           'Higher/lower values of the attribute: \n\tHigher: positive number\n\tLower: negative number\n\n'
                           'e.g: "-10%" == "-0.1" == 10% of the items with lower values\n'
                           'e.g: "5" == 5 items with higher values')
        form.addParam('rankingField', params.StringParam, label="Ranking field: ",
                      condition="operation==%d" % self.CHOICE_RANKED,
                      help='Attribute to sort the set by.')

    def _insertAllSteps(self):
        operation = self.operation.get()
        if operation == self.CHOICE_FORMULA:
            self._insertFunctionStep(self.formulaStep)
        elif operation == self.CHOICE_DISTANCE_CENTER:
            self._insertFunctionStep(self.distanceCenterStep)
        elif operation == self.CHOICE_DISTANCE_BETWEEN_COORDS:
            self._insertFunctionStep(self.distanceBetweenCoorStep)
        elif operation == self.CHOICE_RANKED:
            self._insertFunctionStep(self.rankingStep)

    def createOutput(self, modifiedSet):
        # TODO: copyInfo does not copy the set of micrographs
        # associate to the setOfCoordinates
        if hasattr(modifiedSet, "setMicrographs"):
            modifiedSet.setMicrographs(self.inputSet.get().getMicrographs())
        outputArgs = {self.inputSet.getExtended(): modifiedSet}
        self._defineOutputs(**outputArgs)

    def distanceCenterStep(self):
        " compute distance of center to coordinates. Filter if it is two close "
        inputSet = self.inputSet.get()
        modifiedSet = inputSet.create(self._getExtraPath())
        modifiedSet.copyInfo(inputSet)
        mic = inputSet.getMicrographs().getFirstItem()
        micXcenter = mic.getDim()[0] // 2
        micYcenter = mic.getDim()[0] // 2
        sampling = mic.getSamplingRate()
        distance = self.distance.get() / sampling
        distance2 = distance * distance
        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            x = item.getX()
            y = item.getY()
            diffX =  (x-micXcenter)
            diffY =  (y-micYcenter)
            if diffX * diffX + diffY * diffY > distance2:
                modifiedSet.append(item)

        self.createOutput(modifiedSet)

    def distanceBetweenCoorStep(self):
        """ filter by distance between coordinates. If they are too
        close do not keep them"""
        inputSet = self.inputSet.get()
        keepFirstNot = not self.keepFirst.get()
        modifiedSet = inputSet.create(self._getExtraPath())
        modifiedSet.copyInfo(inputSet)
        # set of micrographs is not copied in copyInfo
        mic = inputSet.getMicrographs().getFirstItem()
        sampling = mic.getSamplingRate()
        distance = self.distance.get() / sampling
        distance2 = distance * distance

        # duplicate setOfCoordinates in memory
        coordList = []
        oldMicId=-1
        # If coordenates are sorted by micId then
        # we can process much faster
        useBreak = True
        for item in inputSet.iterItems():
            micId = item.getMicId()
            coordList.append([item.getX(), item.getY(), micId])
            if oldMicId > micId:
                useBreak = False
            oldMicId = micId
        keepList = np.full((len(inputSet)), True, dtype=bool)
        # search for close coordinates
        # this one remove both particles
        for start1, coord1 in enumerate(coordList):
            for start2, coord2 in enumerate(coordList[start1+1:]):
                if (coord1[2] != coord2[2]):  # compare micIds
                    if useBreak:
                        break
                    else:
                        continue
                xd = coord1[0] - coord2[0]
                yd = coord1[1] - coord2[1]
                if (xd*xd + yd*yd) < distance2:
                    if keepFirstNot:
                        keepList[start1] = False
                    keepList[start2+start1+1] = False

        for condition, coord in zip(keepList, inputSet):
            if condition:
                modifiedSet.append(item)

        self.createOutput(modifiedSet)


    def formulaStep(self):
        """
        Goes through all items in the input set and applies the formula to each of them using exec.
        Complex python code could be run separating lines with ;  To use numpy you could do
        import numpy; item._resolution.set(numpy.random.randint(10))
        If result is True, item will be transferred to the output set
        """
        inputSet = self.inputSet.get()
        modifiedSet = inputSet.create(self._getExtraPath())
        modifiedSet.copyInfo(inputSet)

        for sourceItem in inputSet.iterItems():
            item = sourceItem.clone()
            exec("item.setEnabled(%s)"% self.formula.get())
            if item.isEnabled():
                modifiedSet.append(item)
        # TODO: copyInfo does not copy the set of micrographs
        # associate to the setOfCoordinates
        if isinstance(modifiedSet, SetOfCoordinates):
            modifiedSet.setMicrographs(inputSet.getMicrographs())
        self.createOutput(modifiedSet)

    def rankingStep(self):
        """
        Goes through all items in the input set and takes the number/proportion of items with a higher/lower value
        of the chosen attribute
        """
        inputSet = self.inputSet.get()
        finalNumber, direc = self.parseTopRankParam()
        modifiedSet = self.getTopRankItems(self.rankingField.get(), inputSet, finalNumber, direc)
        self.createOutput(modifiedSet)

    def parseTopRankParam(self):
        direc = 'DESC'
        inputSet, threshold = self.inputSet.get(), self.threshold.get().strip()
        # Filed with % at the end
        if threshold.endswith('%'):
            perc = float(threshold[:-1])
            finalNumber = round(perc * len(inputSet) / 100)

        # percentage specified in decimal format.
        elif -1 < float(threshold) < 1:
            prop = float(threshold)
            finalNumber = round(prop * len(inputSet))

        # Else integers (positive or negative).
        else :
            finalNumber = int(float(threshold))

        # If negative
        if finalNumber < 0:
            finalNumber = abs(finalNumber)
            direc = 'ASC'

        return finalNumber, direc

    def getTopRankItems(self, attribute, iSet, finalNumber, direc='ASC'):
        modifiedSet = iSet.createCopy(self._getExtraPath(), copyInfo=True)
        for item in iSet.iterItems(orderBy=attribute, direction=direc, limit=finalNumber):
            modifiedSet.append(item.clone())
        return modifiedSet


    def _validate(self):
        errors = []
        inputSet = self.inputSet.get()
        operation = self.operation.get()
        if operation == self.CHOICE_DISTANCE_CENTER or operation == self.CHOICE_DISTANCE_BETWEEN_COORDS:
            if not isinstance(inputSet, SetOfCoordinates):
                errors.append("The input data set is not a set of coordinates")
        elif operation == self.CHOICE_RANKED:
            param = self.threshold.get().strip()
            perc = False
            try:
                if param.endswith('%'):
                    param = param[:-1]
                    perc=True
                param = float(param)
                if perc:
                    if abs(param) > 100:
                        errors.append("Percentage cannot be higher than 100")
            except:
                errors.append("The filter value must be a number or a percentage")
        return errors
