# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import pwem.objects as emobj
from pwem.protocols import ProtCTFMicrographs


class ProtCTFAssign(ProtCTFMicrographs):
    """ This protocol assigns a CTF estimation to a particular
    set of particles producing a new set. """
    _label = 'assign ctf'
    _unionTypes = ['Micrographs',
                   'Particles']
    
    def __init__(self, **kwargs):
        ProtCTFMicrographs.__init__(self, **kwargs)
        # We need to trace the changes of 'inputType' to 
        # dynamically modify the property of pointerClass
        # of the 'inputSets' parameter

        def onChangeInputType():
            pointerClass = 'SetOf' + self.getEnumText('inputType')
            self.inputSetsParam.setPointerClass(pointerClass)
        # Initial update
        onChangeInputType()
        # Now keep track of changes and update
        self.inputType.trace(onChangeInputType)
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputType', params.EnumParam, choices=self._unionTypes,
                      default=0,  # Micrographs
                      label='Input type:',
                      help='Select the type of objects that you want to assign'
                           ' the CTF.')
        self.inputSetsParam = form.addParam('inputSet', params.PointerParam,
                                            pointerClass='EMSet',
                                            label='Input set',
                                            help='Select the images '
                                                 '(micrographs or particles) '
                                                 'that you want to update the '
                                                 'CTF parameters.')        
        form.addParam('inputCTF', params.PointerParam,
                      pointerClass='SetOfCTF, SetOfParticles',
                      label="Input CTF",
                      help='Select the CTFs that will be used to update '
                           'particles. It can be another set of particles')
        form.addParallelSection(threads=0, mpi=0) 
    
    # --------------------------- INSERT steps functions ----------------------
                                
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')
    
    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self):
        inputSet = self.inputSet.get()
        inputCTF = self.inputCTF.get()
        
        if isinstance(inputSet, emobj.SetOfParticles):
            self._particlesOutputStep(inputSet, inputCTF)
        else:
            self._micrographsOutputStep(inputSet, inputCTF)
            
    def _particlesOutputStep(self, inputSet, inputCTF):
        outputParts = self._createSetOfParticles()
        outputParts.copyInfo(inputSet)
        outputParts.setHasCTF(True)

        if isinstance(inputCTF, emobj.SetOfParticles):
            outputParts.copyItems(inputSet,
                                  itemDataIterator=inputCTF.iterItems(),
                                  updateItemCallback=self._updateItem)
        else:
            ctfDict = {}
            firstCoord = inputSet.getFirstItem().getCoordinate()
            hasMicName = firstCoord.getMicName() is not None

            for ctf in inputCTF:
                if hasMicName:
                    ctfName = ctf.getMicrograph().getMicName()
                else:
                    ctfName = ctf.getMicrograph().getObjId()
    #             print("ctf: ", ctf.printAll(), ctfName)
                ctfDict[ctfName] = ctf.clone()

            missingSet = set()  # Report missing micrographs only once

            for particle in inputSet:
                if particle.hasCoordinate():
                    coord = particle.getCoordinate()
                    micKey = coord.getMicName() if hasMicName else particle.getMicId()
                else:
                    micKey = particle.getMicId()

                if micKey not in missingSet:
                    ctf = ctfDict.get(micKey, None)

                    if ctf is None:
                        self.warning("Discarding particles from micrograph with"
                                     " micName: %s, CTF not found. " % micKey)
                        missingSet.add(micKey)
                    else:
                        newParticle = particle.clone()
                        newParticle.setCTF(ctf)
                        outputParts.append(newParticle)
        
        self._defineOutputs(outputParticles=outputParts)
        self._defineSourceRelation(self.inputSet, outputParts)
        self._defineSourceRelation(self.inputCTF, outputParts)

    def __findCTF(self, inputSet, outputSet, ctfDict, keyFunc):
        for mic in inputSet:
            micKey = keyFunc(mic)
            ctf = ctfDict.get(micKey, None)
            if ctf is None:
                self.warning("Discarding micrographs with micName: %s, CTF not found. " % micKey)
            else:
                newMic = mic.clone()
                outputSet.append(newMic)
        
    def _micrographsOutputStep(self, inputSet, inputCTF):
        outputMics = self._createSetOfMicrographs()
        outputMics.copyInfo(inputSet)
        ctfDict = {}
        
        for ctf in inputCTF:
            mic = ctf.getMicrograph()
            ctfDict[mic.getMicName()] = ctf
            ctfDict[mic.getObjId()] = ctf
        
        # Try first to find the ctf by the micrograph micName
        self.__findCTF(inputSet, outputMics, ctfDict, lambda mic: mic.getMicName())
        
        # Now, if no ctf was found, try using the micId
        if outputMics.getSize() == 0:
            self.warning("No CTF found using micName, now trying with micId")
            self.__findCTF(inputSet, outputMics, ctfDict, lambda mic: mic.getObjId())            
        
        self._defineOutputs(outputMicrographs=outputMics)
        self._defineSourceRelation(self.inputSet, outputMics)
        self._defineCtfRelation(outputMics, self.inputCTF)
    
    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary    
    
    def _methods(self):
        return []
    
    def _validate(self):
        """ The function of this hook is to add some validation before the
        protocol is launched to be executed. It should return a list of errors.
        If the list is empty the protocol can be executed.
        """
        errors = []
        # Add some errors if input is not valid
        inputSet = self.inputSet.get()
        if isinstance(inputSet, emobj.SetOfParticles):
            part = inputSet.getFirstItem()
            if not part.hasMicId():
                errors.append("The input particles doesn't have any micrograph"
                              " assigned.")

        # same micrographs in both CTF??
        return errors
    
    # --------------------------- UTILS functions -----------------------------
    def _stepsCheck(self):
        pass

    def _updateItem(self, item, row):
        item.setCTF(row.getCTF())
