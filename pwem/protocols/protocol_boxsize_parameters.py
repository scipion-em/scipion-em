# **************************************************************************
# *
# * Authors:     Daniel Marchán (da.marchan@cnb.csic.es)
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
import enum

import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow.object import Integer
from collections import OrderedDict

RADIUS_GAUTOMATCH ='radiusGautomatch'
EXTRACTION = 'boxSizeExtraction'

class ProtBoxSizeParameters(EMProtocol):
    """
    Protocol to make mathematical operations on particle picking boxsize
    """

    _label = 'box size related parameters'
    _possibleOutputs = {RADIUS_GAUTOMATCH: Integer,
                        EXTRACTION: Integer}
    outputsToDefine = {}

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _defineParams(self, form):
        """
        Defines the parameters the protocol form will show and its behaviour
        :param form:
        """
        form.addSection(label='Input')

        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      important=True,
                      label="Input micrographs",
                      help='Select the SetOfMicrographs from where we extract the sampling rate')

        form.addParam('boxSize', params.IntParam,
                      label='Particle box size (px)',
                      allowsPointers=True,
                      important=True,
                      help='This is size of the boxed particles (in pixels)')

        # Extraction box size
        form.addParam('boolExtractPartBx', params.BooleanParam, default=False,
                      label='Extract particles boxsize?',
                      help='Select yes if you want to generate these parameters.')
        groupEx = form.addGroup('Extract particles boxsize',
                                condition="boolExtractPartBx==%d" % True)
        groupEx.addParam('factorExtractPartBx', params.FloatParam, default=1.5,
                         label='Factor to multiply the box size (px)',
                         condition="boolExtractPartBx==%d" % True,
                         help='Extraction box size (px) = picking box size (px) * *factor*')

        # Gautomatch picking parameters
        form.addParam('boolGautomatchParams', params.BooleanParam, default=False,
                      label='Calculate gautomatch picking parameters?',
                      help='Select yes if you want to generate these parameters.')
        groupGaut = form.addGroup('Gautomatch picking parameters',
                                  condition="boolGautomatchParams==%d" % True)
        groupGaut.addParam('factorGautRadius', params.FloatParam, default=0.75,
                        label='Factor to obtain the particle radius (A)',
                        condition="boolGautomatchParams==%d" % True,
                        help='Particle radius in Angstrom. Default will be equal to 75% of reference size (box size).'
                             '\nParticle radius (A) = picking box size (px) * sampling_rate (A/px) * *factor*')
        groupGaut.addParam('factorGautMinInterPartDist', params.FloatParam, default=0.9,
                           label='Factor to obtain the min inter-particle distance (A)',
                           condition="boolGautomatchParams==%d" % True,
                           help='Minimum distance between particles in Angstrom\n '
                                'Use value of 0.9~1.1X diameter; '
                                'can be 0.3~0.5X for filament-like particle.'
                                '\nMin inter-particle distnace (A) = picking box size (px) * sampling_rate (A/px)'
                                ' * *factor*')
        groupGaut.addParam('factorGautSigmaDiameter', params.FloatParam, default=1.2,
                           label='Factor to obtain the local sigma diameter (A)',
                           condition="boolGautomatchParams==%d" % True,
                           help='Diameter for estimation of local sigma, '
                                'in Angstrom.\n'
                                'Usually this diameter could be 0.5-2x of your '
                                'particle diameter according to several factors. '
                                'When using bigger values, normally you should '
                                'decrease *Local sigma cut-off*. For smaller and '
                                'sharper high density contamination/ice/metal '
                                'particles you could use a smaller diameter and '
                                'larger *Local sigma cut-off*.'
                                '\nLocal sigma diameter (A) = picking box size (px) * sampling_rate (A/px) * *factor*')
        groupGaut.addParam('factorGautAvgDiameter', params.FloatParam, default=1.5,
                           label='Factor to obtain the local average diameter (A)',
                           condition="boolGautomatchParams==%d" % True,
                           help='Diameter for estimation of local average, in '
                                'Angstrom. 1.5~2.0X particle diameter suggested.\n'
                                'However, if you have sharp/small ice or any '
                                'dark/bright dots, using a smaller value will be '
                                'much better to get rid of these areas.'
                                '\nLocal average diameter (A) = picking box size (px) * sampling_rate (A/px) '
                                '* *factor*')

        # Relion box size
        form.addParam('boolRelionParams', params.BooleanParam, default=False,
                      label='Calculate relion picking parameters?',
                      help='Select yes if you want to generate these parameters.')
        groupRelion = form.addGroup('Relion picking parameters',
                                    condition="boolRelionParams==%d" % True)
        groupRelion.addParam('factorMinLoGFilter', params.FloatParam, default=0.95,
                             label='Factor to obtain the Min diameter for LoG filter (A)',
                             condition="boolRelionParams==%d" % True,
                             help='This should correspond to the smallest size '
                                  'of your particles projections in Ångstroms.'
                                  '\nMin diameter for LoG filter (A) = picking box size (px) '
                                  '* sampling_rate (A/px) * *factor*')
        groupRelion.addParam('factorMaxLoGFilter', params.FloatParam, default=1.05,
                             label='Factor to obtain the Max diameter for LoG filter (A)',
                             condition="boolRelionParams==%d" % True,
                             help='This should correspond to the largest size '
                                  'of your particles projections in Ångstroms.'
                                  '\nMax diameter for LoG filter (A) = picking box size (px) '
                                  '* sampling_rate (A/px) * *factor*')

        # Topaz picking parameters
        form.addParam('boolTopazParams', params.BooleanParam, default=False,
                      label='Calculate topaz picking parameters?',
                      help='Select yes if you want to generate these parameters.')
        groupGaut = form.addGroup('Topaz picking parameters',
                                  condition="boolTopazParams==%d" % True)
        groupGaut.addParam('factorTopazRadius', params.FloatParam, default=0.45,
                           label='Factor to obtain the particle radius (px)',
                           condition="boolTopazParams==%d" % True,
                           help='Pixel radius around particle centers to consider.'
                                '\nParticle radius (px) = picking box size (px) * *factor*')
        groupGaut.addParam('numPartPerImg', params.IntParam, default=300,
                           label='Number of particles per image',
                           condition="boolTopazParams==%d" % True,
                           help='Expected number of particles per micrograph.\n'
                                ' If -1 it will be estimated for you.')

        # Particle picking consensus parameters
        form.addParam('boolConsensusParams', params.BooleanParam, default=False,
                      label='Calculate picking consensus parameters?',
                      help='Select yes if you want to generate these parameters.')
        groupGaut = form.addGroup('Picking consensus parameters',
                                  condition="boolConsensusParams==%d" % True)
        groupGaut.addParam('factorConsensusRadius', params.FloatParam, default=0.9,
                           label='Factor to obtain the particle radius (px)',
                           condition="boolConsensusParams==%d" % True,
                           help='Pixel radius around particle centers to consider.'
                                '\nParticle radius (px) = picking box size (px) * *factor*')

        form.addParallelSection(threads=1, mpi=1)

    def _insertAllSteps(self):
        self.initParams()
        self._checkNewInput()

    def initParams(self):
        self.dictParams = OrderedDict()
        self.outputDone = False

    def createOutput(self, modifiedSet):
        pass

    def _stepsCheck(self):
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        if hasattr(self, 'samplingRate') and hasattr(self, 'boxSize'):
            return None

        self.boxSize = self.boxSize.get()
        self.samplingRate = self.inputMicrographs.get().getSamplingRate()

        fDeps = self._insertNewOperationsStep(self.boxSize, self.samplingRate)
        self.updateSteps()

    def _insertNewOperationsStep(self, boxSize, samplingRate):
        deps = []
        stepId = self._insertFunctionStep('applyFormulaStep', boxSize, samplingRate, prerequisites=[])
        deps.append(stepId)
        return deps

    def _checkNewOutput(self):
        if self.outputDone:
            self.createResultsOutput(self.dictParams)

    def applyFormulaStep(self, boxSize, samplingRate):
        """
        Applies the formula to each of the parameters selected by the user.
        """

        if self.boolExtractPartBx.get():
            self.calculateParticleExtractionParams(boxSize)

        if self.boolGautomatchParams.get():
            self.dictParams['gautomatch'] = self.calculateGautomatchParams(boxSize, samplingRate)
            print('Gautomatch particle radius: %d (A)' % self.dictParams['gautomatch'][RADIUS_GAUTOMATCH])
            print('Gautomatch min inter particle distance: %d (A)'
                  % self.dictParams['gautomatch']['minIntPartDistanceGautomatch'])
            print('Gautomatch local sigma diameter: %d (A)'
                  % self.dictParams['gautomatch']['sigmaDiameterGautomatch'])
            print('Gautomatch local average diameter: %d (A) \n'
                  % self.dictParams['gautomatch']['averageDiameterGautomatch'])

        if self.boolRelionParams.get():
            self.dictParams['relion'] = self.calculateRelionParams(boxSize, samplingRate)
            print('Relion min diameter for loG filter: %d (A)' % self.dictParams['relion']['minLoGFilterRelion'])
            print('Relion max diameter for loG filter: %d (A) \n' % self.dictParams['relion']['maxLoGFilterRelion'])

        if self.boolTopazParams.get():
            self.dictParams['topaz'] = self.calculateTopazParams(boxSize)
            print('Topaz particle radius: %d (px)' % self.dictParams['topaz']['radiusTopaz'])
            print('Topaz number of particles per image: %d \n' % self.dictParams['topaz']['numPartPerImgTopaz'])

        if self.boolConsensusParams.get():
            self.dictParams['consensus'] = self.calculateConsensusParams(boxSize)
            print('Picking consensus radius: %d (px)' % self.dictParams['consensus']['radiusConsensus'])

        self.outputDone = True


    def registerOutput(self, outputName, value, msg):
        self.outputsToDefine[outputName] = value
        self.info(msg % value)

    def calculateParticleExtractionParams(self, boxSize):

        self.registerOutput(EXTRACTION,
                            Integer(boxSize * self.factorExtractPartBx.get()),
                            'Particle extraction box size: %d (px) \n')


    def calculateGautomatchParams(self, boxSize, samplingRate):
        gautomatchDict = {}
        gautomatchDict[RADIUS_GAUTOMATCH] = Integer(boxSize * samplingRate * self.factorGautRadius.get())
        gautomatchDict['minIntPartDistanceGautomatch'] = Integer(boxSize * samplingRate * self.factorGautMinInterPartDist.get())
        gautomatchDict['sigmaDiameterGautomatch'] = Integer(boxSize * samplingRate * self.factorGautSigmaDiameter.get())
        gautomatchDict['averageDiameterGautomatch'] = Integer(boxSize * samplingRate * self.factorGautAvgDiameter.get())

        return gautomatchDict

    def calculateRelionParams(self, boxSize, samplingRate):
        relionDict = {}
        relionDict['minLoGFilterRelion'] = Integer(boxSize * samplingRate * self.factorMinLoGFilter.get())
        relionDict['maxLoGFilterRelion'] = Integer(boxSize * samplingRate * self.factorMaxLoGFilter.get())

        return relionDict

    def calculateTopazParams(self, boxSize):
        topazDict = {}
        topazDict['radiusTopaz'] = Integer(boxSize * self.factorTopazRadius.get())
        # TODO: the numPartPerImg needs to be estimated
        topazDict['numPartPerImgTopaz'] = Integer(self.numPartPerImg.get())

        return topazDict

    def calculateConsensusParams(self, boxSize):
        consensusDict = {}
        consensusDict['radiusConsensus'] = Integer(boxSize * self.factorConsensusRadius.get())

        return consensusDict


    def createResultsOutput(self, dictParams):
        """ The output can be an Integer, Float or String. Other protocols can use it in those
            Params if it has set allowsPointer=True
        """

        self._defineOutputs(**self.outputsToDefine)

    def _summary(self):
        summary = []

        return summary

    def _validate(self):
        errors = []

        if self.inputMicrographs.get().getSamplingRate() == None:
            errors.append('The input micrographs do not have a sampling rate')

        return errors