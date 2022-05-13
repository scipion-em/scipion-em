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
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow.object import Integer

BOXSIZE_EVEN = 'boxSizeEven'
EXTRACTION = 'boxSizeExtraction'
RADIUS_GAUTOMATCH = 'radiusGautomatch'
MIN_INTPAR_DIST_GAUTOMATCH= 'minIntPartDistanceGautomatch'
SIGMA_DIAM_GAUTOMATCH = 'sigmaDiameterGautomatch'
AVG_DIAM_GAUTOMATCH = 'averageDiameterGautomatch'
MIN_LOG_FILTER_RELION = 'minLoGFilterRelion'
MAX_LOG_FILTER_RELION = 'maxLoGFilterRelion'
RADIUS_TOPAZ = 'radiusTopaz'
NUM_PART_IMG_TOPAZ = 'numPartPerImgTopaz'
RADIUS_CONSENSUS = 'radiusConsensus'

class ProtBoxSizeParameters(EMProtocol):
    """
    Protocol to make mathematical operations on particle picking boxsize.
    For sanity check all the generated outputs are even numbers.
    """

    _label = 'box size related parameters'
    _possibleOutputs = {BOXSIZE_EVEN: Integer,
                        EXTRACTION: Integer,
                        RADIUS_GAUTOMATCH: Integer,
                        MIN_INTPAR_DIST_GAUTOMATCH: Integer,
                        SIGMA_DIAM_GAUTOMATCH: Integer,
                        AVG_DIAM_GAUTOMATCH: Integer,
                        MIN_LOG_FILTER_RELION: Integer,
                        MAX_LOG_FILTER_RELION: Integer,
                        RADIUS_TOPAZ: Integer,
                        NUM_PART_IMG_TOPAZ: Integer,
                        RADIUS_CONSENSUS: Integer}

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
                      help='This is size of the boxed particles (in pixels).\n'
                           'For sanity check if it is not, it will be transform to an even number.')

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
            self.createResultsOutput()

    def applyFormulaStep(self, boxSize, samplingRate):
        """
        Applies the formula to each of the parameters selected by the user.
        """

        self.registerEvenBoxSize(boxSize)

        if self.boolExtractPartBx.get():
            self.calculateParticleExtractionParams(boxSize)

        if self.boolGautomatchParams.get():
            self.calculateGautomatchParams(boxSize, samplingRate)

        if self.boolRelionParams.get():
            self.calculateRelionParams(boxSize, samplingRate)

        if self.boolTopazParams.get():
            self.calculateTopazParams(boxSize)

        if self.boolConsensusParams.get():
            self.calculateConsensusParams(boxSize)

        self.outputDone = True

    def registerOutput(self, outputName, value):
        self.outputsToDefine[outputName] = value

    def registerEvenBoxSize(self, boxSize):
        boxSize = transform2EvenNumber(boxSize)
        self.registerOutput(BOXSIZE_EVEN, Integer(boxSize))

    def calculateParticleExtractionParams(self, boxSize):
        self.registerOutput(EXTRACTION,
                            Integer(transform2EvenNumber(boxSize * self.factorExtractPartBx.get())))

    def calculateGautomatchParams(self, boxSize, samplingRate):
        self.registerOutput(RADIUS_GAUTOMATCH,
                            Integer(transform2EvenNumber(boxSize * samplingRate * self.factorGautRadius.get())))
        self.registerOutput(MIN_INTPAR_DIST_GAUTOMATCH,
                            Integer(transform2EvenNumber(boxSize * samplingRate * self.factorGautMinInterPartDist.get())))
        self.registerOutput(SIGMA_DIAM_GAUTOMATCH,
                            Integer(transform2EvenNumber(boxSize * samplingRate * self.factorGautSigmaDiameter.get())))
        self.registerOutput(AVG_DIAM_GAUTOMATCH,
                            Integer(transform2EvenNumber(boxSize * samplingRate * self.factorGautAvgDiameter.get())))

    def calculateRelionParams(self, boxSize, samplingRate):
        self.registerOutput(MIN_LOG_FILTER_RELION,
                            Integer(transform2EvenNumber(boxSize * samplingRate * self.factorMinLoGFilter.get())))
        self.registerOutput(MAX_LOG_FILTER_RELION,
                            Integer(transform2EvenNumber(boxSize * samplingRate * self.factorMaxLoGFilter.get())))

    def calculateTopazParams(self, boxSize):
        self.registerOutput(RADIUS_TOPAZ,
                            Integer(transform2EvenNumber(boxSize * self.factorTopazRadius.get())))
        # TODO: the numPartPerImg needs to be estimated
        self.registerOutput(NUM_PART_IMG_TOPAZ,
                            Integer(transform2EvenNumber(self.numPartPerImg.get())))

    def calculateConsensusParams(self, boxSize):
        self.registerOutput(RADIUS_CONSENSUS,
                            Integer(transform2EvenNumber(boxSize * self.factorConsensusRadius.get())))

    def createResultsOutput(self):
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


def transform2EvenNumber(var):
    if var % 2 != 0:
        var = round(var / 2) * 2

    return var