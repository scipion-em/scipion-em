# **************************************************************************
# *
# * Authors:     R. Marabini (roberto@cnb.csic.es)
# *              Tomas Majtner (tmajtner@cnb.csic.es)   -- added particles
# *              Amaya Jimenez (ajimenez@cnb.csic.es) -- Mic for SetOfMics
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

from collections import OrderedDict
from os.path import basename, splitext
import time
import random

import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STEPS_PARALLEL

from pwem.protocols import EMProtocol
import pwem.objects as emobj
from pwem import emlib


SET_OF_MOVIES = 0
SET_OF_MICROGRAPHS = 1
SET_OF_RANDOM_MICROGRAPHS = 2
SET_OF_PARTICLES = 3
SET_OF_COORDINATES = 4


class ProtCreateStreamData(EMProtocol):
    """ create  setofXXXX in streaming mode.
        micrograph -> read a micrograph in memory and writes it nDim times
        movie      -> read a movie in memory and writes it nDim times
        randomMicrographs -> creates a micrograph with random values
        and applies a random CTF
        particles  -> read nDim particles in memory and writes it in streaming
    """
    _label = "create stream data"
    _singleImageFn = "singleImage.xmp"
    _magnification = 500000
    _voltage = 200
    _sphericalAberration = 2.0
    _amplitudeContrast = 0.07
    object = None
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.dictObj = OrderedDict()

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('setof', params.EnumParam, default=0,
                      choices=['Movies', 'Micrographs', 'RandomMicrographs',
                               'Particles', 'Coordinates'],
                      label='set Of',
                      help='create set of')
        form.addParam('inputMovies', params.PointerParam, pointerClass='SetOfMovies',
                      condition="setof==%d" % SET_OF_MOVIES,
                      label="movie",
                      help='This movie will be copied "number of items" times')
        form.addParam('inputMics', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      condition="setof==%d" % SET_OF_MICROGRAPHS,
                      label="micrograph",
                      help='This micrograph will be copied "number of items"'
                           ' times')
        form.addParam('xDim', params.IntParam, default=1024,
                      condition="setof==%d" % SET_OF_RANDOM_MICROGRAPHS,
                      label="xdim",
                      help="X dim ")
        form.addParam('yDim', params.IntParam, default=1024,
                      condition="setof==%d" % SET_OF_RANDOM_MICROGRAPHS,
                      label="ydim",
                      help="Y dim ")
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      condition="setof==%d" % SET_OF_PARTICLES,
                      label="SetOfParticles",
                      help='These particles will be written in streaming')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates',
                      condition="setof==%d" % SET_OF_COORDINATES,
                      label="SetOfCoordinates",
                      help='These Coordinates will be written in streaming')
        form.addParam('groups', params.IntParam, default=50,
                      condition="setof==%d" % SET_OF_PARTICLES,
                      label="groups",
                      help='How many items will be created every iteration')
        form.addParam('nDim', params.IntParam, default=10,
                      label="number of items",
                      help="setofXX size")
        form.addParam('samplingRate', params.FloatParam, default=4,
                      condition="setof!=%d and setof!=%d" % (
                          SET_OF_MICROGRAPHS, SET_OF_MOVIES),
                      label="samplingRate",
                      help="Sampling rate")
        form.addParam('creationInterval', params.IntParam, default=60,
                      label="Create Object each (sec)",
                      pointerClass='EMProtocol',
                      help="create one object each creationInterval seconds")
        form.addParam('extraRandomInterval', params.IntParam, default=0,
                      label='Extra random seconds interval', expertLevel=params.LEVEL_ADVANCED,
                      help='Each object will be generated in a random time uniformly picked from the interval defined by'
                           '[baseInterval, baseInterval+extraInterval]')
        form.addParam('delay', params.IntParam, default=0,
                      label="delay (sec)",
                      help="wait this seconds before creating stream data")

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self.counter = 0
        deps = []

        if self.delay.get() != 0:
            delayId = self._insertFunctionStep('delayStep')
            deps.append(delayId)

        step = None
        if self.setof == SET_OF_MOVIES:
            step = 'createStep'
        elif self.setof == SET_OF_MICROGRAPHS:
            step = 'createStep'
        elif self.setof == SET_OF_RANDOM_MICROGRAPHS:
            step = 'createRandomMicStep'
        elif self.setof == SET_OF_PARTICLES:
            step = 'createParticlesStep'
        elif self.setof == SET_OF_COORDINATES:
            step = 'createCoordinatesStep'
        else:
            raise Exception('Unknown data type')

        if self.setof == SET_OF_PARTICLES:
            self.nDims = int(min(self.nDim, len(self.inputParticles.get())))
            self.group = int(min(self.nDims, self.groups.get()))
            for mic in range(1, int((self.nDims / self.group) +
                                    (self.nDims % self.group > 0) + 1)):
                self._insertFunctionStep(step, prerequisites=deps)

        elif self.setof == SET_OF_COORDINATES:
            self.inputMicrographs = self.inputCoordinates.get().getMicrographs()
            self.nDims = min(self.nDim.get(), len(self.inputMicrographs))
            for micIdx in range(1, self.nDims + 1):
                self._insertFunctionStep(step, micIdx, prerequisites=deps)

        else:
            for mic in range(1, self.nDim.get() + 1):
                self._insertFunctionStep(step, mic, prerequisites=deps)

    # -------------------------- STEPS functions ------------------------
    def delayStep(self):
        time.sleep(self.delay.get())

    def _checkNewItems(self, objSet):
        """ Check for already computed micrograph/movie and
        update the output set. """
        objDict = {}
        newObj = False
        for obj in objSet:
            objDict[obj.getFileName()] = True

        if objDict:
            if objSet.getSize():
                objSet.enableAppend()
                objSet.loadAllProperties()
        else:
            objSet.setStreamState(objSet.STREAM_OPEN)
            acquisition = emobj.Acquisition()
            if self.setof == SET_OF_MICROGRAPHS:
                acquisition.setMagnification(
                    self.inputMics.get().getAcquisition().getMagnification())
                acquisition.setVoltage(
                    self.inputMics.get().getAcquisition().getVoltage())
                acquisition.setSphericalAberration(
                    self.inputMics.get().getAcquisition().getSphericalAberration())
                acquisition.setAmplitudeContrast(
                    self.inputMics.get().getAcquisition().getAmplitudeContrast())
                objSet.setAcquisition(acquisition)
                objSet.setSamplingRate(self.inputMics.get().getSamplingRate())
            elif self.setof == SET_OF_MOVIES:
                acquisition.setMagnification(
                    self.inputMovies.get().getAcquisition().getMagnification())
                acquisition.setVoltage(
                    self.inputMovies.get().getAcquisition().getVoltage())
                acquisition.setSphericalAberration(
                    self.inputMovies.get().getAcquisition().getSphericalAberration())
                acquisition.setAmplitudeContrast(
                    self.inputMovies.get().getAcquisition().getAmplitudeContrast())
                objSet.setAcquisition(acquisition)
                objSet.setSamplingRate(
                    self.inputMovies.get().getSamplingRate())
            else:
                acquisition.setMagnification(self._magnification)
                acquisition.setVoltage(self._voltage)
                acquisition.setSphericalAberration(self._sphericalAberration)
                acquisition.setAmplitudeContrast(self._amplitudeContrast)
                objSet.setAcquisition(acquisition)
                if self.setof == SET_OF_PARTICLES:
                    objSet.setSamplingRate(
                        self.inputParticles.get().getSamplingRate())
                else:
                    objSet.setSamplingRate(self.samplingRate.get())

        if self.setof == SET_OF_MOVIES:
            obj = emobj.Movie()
        elif self.setof == SET_OF_MICROGRAPHS:
            obj = emobj.Micrograph()
        elif self.setof == SET_OF_RANDOM_MICROGRAPHS:
            obj = emobj.Micrograph()
        elif self.setof == SET_OF_PARTICLES:
            obj = emobj.Particle()
        else:
            raise Exception('Unknown data type')

        for k, v in self.dictObj.items():
            if k not in objDict:
                self.counter += 1
                obj.setFileName(k)
                if self.setof != SET_OF_PARTICLES:
                    obj.setMicName(basename(k))
                obj.setObjId(self.counter)
                objSet.append(obj)
                newObj = True

        return objSet, newObj  # why a dictionary, a boolean may be enough

    def _updateOutput(self, objSet):
        if self.setof == SET_OF_MOVIES:
            self._defineOutputs(outputMovies=objSet)
            self._defineTransformRelation(self.inputMovies, objSet)
        elif self.setof == SET_OF_MICROGRAPHS:
            self._defineOutputs(outputMicrographs=objSet)
            self._defineTransformRelation(self.inputMics, objSet)
        elif self.setof == SET_OF_RANDOM_MICROGRAPHS:
            self._defineOutputs(outputMicrographs=objSet)
        elif self.setof == SET_OF_PARTICLES:
            self._defineOutputs(outputParticles=objSet)
            self._defineTransformRelation(self.inputParticles, objSet)

    def _checkProcessedData(self):
        if self.setof == SET_OF_MOVIES:
            objSet = emobj.SetOfMovies(filename=self._getPath('movies.sqlite'))
        elif self.setof == SET_OF_MICROGRAPHS:
            objSet = emobj.SetOfMicrographs(filename=self._getPath('micrographs.sqlite'))
        elif self.setof == SET_OF_RANDOM_MICROGRAPHS:
            objSet = emobj.SetOfMicrographs(filename=self._getPath('micrographs.sqlite'))
        elif self.setof == SET_OF_PARTICLES:
            objSet = emobj.SetOfParticles(filename=self._getPath('particles.sqlite'))
        else:
            raise Exception('Unknown data type')

        newObjSet, newObj = self._checkNewItems(objSet)

        if self.setof == SET_OF_MOVIES:
            newObjSet.setFramesRange(self.inputMovies.get().getFramesRange())

        # check if end ....
        if self.setof != SET_OF_PARTICLES:
            self.nDims = self.nDim.get()
        endObjs = newObjSet.getSize() == self.nDims

        if newObj:
            if endObjs:
                newObjSet.setStreamState(newObjSet.STREAM_CLOSED)
            self._updateOutput(newObjSet)
        newObjSet.close()

    def createStep(self, counter):
        time.sleep(self.getTimeInterval())
        if not ProtCreateStreamData.object or self.setof == \
                SET_OF_MICROGRAPHS or self.setof == SET_OF_MOVIES:

            if self.setof == SET_OF_MOVIES:
                setDim = self.inputMovies.get().getSize()
                for idx, mov in enumerate(self.inputMovies.get()):
                    if idx == (counter - 1) % setDim:
                        newMov = mov.clone()
                        break
                _, ext = splitext(newMov.getFileName())
                ProtCreateStreamData.object = \
                    emlib.image.ImageHandler().read(newMov.getLocation())
                self.name = "movie"

            elif self.setof == SET_OF_MICROGRAPHS:
                setDim = self.inputMics.get().getSize()
                for idx, mic in enumerate(self.inputMics.get()):
                    if idx == (counter-1) % setDim:
                        newMic = mic.clone()
                        break
                _, ext = splitext(newMic.getFileName())
                ProtCreateStreamData.object = \
                    emlib.image.ImageHandler().read(newMic.getLocation())
                self.name = "micro"

        # save file
        destFn = self._getExtraPath("%s_%05d%s" % (self.name, counter, ext))
        ProtCreateStreamData.object.write(destFn)
        self.dictObj[destFn] = True
        self._checkProcessedData()

    def createParticlesStep(self):
        self.name = "particle"
        time.sleep(self.getTimeInterval())

        for idx, p in enumerate(self.inputParticles.get()):
            if ((idx > self.counter-1) and (idx < self.nDims) and
                    (idx <= self.counter-1 + self.group)):
                newP = p.clone()
                ProtCreateStreamData.object = \
                    emlib.image.ImageHandler().read(newP.getLocation())
                destFn = self._getExtraPath("%s_%05d" % (self.name, idx))
                ProtCreateStreamData.object.write(destFn)
                self.dictObj[destFn] = True
        self._checkProcessedData()

    def createCoordinatesStep(self, micIdx):
        self.name = "coordinate"
        time.sleep(self.getTimeInterval())
        inputCoordinates = self.inputCoordinates.get()
        micrographs = inputCoordinates.getMicrographs()
        if not hasattr(self, "outputCoordinates"):
            self.outputCoordinates = self._createSetOfCoordinates(micrographs)
            self.outputCoordinates.copyInfo(inputCoordinates)
            self.outputCoordinates.setBoxSize(inputCoordinates.getBoxSize())
            self.outputCoordinates.setStreamState(emobj.SetOfParticles.STREAM_OPEN)
            self._defineOutputs(outputCoordinates=self.outputCoordinates)
            self._defineSourceRelation(inputCoordinates, self.outputCoordinates)

        for idx, mic in enumerate(micrographs):
            if idx == micIdx - 1:
                newMic = mic.clone()
                for newCoord in self.inputCoordinates.get().iterCoordinates(newMic):
                    self.outputCoordinates.append(newCoord)
                if micIdx == self.nDims:
                    self.outputCoordinates.setStreamState(emobj.SetOfParticles.STREAM_CLOSED)
                self.outputCoordinates.write()
                self._store(self.outputCoordinates)
                break


    def createRandomMicStep(self, mic):
        from pwem import Domain
        time.sleep(self.getTimeInterval())
        getEnviron = Domain.importFromPlugin('xmipp3', 'Plugin',
                                             doRaise=True).getEnviron

        # create image
        img = emlib.Image()
        img.setDataType(emlib.DT_FLOAT)
        img.resize(self.xDim, self.yDim)
        img.initRandom(0., 1., emlib.XMIPP_RND_UNIFORM)
        baseFn = self._getExtraPath(self._singleImageFn)
        img.write(baseFn)

        md1 = emlib.MetaData()
        md1.setColumnFormat(False)
        idctf = md1.addObject()

        baseFnCtf = self._getTmpPath("ctf_%d.param" % mic)
        baseFnImageCTF = self._getExtraPath("imageCTF_%d.xmp" % mic)

        md1.setValue(emlib.MDL_CTF_SAMPLING_RATE, 1., idctf)
        md1.setValue(emlib.MDL_CTF_VOLTAGE, 200., idctf)
        defocus = 20000 + 10000 * random.random()
        udefocus = defocus + 1000 * random.random()
        vdefocus = defocus + 1000 * random.random()
        if udefocus < vdefocus:
            aux = vdefocus
            vdefocus = udefocus
            udefocus = aux
        md1.setValue(emlib.MDL_CTF_DEFOCUSU, udefocus, idctf)
        md1.setValue(emlib.MDL_CTF_DEFOCUSV, vdefocus, idctf)
        md1.setValue(emlib.MDL_CTF_DEFOCUS_ANGLE, 180.0 * random.random(),
                     idctf)
        md1.setValue(emlib.MDL_CTF_CS, 2., idctf)
        md1.setValue(emlib.MDL_CTF_Q0, 0.07, idctf)
        md1.setValue(emlib.MDL_CTF_K, 1., idctf)

        md1.write(baseFnCtf)

        # apply ctf
        args = " -i %s" % baseFn
        args += " -o %s" % baseFnImageCTF
        args += " -f ctf %s" % baseFnCtf
        args += " --sampling %f" % self.samplingRate
        self.runJob("xmipp_transform_filter", args, env=getEnviron())
        self.dictObj[baseFnImageCTF] = True
        self._checkProcessedData()

    # -------------------------- UTILS functions -----------------------------
    def getTimeInterval(self):
        baseInterval = self.creationInterval.get()
        interval = random.randint(baseInterval, baseInterval+self.extraRandomInterval.get())
        return interval


    # -------------------------- INFO functions ------------------------------
    def _validate(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []
