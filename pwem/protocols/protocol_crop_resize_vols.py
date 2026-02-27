# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
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
import logging
from typing import Tuple, Optional, Union
import numpy as np
import scipy.fft as fft
import mrcfile
import pyworkflow.protocol.params as params
from pwem.convert.headers import setMRCSamplingRate
from pwem.emlib.image.image_readers import ImageReadersRegistry, MRCImageReader
from pwem.objects import Volume
from pwem.protocols import EMProtocol
from pyworkflow.utils import removeBaseExt, createLink, cyanStr

logger = logging.getLogger(__name__)
OUTPUT_VOLUME = 'volume'
MRC_EXT = '.mrc'


class ProtCropResizeVols(EMProtocol):
    """ Cropping and resizing volumes. The protocol can resize the volumes in to manners:\n
    1) By setting a target sampling rate.\n
    2) By setting the output box dimensions.\n
    The protocol also allows to crop or pad the volume, setting the output boxsize. This
    option is applied after the previous resizing option.
    """
    _label = 'crop and/or resize volumes'
    _program = ""
    _possibleOutputs = {OUTPUT_VOLUME: Volume}

    RESIZE_SAMPLINGRATE = 0
    RESIZE_DIMENSIONS = 1

    WINDOW_OP_CROP = 0
    WINDOW_OP_WINDOW = 1

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inVolume', params.PointerParam, label="Input Volume", important=True,
                      pointerClass='Volume',
                      help='This is the volume to be rezised')
        form.addParam('doResize', params.BooleanParam, label="Scale volume?",
                      default=True,
                      help='This option allows to resize the input volume according to a given sampling rate' \
                           'or a given boxsize.')
        form.addParam('resizeOption', params.EnumParam,
                      choices=['Sampling Rate', 'Dimensions'],
                      default=self.RESIZE_SAMPLINGRATE,
                      label="Resize option", display=params.EnumParam.DISPLAY_COMBO,
                      help='Select an option to resize the images: \n '
                           '_Sampling Rate_: Set the desire sampling rate to resize. \n'
                           '_Dimensions_: Set the output dimensions. Resize operation can be done in Fourier space.\n')
        form.addParam('resizeSamplingRate', params.FloatParam, default=1.0,
                      condition='doResize and resizeOption==%d' % self.RESIZE_SAMPLINGRATE,
                      label='Target sampling rate (Å/px)',
                      help='This is the output sampling rate.')
        form.addParam('resizeDim', params.IntParam, default=0,
                      condition='doResize and resizeOption==%d' % self.RESIZE_DIMENSIONS,
                      allowsPointers=True,
                      label='New size (px)',
                      help='Size in pixels of the particle images <x> <y=x> <z=x>.')
        form.addParam('doCropPad', params.BooleanParam, label="Cropping/padding volume?",
                      default=False,
                      help='This option allows to crop or pad the volumen. If yes, a new' \
                           'form parameter will apper to define the boxsize.')
        form.addParam('cropPadDim', params.IntParam, default=0,
                      condition='doCropPad',
                      allowsPointers=True,
                      label='New boxsize (px)',
                      help='Size in pixels of the output volume. The volume will be cropped'
                           'in real space to have the set box size. If the boxsize is larger' \
                           'than the original boxsize, the padding will be carried out with' \
                           'zeros.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.cropResizeStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        inVolume = self.inVolume.get()
        inVolFName = inVolume.getFileName()
        outVolFName = self._getConvertedOrLinkedName(inVolFName)
        if inVolFName.endswith(MRC_EXT):
            logger.info(cyanStr(f'Making a link {inVolFName} -> {outVolFName}...'))
            createLink(inVolFName, outVolFName)
            if inVolume.hasHalfMaps():
                fnHalf1, fnHalf2 = inVolume.getHalfMaps().split(',')
                # Even
                outHalf1 = self._getConvertedOrLinkedName(fnHalf1)
                createLink(fnHalf1, outHalf1)
                # Odd
                outHalf2 = self._getConvertedOrLinkedName(fnHalf2)
                createLink(fnHalf2, outHalf2)
        else:
            self._convertToMrc(inVolFName)
            if inVolume.hasHalfMaps():
                fnHalf1, fnHalf2 = inVolume.getHalfMaps().split(',')
                # Even
                self._convertToMrc(fnHalf1)
                # Odd
                self._convertToMrc(fnHalf2)

    def cropResizeStep(self):
        logger.info(cyanStr(f'Cropping and/or resizing...'))
        inVol = self.inVolume.get()
        if self.doResize.get():
            targetShape = self._getTargetShape(inVol)
            self._scaleCropOperation(inVol, targetShape=targetShape)
        else:
            boxsize = self.cropPadDim.get()
            self._readCropPadWrite(inVol, boxsize)

    def createOutputStep(self):
        logger.info(cyanStr(f'Registering the results...'))
        inVol = self.inVolume.get()
        outVol = Volume()
        outVol.copy(inVol)

        oldSampling = inVol.getSamplingRate()
        if self.doResize.get():
            if self.resizeOption.get() == self.RESIZE_SAMPLINGRATE:
                samplingRate = self.resizeSamplingRate.get()
            else:
                inBoxsize = inVol.getDim()[0]
                newboxsize = self.resizeDim.get()
                samplingRate = (inBoxsize / newboxsize) * oldSampling
        else:
            samplingRate = oldSampling

        outVol.setSamplingRate(samplingRate)
        outVol.setOrigin()

        fnOut = self._getResultVolumeFName(inVol.getFileName())
        setMRCSamplingRate(fnOut, samplingRate)
        outVol.setFileName(fnOut)

        if inVol.hasHalfMaps():
            fnHalf1, fnHalf2 = inVol.getHalfMaps().split(',')
            fnHalf1 = self._getResultVolumeFName(fnHalf1)
            fnHalf2 = self._getResultVolumeFName(fnHalf2)
            outVol.setHalfMaps([fnHalf1, fnHalf2])

        self._defineOutputs(**{OUTPUT_VOLUME: outVol})
        self._defineSourceRelation(inVol, outVol)

    # --------------------------- UTILS functions -----------------------------
    def _getTargetShape(self, inVol: Volume) -> Tuple[int, int, int]:
        if self.resizeOption.get() == self.RESIZE_SAMPLINGRATE:
            inBoxsize = inVol.getDim()
            newSampling = self.resizeSamplingRate.get()
            sampling = inVol.getSamplingRate()
            scaleFactor = sampling / newSampling
            boxsize = round(inBoxsize[0] * scaleFactor)
            return boxsize, boxsize, boxsize
        else:  # self.resizeOption.get() == self.RESIZE_DIMENSIONS:
            newboxsize = self.resizeDim.get()
            return newboxsize, newboxsize, newboxsize

    def _scaleCropOperation(self,
                            inVol: Volume,
                            targetShape: Tuple[int, int, int]) -> None:

        inVolFName, outVolFName = self._getInOutFileNames(inVol)
        self._scaleCrop(inVolFName, targetShape, outVolFName)

        if inVol.hasHalfMaps():
            # Even
            inVolFName, outVolFName = self._getInOutFileNames(inVol, even=True)
            self._scaleCrop(inVolFName, targetShape, outVolFName)
            # Odd
            inVolFName, outVolFName = self._getInOutFileNames(inVol, even=False)
            self._scaleCrop(inVolFName, targetShape, outVolFName)

    def _scaleCrop(self,
                   fnVol: str,
                   targetShape: Tuple[int, int, int],
                   fnOut: str) -> None:
        resizedVol = self._fourierRescale3D(fnVol, targetShape)
        if self.doCropPad.get():
            boxsize = self.cropPadDim.get()
            self._cropPad(resizedVol, boxsize, fnOut)
        else:
            self.writeMRCfile(resizedVol, fnOut)

    def _fourierRescale3D(self,
                          fnVol: str,
                          targetShape: Tuple[int, int, int]) -> np.ndarray:
        arr = self.readMRCfile(fnVol)
        oldShape = arr.shape
        targetShape = tuple(int(dim) for dim in targetShape)

        if oldShape == targetShape:
            return np.copy(arr)

        myfft = fft.fftn(arr)
        myfftShifted = fft.fftshift(myfft)

        newfftShifted = np.zeros(targetShape, dtype=myfft.dtype)

        slicesOld = []
        slicesNew = []

        for o, n in zip(oldShape, targetShape):
            min_dim = min(o, n)

            c_o = o // 2
            c_n = n // 2

            halflLeft = min_dim // 2
            halfRight = min_dim - halflLeft

            slicesOld.append(slice(c_o - halflLeft, c_o + halfRight))
            slicesNew.append(slice(c_n - halflLeft, c_n + halfRight))

        newfftShifted[tuple(slicesNew)] = myfftShifted[tuple(slicesOld)]

        newfft = fft.ifftshift(newfftShifted)
        resisedVol = fft.ifftn(newfft)

        # It is neccesary to scale the amplitudes due to the total energy changes
        scaleFactor = np.prod(targetShape) / np.prod(oldShape)
        resisedVol *= scaleFactor
        resisedVol = np.real(resisedVol)

        return resisedVol

    def _cropPad(self,
                 inVolData: Union[str, np.ndarray],
                 boxsize: int,
                 outVolFName: str) -> None:
        volData = self.readMRCfile(inVolData) if type(inVolData) == str else inVolData
        oldShape = volData.shape
        newShape = (boxsize, boxsize, boxsize)
        if oldShape == newShape:
            self.writeMRCfile(np.copy(volData), outVolFName)
            return

        newVol = np.zeros(newShape, dtype=volData.dtype)
        slicesOld = []
        slicesNew = []

        for oldDim, newDim in zip(oldShape, newShape):
            copySize = min(oldDim, newDim)

            startOld = max(0, (oldDim - newDim) // 2)
            slicesOld.append(slice(startOld, startOld + copySize))

            startNew = max(0, (newDim - oldDim) // 2)
            slicesNew.append(slice(startNew, startNew + copySize))

        newVol[tuple(slicesNew)] = volData[tuple(slicesOld)]
        self.writeMRCfile(newVol, outVolFName)

    def _readCropPadWrite(self, inVol: Volume, boxsize: int):
        inVolFName, outVolFName = self._getInOutFileNames(inVol)
        self._cropPad(inVolFName, boxsize, outVolFName)
        if inVol.hasHalfMaps():
            # Even
            inVolFName, outVolFName = self._getInOutFileNames(inVol, even=True)
            self._cropPad(inVolFName, boxsize, outVolFName)
            # Odd
            inVolFName, outVolFName = self._getInOutFileNames(inVol, even=False)
            self._cropPad(inVolFName, boxsize, outVolFName)

    @staticmethod
    def readMRCfile(inVolFName: str) -> np.ndarray:
        return MRCImageReader.open(inVolFName)

    @staticmethod
    def writeMRCfile(vol: np.ndarray, fnOut: str) -> None:
        with mrcfile.new(fnOut) as mrc:
            mrc.set_data(vol)

    def _getConvertedOrLinkedName(self, inVolFName: str) -> str:
        return self._getTmpPath(f'{removeBaseExt(inVolFName)}{MRC_EXT}')

    def _getResultVolumeFName(self, inVolFName: str) -> str:
        return self._getExtraPath(f'{removeBaseExt(inVolFName)}{MRC_EXT}')

    def _convertToMrc(self, inFileName: str):
        outVolFName = self._getConvertedOrLinkedName(inFileName)
        logger.info(cyanStr(f'Converting {inFileName} to {outVolFName}...'))
        imgStack = ImageReadersRegistry.open(inFileName)
        ImageReadersRegistry.write(imgStack, outVolFName)

    def _getInOutFileNames(self, inVol: Volume, even: Optional[bool] = None) -> Tuple[str, str]:
        if even is None:
            fnVol = inVol.getFileName()
        elif even:
            fnVol = inVol.getHalfMaps().split(',')[0]
        else:
            fnVol = inVol.getHalfMaps().split(',')[1]
        inVolFName = self._getConvertedOrLinkedName(fnVol)
        outVolFName = self._getResultVolumeFName(fnVol)
        return inVolFName, outVolFName

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        message = []
        if (not self.doResize.get() and not self.doCropPad.get()):
            message.append('An option must be selected for resizing.')
        return message
