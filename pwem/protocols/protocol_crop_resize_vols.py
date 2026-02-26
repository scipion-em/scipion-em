# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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


import os
import numpy as np
import scipy.fft as fft
import mrcfile

import pyworkflow.protocol.params as params
from pyworkflow import VERSION_3_0

from pwem.protocols import EMProtocol
import pwem.objects as emobj
 
OUTPUT_VOLUME = 'volume'


class ProtCropResizeVols(EMProtocol):
    """ Cropping and resizing volumes. The protocol can resize the volumes in to manners:\n
    1) By setting a target samling rate.\n
    2) By setting the output box dimensions.\n
    The protocol also allows to crop or pad the volume, setting the output boxsize. This
    option is applied after the previous resizing option.
    """
    _label = 'crop and/or resize'
    _program = ""
    _possibleOutputs = {OUTPUT_VOLUME: emobj.Volume}

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
        self._insertFunctionStep(self.cropResizeStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------        
    def cropResizeStep(self):
        inVol = self.inVolume.get()
        inBoxsize = inVol.getDim()

        if self.doResize.get():
            if self.resizeOption.get() == self.RESIZE_SAMPLINGRATE:
                newSampling = self.resizeSamplingRate.get()
                sampling = inVol.getSamplingRate()
                targetShape = None
                scaleFactor = sampling/newSampling
                boxsize = round(inBoxsize[0]*scaleFactor)
                targetShape = (boxsize, boxsize, boxsize)
            elif self.resizeOption.get() == self.RESIZE_DIMENSIONS:
                newboxsize = self.resizeDim.get()
                targetShape = (newboxsize, newboxsize, newboxsize)

            self.scaleCropOperation(inVol, targetShape=targetShape)
        else:
            boxsize = self.cropPadDim.get()
            fnVol = inVol.getFileName()
            
            self.readCropPadWrite(fnVol, boxsize)

            if inVol.hasHalfMaps():
                fnHalf1, fnHalf2 = inVol.getHalfMaps().split(',')
                self.readCropPadWrite(fnHalf1, boxsize)
                self.readCropPadWrite(fnHalf2, boxsize)

    def readCropPadWrite(self, fnIn, boxsize)  :
        fnOut = self._getExtraPath(os.path.basename(fnIn))
        vol = self.readMRCfile(fnIn)
        cropVol = self.cropPad(vol, boxsize)
        self.writeMRCfile(cropVol, fnOut)
            
    def scaleCrop(self, fnVol, targetShape, fnOut):
        vol = self.readMRCfile(fnVol)
        resizedVol = self.fourierRescale3D(vol, targetShape)
        if self.doCropPad.get():
            boxsize = self.cropPadDim.get()
            resizedVol = self.cropPad(resizedVol, boxsize)
        self.writeMRCfile(resizedVol, fnOut)

    def scaleCropOperation(self, inVol, targetShape):
        fnVol = inVol.getFileName()
        fnBase = inVol.getBaseName()
        fnOut = self._getExtraPath(fnBase)

        self.scaleCrop(fnVol, targetShape, fnOut)

        if inVol.hasHalfMaps():
            fnHalf1, fnHalf2 = inVol.getHalfMaps().split(',')
            fnOut1 = self._getExtraPath(os.path.basename(fnHalf1))
            fnOut2 = self._getExtraPath(os.path.basename(fnHalf2))

            self.scaleCrop(fnHalf1, targetShape, fnOut1)
            self.scaleCrop(fnHalf2, targetShape, fnOut2)


    def createOutputStep(self):
        inVol = self.inVolume.get()
        outVol = emobj.Volume()
        outVol.copy(inVol)

        oldSampling = inVol.getSamplingRate()
        if self.doResize.get():
            if self.resizeOption.get() == self.RESIZE_SAMPLINGRATE:
                samplingRate = self.resizeSamplingRate.get()
            else:
                inBoxsize = inVol.getDim()[0]
                newboxsize = self.resizeDim.get()
                samplingRate = (inBoxsize/newboxsize)*oldSampling
        else:
            samplingRate = oldSampling

        outVol.setSamplingRate(samplingRate)
        outVol.setOrigin()

        fnIn = inVol.getBaseName()
        fnOut = self._getExtraPath(fnIn)
        outVol.setFileName(fnOut)

        if inVol.hasHalfMaps():
            fnHalf1, fnHalf2 = inVol.getHalfMaps().split(',')
            fnHalf1 = self._getExtraPath(os.path.basename(fnHalf1))
            fnHalf2 = self._getExtraPath(os.path.basename(fnHalf2))
            outVol.setHalfMaps([fnHalf1, fnHalf2])

        self._defineOutputs(**{OUTPUT_VOLUME: outVol})
        self._defineSourceRelation(inVol, outVol)

    def readMRCfile(self, fn):
        return mrcfile.open(fn, mode='r+').data
    
    def writeMRCfile(self, vol, fnOut):
        with mrcfile.new(fnOut) as mrc:
            mrc.set_data(vol)

    def fourierRescale3D(self, arr, targetShape):

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

    def cropPad(self, vol, boxsize):
        boxsize = int(boxsize)
        
        oldShape = vol.shape
        newShape = (boxsize, boxsize, boxsize)
        
        if oldShape == newShape:
            return np.copy(vol)
            
        newVol = np.zeros(newShape, dtype=vol.dtype)
        
        slicesOld = []
        slicesNew = []
        
        for oldDim, newDim in zip(oldShape, newShape):
            copySize = min(oldDim, newDim)
            
            startOld = max(0, (oldDim - newDim) // 2)
            slicesOld.append(slice(startOld, startOld + copySize))
            
            startNew = max(0, (newDim - oldDim) // 2)
            slicesNew.append(slice(startNew, startNew + copySize))

        newVol[tuple(slicesNew)] = vol[tuple(slicesOld)]
        
        return newVol

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        message = []
        if (not self.doResize.get() and not self.doCropPad.get()):
            message.append('An option must be selected for resizing.')
        return message
