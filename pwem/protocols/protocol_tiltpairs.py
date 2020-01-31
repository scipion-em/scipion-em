# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano         (ldelcano@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
"""
This module contains protocols classes related to Random Conical Tilt.
"""

import sys
from os.path import basename
from glob import glob

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params

import pwem.constants as emcts
from pwem import emlib
import pwem.objects as emobj

from .protocol_import import ProtImportFiles


class ProtImportMicrographsTiltPairs(ProtImportFiles):
    """Protocol to import untilted-tilted pairs of micrographs in the project"""
    
    _className = 'Micrograph'
    _label = 'import tilted micrographs'
    OUTPUT_NAME = "outputMicrographsTiltPair"
        
    # --------------------------- DEFINE param functions ----------------------
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('patternUntilted', params.PathParam,
                      label=pwutils.Message.LABEL_PATTERNU,
                      help=pwutils.Message.TEXT_PATTERN)
        form.addParam('patternTilted', params.PathParam,
                      label=pwutils.Message.LABEL_PATTERNT,
                      help=pwutils.Message.TEXT_PATTERN)
        form.addParam('copyToProj', params.BooleanParam,
                      label=pwutils.Message.LABEL_COPYFILES, default=False)
        form.addParam('voltage', params.FloatParam, default=200,
                      label=pwutils.Message.LABEL_VOLTAGE, help=pwutils.Message.TEXT_VOLTAGE)
        form.addParam('sphericalAberration', params.FloatParam, default=2.26,
                      label=pwutils.Message.LABEL_SPH_ABERRATION,
                      help=pwutils.Message.TEXT_SPH_ABERRATION)
        form.addParam('ampContrast', params.FloatParam, default=0.1,
                      label=pwutils.Message.LABEL_AMPLITUDE,
                      help=pwutils.Message.TEXT_AMPLITUDE)
        form.addParam('samplingRateMode', params.EnumParam,
                      default=emcts.SAMPLING_FROM_IMAGE,
                      label=pwutils.Message.LABEL_SAMP_MODE,
                      help=pwutils.Message.TEXT_SAMP_MODE,
                      choices=[pwutils.Message.LABEL_SAMP_MODE_1,
                               pwutils.Message.LABEL_SAMP_MODE_2])
        form.addParam('samplingRate', params.FloatParam, default=1,
                      label=pwutils.Message.LABEL_SAMP_RATE,
                      help=pwutils.Message.TEXT_SAMP_RATE,
                      condition='samplingRateMode==%d' % emcts.SAMPLING_FROM_IMAGE)
        form.addParam('magnification', params.IntParam, default=50000,
                      label=pwutils.Message.LABEL_MAGNI_RATE,
                      help=pwutils.Message.TEXT_MAGNI_RATE,
                      condition='samplingRateMode==%d' % emcts.SAMPLING_FROM_SCANNER)
        form.addParam('scannedPixelSize', params.FloatParam, default=7.0,
                      label=pwutils.Message.LABEL_SCANNED,
                      condition='samplingRateMode==%d' % emcts.SAMPLING_FROM_SCANNER)

    # --------------------------- INSERT steps functions ----------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('createTiltedPairsStep')
        
    # --------------------------- STEPS functions -----------------------------
    
    def importMicrographs(self, pattern, suffix, voltage, sphericalAberration, amplitudeContrast):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        filePaths = glob(pwutils.expandPattern(pattern))
        
        # imgSet = SetOfMicrographs(filename=self.micsPairsSqlite, prefix=suffix)
        imgSet = self._createSetOfMicrographs(suffix=suffix)
        acquisition = imgSet.getAcquisition()
        # Setting Acquisition properties
        acquisition.setVoltage(voltage)
        acquisition.setSphericalAberration(sphericalAberration)
        acquisition.setAmplitudeContrast(amplitudeContrast)
        
        # Call a function that should be implemented by each subclass
        self._setOtherPars(imgSet)
        
        outFiles = [imgSet.getFileName()]
        imgh = emlib.image.ImageHandler()
        img = imgSet.ITEM_TYPE()
        n = 1
        size = len(filePaths)
        
        filePaths.sort()
        
        for i, fn in enumerate(filePaths):
            # ext = os.path.splitext(basename(f))[1]
            dst = self._getExtraPath(basename(fn))
            if self.copyToProj:
                pwutils.copyFile(fn, dst)
            else:
                pwutils.createLink(fn, dst)
            
            if n > 1:
                for index in range(1, n+1):
                    img.cleanObjId()
                    img.setFileName(dst)
                    img.setIndex(index)
                    imgSet.append(img)
            else:
                img.cleanObjId()
                img.setFileName(dst)
                # Fill the micName if img is a Micrograph.
                self._fillMicName(img, fn, pattern)
                imgSet.append(img)
            outFiles.append(dst)
            
            sys.stdout.write("\rImported %d/%d" % (i+1, size))
            sys.stdout.flush()
            
        print("\n")
        
        imgSet.write()

        return imgSet
    
    def createTiltedPairsStep(self):
        args = {} 
        self.micsPairsSqlite = self._getPath('micrographs_pairs.sqlite')
        pwutils.cleanPath(self.micsPairsSqlite)  # Delete if exists
        
        micsTiltPair = emobj.MicrographsTiltPair(filename=self.micsPairsSqlite)
        micsU = self.importMicrographs(self.patternUntilted.get(), 'Untilted',
                                       self.voltage.get(),
                                       self.sphericalAberration.get(),
                                       self.ampContrast.get())
        micsT = self.importMicrographs(self.patternTilted.get(), 'Tilted',
                                       self.voltage.get(),
                                       self.sphericalAberration.get(),
                                       self.ampContrast.get())
        
        micsTiltPair.setUntilted(micsU)
        micsTiltPair.setTilted(micsT)
        
        for micU, micT in zip(micsU, micsT):
            micsTiltPair.append(emobj.TiltPair(micU, micT))
        
        args[self.OUTPUT_NAME] = micsTiltPair
        self._defineOutputs(**args)
    
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        if not self.patternUntilted.get() or not self.patternTilted.get():
            errors.append(pwutils.Message.ERROR_PATTERN_EMPTY)
        else:
            filePathsUntilted = glob(pwutils.expandPattern(self.patternUntilted.get()))
            filePathsTilted = glob(pwutils.expandPattern(self.patternTilted.get()))
        
            if len(filePathsUntilted) == 0 or len(filePathsTilted) == 0:
                errors.append(pwutils.Message.ERROR_PATTERN_FILES)

        return errors
    
    def _summary(self):
        summary = []
        output = self._getOutput()
        if not output:
            summary.append("Output not ready yet.") 
            if self.copyToProj:
                summary.append("*Warning*: Import step could take a long "
                               "time due to the images are copying in the project.")
        else:
            summary.append("Import of %d " % output.getTilted().getSize() + self._className + "s tilted from %s" % self.patternTilted.get())
            summary.append("Import of %d " % output.getUntilted().getSize() + self._className + "s untilted from %s" % self.patternUntilted.get())
            summary.append("Sampling rate : %0.2f Å/px" % self.samplingRate.get())
        
        return summary
    
    def _methods(self):
        methods = []
        output = self._getOutput()
        if output:
            methods.append("Imported %d micrographs tilted pairs %s " % (output.getTilted().getSize(), self.getObjectTag(output)))
            methods.append("with a sampling rate : %0.2f Å/px." % self.samplingRate.get())
            
        return methods
    
    # --------------------------- UTILS functions -----------------------------
    def getFiles(self):
        return getattr(self, self._getOutput().getFiles())
    
    def _getOutput(self):
        return getattr(self, self.OUTPUT_NAME, None)
        
    def _setOtherPars(self, micSet):
        micSet.getAcquisition().setMagnification(self.magnification.get())
        
        if self.samplingRateMode == emcts.SAMPLING_FROM_IMAGE:
            micSet.setSamplingRate(self.samplingRate.get())
        else:
            micSet.setScannedPixelSize(self.scannedPixelSize.get())
            
    def _fillMicName(self, img, filename, pattern):
        from pwem.objects import Micrograph
        if isinstance(img, Micrograph):
            filePaths = self.getMatchFiles(pattern=pattern)
            commPath = pwutils.commonPath(filePaths)
            micName = filename.replace(commPath + "/", "").replace("/", "_")
            img.setMicName(micName)
