# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Vahid Abrishami (vabrishami@cnb.csic.es)
# *              Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
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
import os
import warnings
from math import ceil

import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.gui.plotter import Plotter
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as pwcts

import pwem.objects as emobj
from pwem import emlib

from pwem.protocols import ProtProcessMovies

OUT_MICS = "outputMicrographs"
OUT_MICS_DW = "outputMicrographsDoseWeighted"
OUT_MOVIES = "outputMovies"
OUT_MICS_ODD = 'outputMicrographsOdd'
OUT_MICS_EVEN = 'outputMicrographsEven'


class ProtAlignMovies(ProtProcessMovies):
    """
    Base class for movie alignment protocols such as:
    motioncorr, crosscorrelation and optical flow

    Alignment parameters are defined in common. For example,
    the frames range used for alignment and final sum, the binning factor
    or the cropping options (region of interest)
    """
    _possibleOutputs = {OUT_MICS: emobj.SetOfMicrographs,
                        OUT_MICS_DW: emobj.SetOfMicrographs,
                        OUT_MICS_ODD: emobj.SetOfMicrographs,
                        OUT_MICS_EVEN: emobj.SetOfMicrographs,
                        OUT_MOVIES: emobj.SetOfMovies}

    # Even / odd functionality
    evenOddCapable = False

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)
        self._defineAlignmentParams(form)

    def _defineAlignmentParams(self, form):
        group = form.addGroup('Alignment')
        line = group.addLine('Frames to ALIGN',
                             help='Frames range to ALIGN on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to align, it means that you will '
                                  'align until the last frame of the movie.')
        line.addParam('alignFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('alignFrameN', params.IntParam, default=0,
                      label='to')
        group.addParam('useAlignToSum', params.BooleanParam, default=True,
                       label='Use ALIGN frames range to SUM?',
                       help="If *Yes*, the same frame range will be used to "
                            "ALIGN and to SUM. If *No*, you can selected a "
                            "different range for SUM (must be a subset).")
        line = group.addLine('Frames to SUM', condition="not useAlignToSum",
                             help='Frames range to SUM on each movie. The '
                                  'first frame is 1. If you set 0 in the final '
                                  'frame to sum, it means that you will sum '
                                  'until the last frame of the movie.')
        line.addParam('sumFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', params.IntParam, default=0,
                      label='to')
        group.addParam('binFactor', params.FloatParam, default=1.,
                       label='Binning factor',
                       help='1x or 2x. Bin stack before processing.')

        line = group.addLine('Crop offsets (px)')
        line.addParam('cropOffsetX', params.IntParam, default=0, label='X')
        line.addParam('cropOffsetY', params.IntParam, default=0, label='Y')

        line = group.addLine('Crop dimensions (px)',
                             help='How many pixels to crop from offset\n'
                                  'If equal to 0, use maximum size.')
        line.addParam('cropDimX', params.IntParam, default=0, label='X')
        line.addParam('cropDimY', params.IntParam, default=0, label='Y')

        form.addParam('doSaveAveMic', params.BooleanParam, default=True,
                      label="Save aligned micrograph",
                      expertLevel=pwcts.LEVEL_ADVANCED)

        form.addParam('doSaveMovie', params.BooleanParam, default=False,
                      label="Save movie", expertLevel=pwcts.LEVEL_ADVANCED,
                      help="Save Aligned movie")

        if self.evenOddCapable:
            form.addParam('splitEvenOdd', params.BooleanParam,
                          default=False,
                          label='Split & sum odd/even frames?',
                          expertLevel=params.LEVEL_ADVANCED,
                          help='Generate odd and even sums using odd and even frames '
                               'respectively when this option is enabled.')

    # --------------------------- STEPS functions ----------------------------
    def createOutputStep(self):
        # validate that we have some output movies
        for outName, out in self.iterOutputAttributes():
            output = out
            break

        if output.getSize() == 0 and len(self.listOfMovies) != 0:
            raise Exception("All movies failed, didn't create outputMicrographs."
                            "Please review movie processing steps above.")
        elif output.getSize() < len(self.listOfMovies):
            self.warning(pwutils.yellowStr("WARNING - Failed to align %d movies."
                                           % (len(self.listOfMovies) - output.getSize())))

    def _loadOutputSet(self, SetClass, baseName, fixSampling=True):
        """
        Load the output set if it exists or create a new one.
        fixSampling: correct the output sampling rate if binning was used,
        except for the case when the original movies are kept and shifts
        refers to that one.
        """
        setFile = self._getPath(baseName)

        if os.path.exists(setFile) and os.path.getsize(setFile) > 0:
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

            inputMovies = self.inputMovies.get()
            outputSet.copyInfo(inputMovies)

            if fixSampling:
                newSampling = inputMovies.getSamplingRate() * self._getBinFactor()
                outputSet.setSamplingRate(newSampling)

        return outputSet

    def _updateOutputMicSet(self, newDone, sqliteFn, getOutputMicName,
                            outputName, streamMode):
        """ Updated the output micrographs set with new items found. """
        micSet = self._loadOutputSet(emobj.SetOfMicrographs, sqliteFn)
        doneFailed = []

        for movie in newDone:
            mic = micSet.ITEM_TYPE()
            mic.copyObjId(movie)
            mic.setMicName(movie.getMicName())
            # The subclass protocol is responsible for generating the output
            # micrograph file in the extra path with the required name
            extraMicFn = self._getExtraPath(getOutputMicName(movie))
            mic.setFileName(extraMicFn)
            if not os.path.exists(extraMicFn):
                print(pwutils.yellowStr("WARNING: Micrograph %s was not generated, "
                                        "can't add it to output set." % extraMicFn))
                doneFailed.append(movie)
                continue
            # Tolerate errors here. Usually here some plots are generated.
            try:
                self._preprocessOutputMicrograph(mic, movie)
            except Exception as e:
                self.error("Couldn't prepare output details: %s" % e)
                doneFailed.append(movie)
                continue

            micSet.append(mic)

        self._updateOutputSet(outputName, micSet, streamMode)
        if doneFailed:
            self._writeFailedList(doneFailed)

        if self._firstTimeOutput:
            # We consider that Movies are 'transformed' into the Micrographs
            # This will allow to extend the CTF associated to a set of
            # micrographs to another set of micrographs generated from a
            # different movie alignment
            self._defineTransformRelation(self.inputMovies, micSet)

    def _updateOutputMovieSet(self, newDone, streamMode):
        saveMovie = self.getAttributeValue('doSaveMovie', False)
        movieSet = self._loadOutputSet(emobj.SetOfMovies, 'movies.sqlite',
                                       fixSampling=saveMovie)

        # If we need to save the movies
        if saveMovie:
            movieSet.setGain(None)
            movieSet.setDark(None)

        for movie in newDone:
            try:
                newMovie = self._createOutputMovie(movie)
                if newMovie.getAlignment().getShifts()[0]:
                    movieSet.append(newMovie)
                else:
                    print(pwutils.yellowStr("WARNING: Movie %s has empty alignment "
                                            "data, can't add it to output set."
                                            % movie.getFileName()))

            # Warn about any exception creating the movie
            except Exception as e:
                print(pwutils.redStr("ERROR: Movie %s couldn't be "
                                     "added to the output set.\n%s"
                                     % (movie.getFileName(), e)))

        self._updateOutputSet(OUT_MOVIES, movieSet, streamMode)

        if self._firstTimeOutput:
            # Probably is a good idea to store a cached summary for the
            # first resulting movie of the processing.
            self._storeSummary(newDone[0])
            # If the movies are not written out, then dimensions can be
            # copied from the input movies
            if not saveMovie:
                movieSet.setDim(self.inputMovies.get().getDim())
            self._defineTransformRelation(self.inputMovies, movieSet)

    def _updateOutputSets(self, newDone, streamMode):
        if self._createOutputMovies():
            self._updateOutputMovieSet(newDone, streamMode)

        if self._createOutputMicrographs():
            self._updateOutputMicSet(newDone, 'micrographs.sqlite',
                                self._getOutputMicName,
                                OUT_MICS, streamMode)

        if self._createOutputWeightedMicrographs():
            self._updateOutputMicSet(newDone, 'micrographs_dose-weighted.sqlite',
                                self._getOutputMicWtName,
                                OUT_MICS_DW, streamMode)

        if self._doSplitEvenOdd():
            self._updateOutputMicSet(newDone, 'micrographs_even.sqlite',
                                self._getOutputMicEvenName,
                                OUT_MICS_EVEN, streamMode)
            self._updateOutputMicSet(newDone, 'micrographs_odd.sqlite',
                                self._getOutputMicOddName,
                                OUT_MICS_ODD, streamMode)

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return

        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m for m in self.listOfMovies
                   if m.getObjId() not in doneList and self._isMovieDone(m)]

        # Update the file with the newly done movies
        # or exit from the function if no new done movies
        self.debug('_checkNewOutput: ')
        self.debug('   listOfMovies: %s, doneList: %s, newDone: %s'
                   % (len(self.listOfMovies), len(doneList), len(newDone)))

        self._firstTimeOutput = len(doneList) == 0
        allDone = len(doneList) + len(newDone)
        # We have finished when there is not more input movies (stream closed)
        # and the number of processed movies is equal to the number of inputs
        self.finished = self.streamClosed and allDone == len(self.listOfMovies)
        streamMode = pwobj.Set.STREAM_CLOSED if self.finished else pwobj.Set.STREAM_OPEN

        if newDone:
            self._writeDoneList(newDone)

        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        self.debug('   finished: %s ' % self.finished)
        self.debug('        self.streamClosed (%s) AND' % self.streamClosed)
        self.debug('        allDone (%s) == len(self.listOfMovies (%s)'
                   % (allDone, len(self.listOfMovies)))
        self.debug('   streamMode: %s' % streamMode)

        self._updateOutputSets(newDone, streamMode)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(pwcts.STATUS_NEW)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        # Only validate about cropDimensions if the protocol supports them
        if (hasattr(self, 'cropDimX') and hasattr(self, 'cropDimY')
                and (self.cropDimX > 0 >= self.cropDimY
                     or self.cropDimX <= 0 < self.cropDimY)):
            errors.append("If you give cropDimX, you should also give "
                          "cropDimY and vice versa")

        inputMovies = self.inputMovies.get()

        # Do not continue if there ar no movies. Validation message will
        # take place since attribute is a Pointer.
        if inputMovies is None:
            return errors

        firstItem = inputMovies.getFirstItem()

        firstFrame, lastFrame, _ = inputMovies.getFramesRange()
        if lastFrame == 0:
            # Although getFirstItem is not recommended in general, here it is
            # used only once, for validation purposes, so performance
            # problems should not appear.
            frames = firstItem.getNumberOfFrames()
            lastFrame = frames
        else:
            frames = lastFrame - firstFrame + 1

        if frames is not None:
            def _validateRange(prefix):
                # Avoid validation when the range is not defined
                if not hasattr(self, '%sFrame0' % prefix):
                    return

                f0, fN = self._getFrameRange(frames, prefix)
                if fN < firstFrame or fN > lastFrame:
                    errors.append("Check the selected last frame to *%s*. "
                                  "Last frame (%d) should be in range: %s "
                                  % (prefix.upper(), fN, (firstFrame,
                                                          lastFrame)))
                if f0 < firstFrame or f0 > lastFrame:
                    errors.append("Check the selected first frame to *%s*. "
                                  "First frame (%d) should be in range: %s "
                                  % (prefix.upper(), f0, (firstFrame,
                                                          lastFrame)))
                if fN < f0:
                    errors.append("Check the selected frames range to *%s*. "
                                  "Last frame (%d) should be greater or equal "
                                  "than first frame (%d)"
                                  % (prefix.upper(), fN, f0))

            _validateRange("align")
            _validateRange("sum")

        if not emlib.image.ImageHandler().existsLocation(firstItem.getLocation()):
            errors.append("The input movie files do not exist!!! "
                          "Since usually input movie files are symbolic links, "
                          "please check that links are not broken if you "
                          "moved the project folder. ")

        return errors

    # --------------------------- INFO functions -------------------------------
    def _summary(self):
        return [self.summaryVar.get('')]

    # --------------------------- UTILS functions ----------------------------
    def _doSplitEvenOdd(self):
        """ Returns if even/odd stuff has to be done"""
        if not self.evenOddCapable:
            return False
        else:
            return self.splitEvenOdd.get()

    def _useAlignToSum(self):
        return self.getAttributeValue('useAlignToSum', False)

    def _getFrameRange(self, n, prefix):
        """
        Params:
        :param n: Number of frames of the movies
        :param prefix: what range we want to consider, either 'align' or 'sum'
        :return: (i, f) initial and last frame range
        """
        # In case that the user select the same range for ALIGN and SUM
        # we also use the 'align' prefix
        if self._useAlignToSum():
            prefix = 'align'

        first = self.getAttributeValue('%sFrame0' % prefix)
        last = self.getAttributeValue('%sFrameN' % prefix)

        if first <= 1:
            first = 1

        if last <= 0:
            last = n

        return first, last

    def _createOutputMovie(self, movie):
        # Parse the alignment parameters and store the log files
        alignedMovie = movie.clone()
        n = movie.getNumberOfFrames()
        first, last = self._getFrameRange(n, 'align')
        framesRange = alignedMovie.getFramesRange()
        framesRange.setFirstFrame(first)
        framesRange.setLastFrame(last)
        # Check if user selected to save movie, use the getAttributeValue
        # function for allow the protocol to not define this flag
        # and use False as default
        if self.getAttributeValue('doSaveMovie', False):
            # The subclass protocol is responsible for generating the output
            # movie file in the extra path with the required name
            extraMovieFn = self._getExtraPath(self._getOutputMovieName(movie))
            alignedMovie.setFileName(extraMovieFn)
            # When the output movies are saved, the shifts
            # will be set to zero since they are aligned
            totalFrames = last - first + 1
            xshifts = [0] * totalFrames
            yshifts = xshifts
            # If we save the movies, we need to modify which are the index
            # of the first frame in the stack, now is 1 since the stack is
            # written only with the given frames
            firstFrameIndex = 1
        else:
            xshifts, yshifts = self._getMovieShifts(movie)
            firstFrameIndex = first

        framesRange.setFirstFrameIndex(firstFrameIndex)
        alignment = emobj.MovieAlignment(first=first, last=last, xshifts=xshifts,
                                         yshifts=yshifts)

        roiList = [self.getAttributeValue(s, 0) for s in
                   ['cropOffsetX', 'cropOffsetY', 'cropDimX', 'cropDimY']]
        alignment.setRoi(roiList)
        alignedMovie.setAlignment(alignment)

        return alignedMovie

    # ---------- Hook functions that need to be implemented in subclasses ------
    def _getBinFactor(self):
        return self.getAttributeValue('binFactor', 1.0)

    def _getMovieRoot(self, movie):
        # Try to use the 'original' fileName in case it is present
        # the original could be different from the current filename if
        # we are dealing with compressed movies (e.g., movie.mrc.bz2)
        fn = movie.getAttributeValue('_originalFileName',
                                     movie.getFileName())
        # Remove the first extension
        fnRoot = pwutils.removeBaseExt(fn)
        # Check if there is a second extension
        # (Assuming it is only a dot and 3 or 4 characters after it
        # Do not perform this check if the file name is short
        if len(fnRoot) > 5:
            if fnRoot[-4] == '.' or fnRoot[-5] == '.':
                fnRoot = pwutils.removeExt(fnRoot)

        return fnRoot

    def _getOutputMovieName(self, movie):
        """ Returns the name of the output movie.
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_movie.mrcs'

    def _getOutputMicName(self, movie):
        """ Returns the name of the output micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_mic.mrc'

    def _getOutputMicWtName(self, movie):
        """ Returns the name of the output dose-weighted micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_mic_DW.mrc'

    def _getOutputMicEvenName(self, movie):
        """ Returns the name of the output EVEN micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_mic_EVN.mrc'

    def _getOutputMicOddName(self, movie):
        """ Returns the name of the output EVEN micrograph
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_mic_ODD.mrc'

    def _getOutputMicThumbnail(self, movie):
        return self._getExtraPath(self._getMovieRoot(movie) + '_thumbnail.png')

    def _getMovieShifts(self, movie):
        """ Returns the x and y shifts for the alignment of this movie.
         The shifts should refer to the original micrograph without any binning.
         In case of a binning greater than 1, the shifts should be scaled.
        """
        return [], []

    def _createOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return True

    def _createOutputMicrographs(self):
        """ By default check if the user have selected 'doSaveAveMic'
        property. Subclasses can override this method to implement different
        behaviour.
        """
        return self.getAttributeValue('doSaveAveMic', True)

    def _createOutputWeightedMicrographs(self):
        return False

    def _preprocessOutputMicrograph(self, mic, movie):
        """ Hook function that will be call before adding the micrograph
        to the output set of micrographs.
        """
        pass

    def _doComputeMicThumbnail(self):
        """ Should be implemented in sub-classes if want to check
        the generation of thumbnails.
        """
        return False

    def _storeSummary(self, movie):
        """ Implement this method if you want to store the summary. """
        pass

    def _getCorrectedDose(self, movieSet):
        """get and correct the pre-exposure dose. It is important for cases
        in which the first frame is different of one. The method support both
        movie and sets of movies"""

        firstFrame, _, _ = movieSet.getFramesRange()
        preExp = movieSet.getAcquisition().getDoseInitial()
        dose = movieSet.getAcquisition().getDosePerFrame()
        preExp += dose * (firstFrame - 1)

        return preExp, dose

    def __runXmippProgram(self, program, args):
        """ Internal shortcut function to launch a Xmipp program. """
        from pwem import Domain
        xmipp3 = Domain.importFromPlugin('xmipp3')
        xmipp3.Plugin.runXmippProgram(program, args)

    def __runEman2Program(self, program, args):
        """ Internal workaround to launch an EMAN2 program. """
        from pwem import Domain
        eman2 = Domain.importFromPlugin('eman2')
        from pyworkflow.utils.process import runJob
        runJob(self._log, eman2.Plugin.getProgram(program), args,
               env=eman2.Plugin.getEnviron())

    def computePSD(self, inputMic, oroot, dim=384,  # 384 = 128 + 256, which should be fast for any Fourier Transformer
                   overlap=0.4):
        warnings.warn("Use psd = image.computePSD(overlap=0.4, xdim=384, ydim=384, fftthreads=1) instead",
                      DeprecationWarning)
        ih = emlib.image.ImageHandler()
        psdImg1 = ih.read(inputMic)
        res = psdImg1.computePSD(overlap, dim, dim)
        res.write(oroot + ".psd")

    def composePSDImages(self, psdImg1, psdImg2, outputFn,
                         outputFnUncorrected=None, outputFnCorrected=None):
        """ Compose a single PSD image:
         left part from psd1 (uncorrected PSD),
         right-part from psd2 (corrected PSD)
        """
        data1 = psdImg1.getData()  # get the data now, as conversion would change them
        if outputFnUncorrected is not None:
            psdImg1.convertPSD()
            psdImg1.write(outputFnUncorrected)

        data2 = psdImg2.getData()  # get the data now, as conversion would change them
        if outputFnCorrected is not None:
            psdImg2.convertPSD()
            psdImg2.write(outputFnCorrected)

        # Compute middle index
        x, _, _, _ = psdImg1.getDimensions()
        m = int(round(x / 2.))
        data1[:, :m] = data2[:, :m]
        psdImg1.setData(data1)
        psdImg1.write(outputFn)

    def composePSD(self, psd1, psd2, outputFn,
                   outputFnUncorrected=None, outputFnCorrected=None):
        import warnings
        warnings.warn("Use composePSDImages() instead", DeprecationWarning)
        """ Compose a single PSD image:
         left part from psd1 (uncorrected PSD),
         right-part from psd2 (corrected PSD)
        """
        ih = emlib.image.ImageHandler()
        self.composePSDImages(ih.read(psd1), ih.read(psd2), outputFn,
                              outputFnUncorrected, outputFnCorrected)

    def computePSDImages(self, movie, fnUncorrected, fnCorrected,
                         outputFnUncorrected=None, outputFnCorrected=None):
        self.composePSDImages(
            emlib.image.ImageHandler().read(fnUncorrected).computePSD(),
            emlib.image.ImageHandler().read(fnCorrected).computePSD(),
            self._getPsdCorr(movie),
            outputFnUncorrected,
            outputFnCorrected)

    def computePSDs(self, movie, fnUncorrected, fnCorrected,
                    outputFnUncorrected=None, outputFnCorrected=None):
        import warnings
        warnings.warn("Use computePSDImages() instead", DeprecationWarning)
        self.computePSDImages(movie, fnUncorrected, fnCorrected,
                              outputFnUncorrected, outputFnCorrected)

    def computeThumbnail(self, inputFn, scaleFactor=6, outputFn=None):
        """ Generates a thumbnail of the input file"""
        outputFn = outputFn or self.getThumbnailFn(inputFn)
        args = "%s %s " % (inputFn, outputFn)
        args += "--meanshrink %s --fixintscaling=sane" % scaleFactor

        self.__runEman2Program('e2proc2d.py', args)

        return outputFn

    def correctGain(self, movieFn, outputFn, gainFn=None, darkFn=None):
        """correct a movie with both gain and dark images"""
        ih = emlib.image.ImageHandler()
        _, _, z, n = ih.getDimensions(movieFn)
        numberOfFrames = max(z, n)  # in case of wrong mrc stacks as volumes

        def _readImgFloat(fn):
            img = None
            if fn:
                img = ih.read(fn)
                img.convert2DataType(emlib.DT_FLOAT)
            return img

        gainImg = _readImgFloat(gainFn)
        darkImg = _readImgFloat(darkFn)

        img = ih.createImage()

        for i in range(1, numberOfFrames + 1):
            img.read((i, movieFn))
            img.convert2DataType(emlib.DT_FLOAT)

            if darkImg:
                img.inplaceSubtract(darkImg)
            if gainImg:
                img.inplaceMultiply(gainImg)

            img.write((i, outputFn))

    def getThumbnailFn(self, inputFn):
        """ Returns the default name for a thumbnail image"""
        return pwutils.replaceExt(inputFn, "thumb.png")

    def _getPsdCorr(self, movie):
        """ This should be implemented in subclasses."""
        pass

    def _writeFailedList(self, movieList):
        """ Write to a text file the items that have failed. """
        with open(self._getAllFailed(), 'a') as f:
            for movie in movieList:
                f.write('%d\n' % movie.getObjId())

    def _readFailedList(self):
        """ Read from a text file the id's of the items that have failed. """
        failedFile = self._getAllFailed()
        failedList = []
        if os.path.exists(failedFile):
            with open(failedFile) as f:
                failedList += [int(line.strip()) for line in f]

        return failedList


def createAlignmentPlot(meanX, meanY):
    """ Create a plotter with the cumulative shift per frame. """
    figureSize = (8, 6)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()

    ax = figure.add_subplot(111)
    ax.grid()
    ax.axis('equal')
    ax.set_title('Cartesian representation')
    ax.set_xlabel('Drift x (pixels)')
    ax.set_ylabel('Drift y (pixels)')

    # Max range of the plot of the two coordinates
    plotRange = max(max(meanX) - min(meanX), max(meanY) - min(meanY))
    i = 1
    skipLabels = ceil(len(meanX) / 10.0)
    for x, y in zip(meanX, meanY):
        if i % skipLabels == 0:
            ax.text(x - 0.02 * plotRange, y + 0.02 * plotRange, str(i))
        i += 1

    ax.plot(meanX, meanY, color='b')
    ax.plot(meanX, meanY, 'yo')

    # setting the plot windows to properly see the data
    ax.axis([min(meanX) - 0.1 * plotRange, max(meanX) + 0.1 * plotRange,
             min(meanY) - 0.1 * plotRange, max(meanY) + 0.1 * plotRange])

    plotter.tightLayout()

    return plotter


class ProtAverageFrames(ProtAlignMovies):
    """
    Very simple protocol to align all the frames of a given data collection
    session. It can be used as a sanity check.
    """
    _label = 'average frames'

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):
        pass

    def _processMovie(self, movie):
        allFramesSum = self._getPath('all_frames_sum.mrc')
        ih = emlib.image.ImageHandler()
        sumImg = ih.createImage()
        img = ih.createImage()

        n = movie.getNumberOfFrames()
        fn = movie.getFileName()

        sumImg.read((1, fn))

        for frame in range(2, n + 1):
            img.read((frame, fn))
            sumImg.inplaceAdd(img)

        if os.path.exists(allFramesSum):
            img.read(allFramesSum)
            sumImg.inplaceAdd(img)

        sumImg.write(allFramesSum)

    # FIXME: Methods will change when using the streaming for the output
    def createOutputStep(self):
        # Really load the input, since in the streaming case we can not
        # use the self.inputMovies directly
        allFramesSum = self._getPath('all_frames_sum.mrc')
        allFramesAvg = self._getPath('all_frames_avg.mrc')
        self._loadInputList()
        n = len(self.listOfMovies)

        ih = emlib.image.ImageHandler()
        sumImg = ih.read(allFramesSum)
        sumImg.inplaceDivide(float(n))
        sumImg.write(allFramesAvg)

        outputAvg = emobj.Image()
        outputAvg.setFileName(allFramesAvg)
        outputAvg.setSamplingRate(self.listOfMovies[0].getSamplingRate())
        self._defineOutputs(outputAverage=outputAvg)
        self._defineSourceRelation(self.inputMovies, outputAvg)

    def _validate(self):
        return []

    def _summary(self):
        return []

    def _createOutputMovies(self):
        """ Returns True if an output set of movies will be generated.
        The most common case is to always generate output movies,
        either with alignment only or the binary aligned movie files.
        Subclasses can override this function to change this behavior.
        """
        return False

    def _createOutputMicrographs(self):
        """ By default check if the user have selected 'doSaveAveMic'
        property. Subclasses can override this method to implement different
        behaviour.
        """
        return False

    def _createOutputWeightedMicrographs(self):
        return False
