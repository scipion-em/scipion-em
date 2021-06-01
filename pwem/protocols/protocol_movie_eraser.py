import os
import time

from pwem.objects import SetOfMovies
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import params
import pyworkflow.utils as pwutils

class ProtMovieEraser(EMProtocol):
    """ It will REMOVE movies that already have been aligned into micrographs.
    WARNING: There is no way back. Be sure you understand the consequences."""
    _label = "movie eraser"
    _devStatus = BETA
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.deletedMovies = []

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label=pwutils.Message.LABEL_INPUT_MIC, important=True,
                      help='Select the SetOfMicrographs to drive the movie DELETION process.')

        form.addParam('movieSet', params.PointerParam,
                      pointerClass='SetOfMovies',
                      label="Movies source", important=True,
                      help='Select the SetOfMovies to locate the movie files. Match will be done by ids!.')

        form.addParam("dryMode", params.BooleanParam, label="Dry mode", default=True, expertLevel=params.LEVEL_ADVANCED)

    # -------------------------- STEP methods ---------------------------------
    def _insertAllSteps(self):

        self._insertFunctionStep(self.deleteMoviesStep)

    def deleteMoviesStep(self):
        """ Delete all movies origin of the micset until set is closed"""
        lastMicId = 0
        while not self._isAllDone():
            # For each mic not done yet
            for mic in self._getMicSet().iterItems(where="id > %s" % lastMicId):
                self._deleteMovieFromMic(mic)
                lastMicId = mic.getObjId()

            # wait some time
            time.sleep(30)


    def _deleteMovieFromMic(self, mic):
        """ Deletes the movie source of the mic object"""
        movieFile = self._getMovieFileFromMic(mic)

        # Check existence
        if not os.path.exists(movieFile):
            print("%s is gone already!." % movieFile)

        # if dry mode
        elif self.dryMode.get():
            print("%s would have been deleted" % movieFile, flush=True)

        # Delete!!!
        else:
            # Delete the move
            os.remove(movieFile)
            print("%s deleted." % movieFile, flush=True)

        # Annotate it
        self.deletedMovies.append(movieFile)

    def _getMovieFileFromMic(self, mic):
        """ Returns the movie path that was used in the micrograph movie alignment"""
        movieSet = self._getMovieSet()

        # To match the movie use the ID:
        # movie alignment protocols are copying the movie id to the mic --> mic.copyObjId(movie)
        movie = movieSet[mic.getObjId()]

        movieFile = movie.getFileName()

        if os.path.islink(movieFile):
            movieFile = os.path.realpath(movieFile)

        return movieFile

    def _getMovieSet(self):
        """ Returns the SetOFMovies related to the input set of micrographs"""
        # This is not possible since the relations are stored in the run.db of the input protocol
        # and there is no easy way to get them from this run
        # parents= self.getProject().getSourceParents(self._getMicSet())
        #
        # # Get the parent protocols
        # for parent in parents:
        #     # We need to update them taking the most recent data form its own run.db
        #     self.getProject()._updateProtocol(parent)
        #     for attr, output in parent.iterOutputAttributes():
        #         if isinstance(output, SetOfMovies):
        #             return output
        #
        # return None
        return self.movieSet.get()

    def _getMicSet(self):
        """ Returns te input set of micrographs"""
        mics = self.inputMicrographs.get()
        mics.loadAllProperties()
        return mics

    def _isAllDone(self):
        """ Returns True if all work is done. All movies deleted """
        mics = self._getMicSet()
        return mics.isStreamClosed() and mics.getSize() == len(self.deletedMovies)

    @classmethod
    def worksInStreaming(cls):
        return True
    # -------------------------- INFO functions -------------------------------
    def _warnings(self):
        if not self.dryMode.get():
            return ["If you continue, all movies related to the set of micrographs will be DELETED!. Try to run it in \"dry mode\" to double check first."]

    def _summary(self):
        summary = []
        summary.append("Dry mode is active. Harmless" if self.dryMode.get() else "Dry mode deactivated, harm can be done!")
        return summary

