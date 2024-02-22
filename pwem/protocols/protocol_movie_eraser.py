import os
import time

from pwem.objects import SetOfMovies, SetOfCoordinates, SetOfMicrographs
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import params
import pyworkflow.utils as pwutils
from pwem.objects.data import Movie
from copy import deepcopy
from pyworkflow.utils import ProgressBar

class ProtMovieEraser(EMProtocol):
    """ 
    Protocol for removing movies based on different conditions:
    - If the input is SetOfMicrographs, it removes movies already aligned into micrographs.
    - If the input is SetOfCoordinates, it removes movies with unselected particles.
    WARNING: There is no way back. Be sure you understand the consequences.
    """
    _label = "movie eraser"
    _devStatus = BETA
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('sourceSet', params.PointerParam,
                      pointerClass='SetOfMicrographs, SetOfCoordinates',
                      label="Input set", important=True,
                      help='Select the SetOfMicrographs or Coordinates to drive the movie DELETION process.')
        # TODO: ADD SetOFCoordinates
        form.addParam('movieSet', params.PointerParam,
                      pointerClass='SetOfMovies',
                      label="Movies source", important=True,
                      help='Select the SetOfMovies to locate the movie files. Match will be done by ids!.')

        form.addParam("dryMode", params.BooleanParam, label="Dry mode",
                      default=True, expertLevel=params.LEVEL_ADVANCED)

    # -------------------------- STEP methods ---------------------------------
    def _insertAllSteps(self):
        # Insert different steps based on the input type
        sourceSet = self.sourceSet.get()
        if isinstance(sourceSet, SetOfMicrographs):
            self._insertFunctionStep(self.deleteMoviesMicStep)
            state = False
        elif isinstance(sourceSet, SetOfCoordinates):
            self._insertFunctionStep(self.deleteMoviesCoorStep)
            state = True
        else:
            raise NotImplementedError(
                "Only SetOfMicrographs or SetOfCoordinates are supported")
        self._insertFunctionStep(self.deleteMoviesStep)
        
        # Initialize a dictionary to keep track of movies and
        # their deletion status
        self.movMicNameList = {}
        self.movieSet = self.movieSet.get()
        self.coordSet = self.sourceSet.get()
        for mov in self.movieSet:
            newMovie = Movie()
            newMovie = deepcopy(mov)
            self.movMicNameList[mov.getMicName()] =\
               {'mov': newMovie,
                'deletethis' : state}

    def deleteMoviesMicStep(self):
        """ Delete all movies origin of the micset until set is closed"""
        lastMicId = 0
        # create an auxiliary list with micNames so we
        # do not need to make a search in the databse for each mic

        # For each mic 
        micSet = self.sourceSet.get()
        for mic in micSet:
            self.movMicNameList[mic.getMicName()]['deletethis'] = True

    def deleteMoviesCoorStep(self):
        """ Delete movies that are not related with the coordinate set"""
        movieSetWithCoordinates = self.coordSet.getUniqueValues("_micName")

        for micName in movieSetWithCoordinates:
            self.movMicNameList[micName]['deletethis'] = False

    def deleteMoviesStep(self):
        """ Delete movies using the movMicNameList"""
        deletedMovies = self._createSetOfMovies(suffix="deleted")
        keptMovies = self._createSetOfMovies(suffix="kept")
        deletedMovies.copyInfo(self.movieSet )
        keptMovies.copyInfo(self.movieSet )
        try: ## PHN: for security
            deletedMovies.setGain(self.movieSet.getGain())
            deletedMovies.setDark(self.movieSet.getDark())
            deletedMovies.setFramesRange(self.movieSet.getFramesRange())
            deletedMovies._firstFramesRange = self.movieSet._firstFramesRange
            keptMovies.setGain(self.movieSet.getGain())
            keptMovies.setDark(self.movieSet.getDark())
            keptMovies.setFramesRange(self.movieSet.getFramesRange())
            keptMovies._firstFramesRange = self.movieSet._firstFramesRange
        except Exception as e:
            raise AttributeError(f"Missing parameter {e}")

        # Iterate through movies, delete and update sets
        progress = ProgressBar(total=len(self.movieSet), fmt=ProgressBar.NOBAR)
        progress.start()
        step = max(100, len(self.movieSet) // 100)

        for i, mov in enumerate(self.movMicNameList.values()):
            if i % step == 0:
                progress.update(i+1)
            newMovie = Movie()
            oldMovie = mov['mov']
            newMovie.copyInfo(oldMovie)
            if mov['deletethis']:
                deletedMovies.append(newMovie)
                self._deleteMovie(oldMovie.getFileName())
            else:
                keptMovies.append(newMovie)

        progress.finish()

        outputArgs = {"deletedMovies": deletedMovies,
                      "keptMovies": keptMovies}
        self._defineOutputs(**outputArgs)    
        self._defineSourceRelation(self.movieSet, deletedMovies)
        self._defineSourceRelation(self.movieSet, keptMovies)
         

    def _deleteMovie(self, movieFile):
        """ Deletes the movie named movieFile"""

        # Check existence
        if not os.path.exists(movieFile):
            print("%s is gone already!." % movieFile)

        # if dry mode
        elif self.dryMode.get():
            print("%s would have been deleted" % movieFile, flush=True)

        # Delete!!!
        else:
            # Delete the movie
            if os.path.islink(movieFile):
                movieFile2 = os.path.realpath(movieFile)
            else:
                movieFile2 = movieFile

            os.remove(movieFile2)
            # print("%s deleted." % movieFile, movieFile2, flush=True)


    # -------------------------- INFO functions -------------------------------
    def _warnings(self):
        if not self.dryMode.get():
            return ["If you continue, all movies related to the set of micrographs will be DELETED!. Try to run it in \"dry mode\" to double check first."]

    def _summary(self):
        summary = []
        summary.append("Dry mode is active. Harmless" if self.dryMode.get() else "Dry mode deactivated, harm can be done!")
        return summary

