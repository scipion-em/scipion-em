# ***************************************************************************
# * Authors:     Roberto Marabini
# *
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
# ***************************************************************************/
from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_movie_eraser import ProtMovieEraser
from pyworkflow.tests import Manager, logger
import os
import pwem.objects as emobj
import numpy as np
from pwem.protocols import EMProtocol
from PIL import Image


class TestProtMovieEraser(BaseTest):
    """ Test wait protocol"""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def createSetOfCoordiantes(self, noCoordinatesMicList):
        # store number of coordenates per movie
        MICSIZE = 1000  # mic size
        # matrix with pixel values for the microgrphs
        matrix = np.zeros((MICSIZE,MICSIZE), np.uint8) 

        movieList=[]
        for i, noCoordinatesMic in enumerate(noCoordinatesMicList):
            movieList.append({
                "coordNo": noCoordinatesMic,
                "index": i+1})
        # create fake protocol so we can store the set of 
        # movies and the set of coordiantes
        # create dummy protocol
        prot = self.newProtocol(EMProtocol)
        prot.setObjLabel('dummy protocol')
        self.launchProtocol(prot)

        # create temporary set of micrographs
        # since setOfCoordinates requires a setOfMicrographs
        micSet = prot._createSetOfMicrographs()

        mic = emobj.Micrograph()
        MICNUMBER = len(movieList)
        for movie in movieList:
            if movie["coordNo"] == 0:
                continue
            i = movie["index"]
            fn = prot._getExtraPath("mic_%06d.tif" % (i))
            Image.fromarray(matrix).save(fn)
            mic.setFileName(fn)
            mic.setMicName("mic_%06d.tif" % (i))
            mic.setSamplingRate(1.)
            micSet.append(mic)
            mic.cleanObjId()
        micSet.setSamplingRate(1.)
        micSet.write()

        # create set of coordinates
        coordSet = prot._createSetOfCoordinates(micSet)
        coordSet.setBoxSize(10)
        coord = emobj.Coordinate()
        for movie in movieList:
            for j in range(movie["coordNo"]):
                i = movie["index"]
                coord.setX(j)
                coord.setY(j)
                coord.setMicId(i)
                coord.setMicName("mic_%06d.tif" % (i))
                coordSet.append(coord)
                coord.cleanObjId()
        coordSet.write()

        # create set of movies
        movSet = prot._createSetOfMovies()

        mov = emobj.Movie()
        for movie in movieList:
            i = movie["index"]
            fn = prot._getExtraPath("mov_%06d.tif" % (i))
            Image.fromarray(matrix).save(fn)
            mov.setFileName(fn)
            mov.setMicName("mic_%06d.tif" % (i))
            mov.setSamplingRate(1.)
            movSet.append(mov)
            mov.cleanObjId()
        movSet.setSamplingRate(1.)
        movSet.write()


        outputArgs = {'outputCoord': coordSet,
                      'outputMic': micSet,
                      'outputMov': movSet}
        prot._defineOutputs(**outputArgs)
        prot._store()
        return prot
    
    def testCoord(self):
        noCoordinatesPerMovie=[1, 2, 0, 4, 0]
        dummyProtocol = self.createSetOfCoordiantes(noCoordinatesPerMovie)
        for mov in dummyProtocol.outputMov:
            fn = mov.getFileName()
            self.assertTrue(os.path.exists(fn))
        protMovieEraser = self.newProtocol(ProtMovieEraser,
                                    sourceSet=dummyProtocol.outputCoord,
                                    movieSet=dummyProtocol.outputMov,
                                    dryMode=False
                                   )

        self.launchProtocol(protMovieEraser)
        for i, mov in enumerate(dummyProtocol.outputMov):
            fn = mov.getFileName()
            if i in [2, 4]:
                self.assertFalse(os.path.exists(fn))
            else:
                self.assertTrue(os.path.exists(fn))

    def testMic(self):
        noCoordinatesPerMovie=[1, 2, 0, 4, 0]
        dummyProtocol = self.createSetOfCoordiantes(noCoordinatesPerMovie)
        for mov in dummyProtocol.outputMov:
            fn = mov.getFileName()
            self.assertTrue(os.path.exists(fn))

        protMovieEraser = self.newProtocol(ProtMovieEraser,
                                    sourceSet=dummyProtocol.outputMic,
                                    movieSet=dummyProtocol.outputMov,
                                    dryMode=False
                                   )
        self.launchProtocol(protMovieEraser)
        for i, mov in enumerate(dummyProtocol.outputMov):
            fn = mov.getFileName()
            if i in [2, 4]:
                self.assertTrue(os.path.exists(fn))
            else:
                self.assertFalse(os.path.exists(fn))
                
        
