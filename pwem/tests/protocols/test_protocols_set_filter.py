#!/usr/bin/env python
# ***************************************************************************
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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


import numpy as np
from PIL import Image

import pyworkflow.tests as pwtests
import pwem.objects as emobj
from pwem.protocols import ProtSetFilter, EMProtocol

OUTPUT_COORDINATES = "outputCoordinates"

class TestSetFilter(pwtests.BaseTest):
    """Run different tests related to the editor set protocol."""
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)


    def createCoorSetProtocol(self):
        MICNUMBER = 2
        MICSIZE = 1000

        # create dummy protocol
        dummyProt = self.newProtocol(EMProtocol)
        dummyProt.setObjLabel('dummy protocol')
        self.launchProtocol(dummyProt)

        # create set of micrographs
        #   create two tif image
        matrix = np.zeros((MICSIZE,MICSIZE), np.uint8)
        R = np.linspace(0, MICSIZE // 2 -2, 20)
        # note we have 0 and 2*pi for each radius
        THETA = np.linspace(0, 2* np.pi, 30)
        radii, thetas = np.meshgrid(R, THETA)
        X = (R * np.cos(thetas)).astype(int) + MICSIZE // 2
        Y = (R * np.sin(thetas)).astype(int) + MICSIZE // 2
        # TODO: next double loop should be equivalen to
        # the more elegant solution
        # pairlist = np.vstack(list(map(np.ravel, (X,Y)))).T
        # matrix[pairList] = 255 but it is not

        for xx, yy in zip(X,Y):
            for x, y in zip(xx, yy):
                matrix[y][x] = 255

        fn = dummyProt._getExtraPath('mic_000001.tif')
        Image.fromarray(matrix).save(fn)
        fn = dummyProt._getExtraPath('mic_000002.tif')
        Image.fromarray(matrix).save(fn)
        micSet = dummyProt._createSetOfMicrographs()
        mic = emobj.Micrograph()

        for i in range(MICNUMBER):
            mic.setFileName\
                (dummyProt._getExtraPath("mic_%06d.tif" % (i%2 + 1 )))
            mic.setMicName("mic_%06d.tif" % (i%2 + 1 ))
            mic.setSamplingRate(1.)
            micSet.append(mic)
            mic.cleanObjId()
        micSet.setSamplingRate(1.)
        micSet.write()

        coordSet = dummyProt._createSetOfCoordinates(micSet)
        coordSet.setBoxSize(10)
        coord = emobj.Coordinate()

        for xx, yy in zip(X,Y):
            for x, y in zip(xx, yy):
                for mic in range(1, MICNUMBER + 1):
                    coord.setX(x)
                    coord.setY(y)
                    coord.setMicId(mic)
                    coord.setMicName("mic_%06d.tif" % (mic))
                    coordSet.append(coord)
                    coord.cleanObjId()
        coordSet.write()
        outputArgs = {OUTPUT_COORDINATES: coordSet,
                      'outputMic': micSet}
        dummyProt._defineOutputs(**outputArgs)
        dummyProt._store()
        return dummyProt

    def testOperation(self):
        """Make a trivial operation. keep coordinates with xcoor > 200"""
        dummyProt = self.createCoorSetProtocol()
        protSetFilter = self.newProtocol(ProtSetFilter,
                                          objLabel="operate")
        protSetFilter.inputSet.set(dummyProt)
        protSetFilter.inputSet.setExtended(OUTPUT_COORDINATES)
        protSetFilter.operation.set(protSetFilter.CHOICE_FORMULA)
        protSetFilter.formula.set('item._x.get()>200')
        self.launchProtocol(protSetFilter)
        for item in protSetFilter.outputCoordinates:
            self.assertGreater(item.getX(), 100)
        self.assertEqual(len(dummyProt.outputCoordinates), 1200)
        self.assertEqual(len(protSetFilter.outputCoordinates), 1100)

    def testCenter(self):
        """Remove coordinates closer to center less than 200"""
        dummyProt = self.createCoorSetProtocol()
        protSetFilter = self.newProtocol(ProtSetFilter,
                                         objLabel="center")
        protSetFilter.inputSet.set(dummyProt)
        protSetFilter.inputSet.setExtended(OUTPUT_COORDINATES)
        protSetFilter.operation.set(protSetFilter.CHOICE_DISTANCE_CENTER)
        protSetFilter.distance.set(200)
        self.launchProtocol(protSetFilter)
        self.assertEqual(len(protSetFilter.outputCoordinates), 720)
        self.assertEqual(len(dummyProt.outputCoordinates), 1200)

    def testDistance(self):
        """Remove coordinates closer to center less than 200"""
        dummyProt = self.createCoorSetProtocol()
        protSetFilter = self.newProtocol(ProtSetFilter,
                                         objLabel="distance- keepfirst true")
        protSetFilter.inputSet.set(dummyProt)
        protSetFilter.inputSet.setExtended(OUTPUT_COORDINATES)
        protSetFilter.operation.set(protSetFilter.CHOICE_DISTANCE_BETWEEN_COORDS)
        protSetFilter.distance.set(10)
        protSetFilter.keepFirst.set(True)
        self.launchProtocol(protSetFilter)
        self.assertEqual(len(protSetFilter.outputCoordinates), 1048)
        self.assertEqual(len(dummyProt.outputCoordinates), 1200)
        protSetFilter = self.newProtocol(ProtSetFilter,
                                         objLabel="distance- keepfirst false")
        protSetFilter.inputSet.set(dummyProt)
        protSetFilter.inputSet.setExtended(OUTPUT_COORDINATES)
        protSetFilter.operation.set(protSetFilter.CHOICE_DISTANCE_BETWEEN_COORDS)
        protSetFilter.distance.set(10)
        protSetFilter.keepFirst.set(False)
        self.launchProtocol(protSetFilter)
        self.assertEqual(len(protSetFilter.outputCoordinates), 1008)
        self.assertEqual(len(dummyProt.outputCoordinates), 1200)

    def testRanking(self):
        """Test ranking filtering"""
        dummyProt = self.createCoorSetProtocol()
        self.addFilterSetProtocol(dummyProt, '10%', dummyProt.outputCoordinates.getSize()*0.1)
        self.addFilterSetProtocol(dummyProt, '-0.3', dummyProt.outputCoordinates.getSize() * 0.3)

        self.addFilterSetProtocol(dummyProt, '10', 10)
        self.addFilterSetProtocol(dummyProt, '-20', 20)

    def addFilterSetProtocol(self, inputProt, threshold, expectedSize):

        protSetFilter = self.newProtocol(ProtSetFilter,
                                          objLabel="ranking %s" % threshold)
        protSetFilter.inputSet.set(inputProt)
        protSetFilter.inputSet.setExtended(OUTPUT_COORDINATES)
        protSetFilter.operation.set(protSetFilter.CHOICE_RANKED)
        protSetFilter.threshold.set(threshold)
        protSetFilter.rankingField.set('_x')
        self.launchProtocol(protSetFilter)
        self.assertSetSize(protSetFilter.outputCoordinates, expectedSize)