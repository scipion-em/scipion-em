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
from math import sqrt
from PIL import Image

import pyworkflow.tests as pwtests
import pwem.objects as emobj
from pwem.protocols import EMProtocol
from pwem.objects.data import Particle, SetOfParticles, Acquisition, CTFModel
from pwem.constants import ALIGN_PROJ
from pwem.convert.trigonometry import FibonacciSphere
from pwem.convert.transformations import rotation_matrix

class TestAngularDistribution(pwtests.BaseTest):
    """Run different tests related to the editor set protocol."""
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.samplingrate = 1.5
        cls.numParticles = 10000
        cls.fibonacci = FibonacciSphere(cls.numParticles)
        # create dummy protocol
        cls.dummyProt = cls.newProtocol(EMProtocol)
        cls.dummyProt.setObjLabel('distribution_1')
        cls.launchProtocol(cls.dummyProt)
        cls.partFn = cls.dummyProt._getExtraPath('particles.stk')

    def _setUp(self, suffix=None) -> None:
        # create set of particles
        self.partSet = self.dummyProt._createSetOfParticles(suffix=suffix)
        self.partSet.setAlignment(ALIGN_PROJ)
        self.partSet.setAcquisition(Acquisition(voltage=300,
                                           sphericalAberration=2,
                                           amplitudeContrast=0.1,
                                           magnification=60000))
        self.partSet.setSamplingRate(self.samplingrate)
    
    def tearDown(self) -> None:
        self.dummyProt._defineOutputs(**self.outputArgs)
        self.dummyProt._store()
        return super().tearDown()
    
    def createParticles(self, points):
        matrices = []
        for point in points:
            if point[2] == 1 or point[2] == -1:
                rot_axis = [1, 0, 0]
            else:
                z_axis = [0, 0, 1]
                rot_axis = np.cross(z_axis, point)
            angle = np.arccos(np.dot(point, z_axis))
            matrices.append(rotation_matrix(angle, rot_axis))

        for i, matrix in enumerate(matrices):
            p = Particle()
            p.setLocation(i + 1, self.partFn)
            p.setTransform(emobj.Transform(matrix))
            self.partSet.append(p)
            
    def populateParticles_few_points(self, partSet, partFn):
        print("""This test picks a few points points in x, y and z axis
              plus (sqrt(1/3), sqrt(1/3), sqrt(1/3)).
              Converts them to transformation matrices.
              Then, two particles are generated with those matrices.
              """)
        points= []
        points.append([1, 0, 0])
        points.append([0, 1, 0])
        points.append([0, 0, 1])
        coord = sqrt(1/3.)
        points.append([coord, coord, coord])
        self.createParticles(points)
        
    def populateParticles_up_down(self, partSet, partFn):
        print("""This test picks  points with abs(z) > 0.5.
              distribution should be almost uniform.
              """)
        points= []
        X, Y , Z = self.fibonacci.sphX, self.fibonacci.sphY, self.fibonacci.sphZ
        for x, y, z in (zip(X, Y, Z)):
            if abs(z) > 0.5:
                points.append([x, y, z])
        self.createParticles(points)
            
    def populateParticles_Z(self, partSet, partFn):
        print("""This test creates points  proportional to |z|.
              """)
        points= []
        X, Y , Z = self.fibonacci.sphX, self.fibonacci.sphY, self.fibonacci.sphZ
        for x, y, z in (zip(X, Y, Z)):
            zz = abs(round(z*10)) + 1
            for i in range(zz):
                points.append([x, y, z])
        self.createParticles(points)

    def populateParticles_45(self, partSet, partFn):
        print("""This test creates points  in 1/8 of the sphere.
              """)
        points= []
        X, Y , Z = self.fibonacci.sphX, self.fibonacci.sphY, self.fibonacci.sphZ
        for x, y, z in (zip(X, Y, Z)):
            if (x > 0) and (y > 0) and (x < y):
                points.append([x, y, z])
        self.createParticles(points)
    
    def test_angularDistribution_few_points(self):
        self._setUp(suffix='few_points')
        self.populateParticles_few_points(self.partSet, self.partFn)
        self.partSet.write()
        self.outputArgs = {"outputTwoParticles": self.partSet}
        
    def test_angularDistribution_up_down(self):
        self._setUp(suffix='up_down')
        self.populateParticles_up_down(self.partSet, self.partFn)
        self.partSet.write()
        self.outputArgs = {"outputUpDownParticles": self.partSet}
        
    def test_angularDistribution_Z(self):
        self._setUp(suffix='Z')
        self.populateParticles_Z(self.partSet, self.partFn)
        self.partSet.write()
        self.outputArgs = {"outputUpZ": self.partSet}

    def test_angularDistribution_45(self):
        self._setUp(suffix='45')
        self.populateParticles_45(self.partSet, self.partFn)
        self.partSet.write()
        self.outputArgs = {"output45": self.partSet}
