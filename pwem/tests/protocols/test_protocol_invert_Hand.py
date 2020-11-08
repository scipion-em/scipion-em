# ***************************************************************************
# * Authors:    Roberto Marabini (roberto@cnb.csic.es)
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
from tempfile import mkstemp

import numpy as np
import os

from pwem.objects import Matrix
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import runJob
from xmipp3 import Plugin
import pwem.protocols as emprot


def createFeatFile(fd):
    # create feat file
    f = os.fdopen(fd, "w")
    f.write("""# XMIPP_STAR_1 *
# Type of feature (sph, blo, gau, Cyl, dcy, cub, ell, con)(Required)
# The operation after adding the feature to the phantom (+/=) (Required)
# The feature density (Required)
# The feature center (Required)
# The vector for special parameters of each vector (Required)
# Sphere: [radius] Blob : [radius alpha m] Gaussian : [sigma]
# Cylinder : [xradius yradius height rot tilt psi]
# DCylinder : [radius height separation rot tilt psi]
# Cube : [xdim ydim zdim rot tilt psi]
# Ellipsoid : [xradius yradius zradius rot tilt psi]
# Cone : [radius height rot tilt psi]
data_block1
 _dimensions3D  '100 100 100'
 _phantomBGDensity  0.
 _scale  1.
data_block2
loop_
 _featureType
 _featureOperation
 _featureDensity
 _featureCenter
 _featureSpecificVector
""")
    theta_max = 8 * np.pi
    nPoints = 100
    thetaStep = theta_max/nPoints/1.2
    for theta in np.arange(-theta_max/2.4, theta_max/2.4, thetaStep):
        z = theta * 3.
        x = np.sin(theta) * theta_max*1.3
        y = np.cos(theta) * theta_max*1.3
        f.write("sph = 1 '%f %f %f' '5'\n" % (x, y, z))
    f.close()


def createParamFile(fd):
    f = os.fdopen(fd, "w")
    f.write("""# XMIPP_STAR_1 *
# Projection Parameters
#_angleFile
# X and Y projection dimensions [Xdim Ydim]
# Rotation range and number of samples [Start Finish NSamples]
# Rotation angle added noise  [noise (bias)]
# Type of range for Rotation (random_group or random or even)
# Tilt range and number of samples for Tilt
# Tilt angle added noise
# Type of range for Tilt
# Psi range and number of samples
# Psi added noise
# Type of range for Psi
# Noise applied to pixels [noise (bias)]
# Noise applied to particle center coordinates [noise (bias)]

data_block1
_dimensions2D   '100 100'
_projRotRange    '0 360 4'
_projRotNoise   '0.0 0.0'
_projRotRandomness   random
_projTiltRange    '0 180 4'
_projTiltNoise   '0.0 0.0'
_projTiltRandomness   random
_projPsiRange    '0 360 4'
_projPsiNoise   '0.0 0.0'
_projPsiRandomness   random
_noisePixelLevel   '0. 0.'
""")
    f.close()


def projectPhantom(featFileName, paramFileName, particlesFileName):
    args = "-i %s -o %s" % (featFileName, particlesFileName)
    args += " --params %s" % paramFileName
    runJob(None, "xmipp_phantom_project", args, env=Plugin.getEnviron())


class TestProtInvertHand(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        (fd, cls.featFileName) = mkstemp(suffix=".feat")
        createFeatFile(fd)
        (fd, cls.paramFileName) = mkstemp(suffix=".param")
        createParamFile(fd)
        (fd, cls.stackFileName) = mkstemp(suffix=".stk")
        cls.xmdFileName = cls.stackFileName.replace(".stk", ".xmd")
        projectPhantom(cls.featFileName, cls.paramFileName,
                       cls.stackFileName)
        (fd, cls.stackFileName2) = mkstemp(suffix=".stk")
        cls.xmdFileName2 = cls.stackFileName2.replace(".stk", ".xmd")
        projectPhantom(cls.featFileName, cls.paramFileName,
                       cls.stackFileName2)

    def test_01_autoaligment(self):
        # import first set of particles
        prot1 = self.newProtocol(emprot.ProtImportParticles,
                                 importFrom=emprot.ProtImportParticles.IMPORT_FROM_XMIPP3,
                                 mdFile=self.xmdFileName,
                                 magnification=10000,
                                 samplingRate=1,
                                 haveDataBeenPhaseFlipped=False
                                 )
        prot1.setObjLabel('import Particles 1')
        self.launchProtocol(prot1)

        # assign transformation from 1 set to second set
        prot3 = self.newProtocol(emprot.ProtAlignmentInvertHand,
                                 inputParticles=prot1.outputParticles,
                                 )
        prot3.setObjLabel('assign angle 1 -> 2')
        self.launchProtocol(prot3)
        for part1, part3 in emprot.izip(
                prot1.outputParticles,
                prot3.outputParticles):
            matrix = Matrix()
            matrix.copy(part1.getTransform())
            m = matrix.getMatrix()
            m[0, 2] *= -1.
            m[1, 2] *= -1.
            m[2, 1] *= -1.
            m[2, 0] *= -1.
            result = np.allclose(m,
                                 part3.getTransform().getMatrix())
            self.assertTrue(result)
            result = np.allclose(part1.getTransform().getMatrix(),
                                 part3.getTransform().getMatrix())
            self.assertFalse(result)
        self.assertEqual(len(prot3.outputParticles),
                         len(prot1.outputParticles))
