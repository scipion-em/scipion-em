# ***************************************************************************
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *
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
# ***************************************************************************/
import os
from tempfile import mkstemp

import pyworkflow.tests as pwtests
import pwem.protocols as emprot

from xmipp3.protocols import XmippProtCreateMask3D


class TestOriginSampling(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        def generate(suffix="feat"):
            (fd, filename) = mkstemp(suffix=suffix)
            f = os.fdopen(fd, "w")
            fileContent = """# XMIPP_STAR_1 *
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
blo = 1 '0 0 0' '50 10.4 2'
"""
            f.write(fileContent)
            f.close()
            return filename
        cls.inFileName = generate()
        pwtests.setupTestProject(cls)

    def test_assignOriginSampling(self):
        """ Create 3D mask
        """
        args = {'source': 2, # mask from feature
                'featureFilePath': self.inFileName,
                'samplingRate': 1.1
                }

        prot = self.newProtocol(XmippProtCreateMask3D, **args)
        prot.setObjLabel('create mask 3d (feat file)')
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputMask,
                             "There was a problem when importing a 3d mask.")

        # execute protocol that modify sampling/origin
        _sampling = 2.1
        _x = 11 ; _y = 22; _z = 33
        args = {'inVolume': prot.outputMask,
                'copyFiles': False,
                'setSampling': True,
                'samplingRate': _sampling,
                'setOrigCoord': True,
                'x': _x,
                'y': _y,
                'z': _z,
                }
        prot2 = self.newProtocol(emprot.ProtOrigSampling, **args)
        prot2.setObjLabel('set sampling and orig')
        self.launchProtocol(prot2)

        # read output and check new sampling
        vol = prot2.outputVolume
        sampling = vol.getSamplingRate()
        x, y, z = vol.getOrigin(force=True).getShifts()
        self.assertAlmostEqual(sampling, _sampling)
        self.assertAlmostEqual(x, _x)
        self.assertAlmostEqual(y, _y)
        self.assertAlmostEqual(z, _z)
        print("sampling=%f"% sampling)
        print("orig=%f %f %f"% (x, y, z))
        self.assertTrue(os.path.islink(vol.getFileName()), "%s is not a link"%vol.getFileName())

    def test_assignOrigin(self):
        """ Create 3D mask
        """
        args = {'source': 2, # mask from feature
                'featureFilePath': self.inFileName,
                'samplingRate': 1.1
                }

        prot = self.newProtocol(XmippProtCreateMask3D, **args)
        prot.setObjLabel('create mask 3d (feat file)')
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputMask,
                             "There was a problem when importing a 3d mask.")

        # execute protocol that modify sampling/origin
        _sampling = 2.1
        _x = 11 ; _y = 22; _z = 33
        args = {'inVolume': prot.outputMask,
                'copyFiles': False,
                'setSampling': False,
                'setOrigCoord': True,
                'x': _x,
                'y': _y,
                'z': _z,
                }
        prot2 = self.newProtocol(emprot.ProtOrigSampling, **args)
        prot2.setObjLabel('set orig')
        self.launchProtocol(prot2)

        # read output and check new sampling
        vol = prot2.outputVolume
        sampling = vol.getSamplingRate()
        x, y, z = vol.getOrigin(force=True).getShifts()
        self.assertAlmostEqual(x, _x)
        self.assertAlmostEqual(y, _y)
        self.assertAlmostEqual(z, _z)
        print("sampling=%f"% sampling)
        print("orig=%f %f %f"% (x, y, z))
        self.assertTrue(os.path.islink(vol.getFileName()))

    def test_assignSampling(self):
        """ Create 3D mask
        """
        args = {'source': 2, # mask from feature
                'featureFilePath': self.inFileName,
                'samplingRate': 1.1
                }

        prot = self.newProtocol(XmippProtCreateMask3D, **args)
        prot.setObjLabel('create mask 3d (feat file)')
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputMask,
                             "There was a problem when importing a 3d mask.")

        # execute protocol that modify sampling/origin
        _sampling = 2.1
        args = {'inVolume': prot.outputMask,
                'copyFiles': False,
                'setSampling': True,
                'samplingRate': _sampling,
                'setOrigCoord': False,
                }
        prot2 = self.newProtocol(emprot.ProtOrigSampling, **args)
        prot2.setObjLabel('set sampling')
        self.launchProtocol(prot2)

        # read output and check new sampling
        vol = prot2.outputVolume
        sampling = vol.getSamplingRate()
        x, y, z = vol.getOrigin(force=True).getShifts()
        self.assertAlmostEqual(sampling, _sampling)
        print("sampling=%f"% sampling)
        print("orig=%f %f %f"% (x, y, z))
        self.assertTrue(os.path.islink(vol.getFileName()))

    def test_assignOriginSamplingCopyFile(self):
        """ Create 3D mask
        """
        args = {'source': 2, # mask from feature
                'featureFilePath': self.inFileName,
                'samplingRate': 1.1
                }

        prot = self.newProtocol(XmippProtCreateMask3D, **args)
        prot.setObjLabel('create mask 3d (feat file)')
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputMask,
                             "There was a problem when importing a 3d mask.")

        # execute protocol that modify sampling/origin
        _sampling = 2.1
        _x = 11 ; _y = 22; _z = 33
        args = {'inVolume': prot.outputMask,
                'copyFiles': True,
                'setSampling': True,
                'samplingRate': _sampling,
                'setOrigCoord': True,
                'x': _x,
                'y': _y,
                'z': _z,
                }
        prot2 = self.newProtocol(emprot.ProtOrigSampling, **args)
        prot2.setObjLabel('set sampling, orig\n and copy file')
        self.launchProtocol(prot2)

        # read output and check new sampling
        vol = prot2.outputVolume
        sampling = vol.getSamplingRate()
        x, y, z = vol.getOrigin(force=True).getShifts()
        self.assertAlmostEqual(sampling, _sampling)
        self.assertAlmostEqual(x, _x)
        self.assertAlmostEqual(y, _y)
        self.assertAlmostEqual(z, _z)
        print("sampling=%f"% sampling)
        print("orig=%f %f %f"% (x, y, z))
        self.assertFalse(os.path.islink(vol.getFileName()))
