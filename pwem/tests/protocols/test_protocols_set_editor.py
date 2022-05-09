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


import unittest
import numpy as np

import pyworkflow.tests as pwtests
import pwem.objects as emobj
from pwem.protocols import ProtImportParticles, ProtSetEditor
from pwem.objects.data import Particle, SetOfParticles, Acquisition, CTFModel
from pwem.constants import ALIGN_PROJ
from pwem.constants import (SYM_I222, SYM_I222r)

# projection matrices
mList = [
    [[1., 0., 0., 0.],  # a1
     [0., 1., 0., 0.],
     [0., 0., 1., 0.],
     [0., 0., 0., 1.]],
    [[0.04341204, -0.82959837, 0.5566704, 7.42774284],  # -50, -40,-30
     [0.90961589, 0.26325835, 0.3213938, -20.82490128],
     [-0.41317591, 0.49240388, 0.76604444, 3.33947946],
     [0., 0., 0., 1.]],
    [[0.04341204, 0.82959837, 0.5566704, -7.42774284],  # a3
     [0.90961589, -0.26325835, 0.3213938, 20.82490128],
     [0.41317591, 0.49240388, -0.76604444, 3.33947946],
     [0., 0., 0., 1.]]
]

# operate rotate vector
tOut = [
    [[ 0.000000e-00,  0.000000e+00,  1.000000e+00,  0.000000e+00],
     [ 0.000000e+00,  1.000000e+00,  0.000000e+00,  0.000000e+00],
     [-1.000000e+00,  0.000000e+00,  0.000000e-00,  0.000000e+00],
     [ 0.000000e+00,  0.000000e+00,  0.000000e+00,  1.000000e+00]],
    [[ -0.41317591,   0.49240388,   0.76604444,   3.33947946],
     [  0.90961589,   0.26325835,   0.3213938,  -20.82490128],
     [ -0.04341204,   0.82959837,  -0.5566704,   -7.42774284],
     [  0.,           0.,           0.,           1.        ]],
    [[ 0.41317591,  0.49240388, -0.76604444,  3.33947946],
     [ 0.90961589, -0.26325835,  0.3213938,  20.82490128],
     [-0.04341204, -0.82959837, -0.5566704,   7.42774284],
     [ 0.,          0.,          0.,          1.        ]]
]
# operate rotate ico symmetry
tIcoOut = [
    [[0.,  1.,  0.,  0.],
     [-1., 0.,  0.,  0.],
     [0.,  0.,  1.,  0.],
     [0., 0., 0., 1.]],
    [[0.90961589,   0.26325835,   0.3213938, -20.82490128],
     [-0.04341204,   0.82959837, -0.5566704, -7.42774284],
     [-0.41317591,   0.49240388,  0.76604444, 3.33947946],
     [0., 0., 0., 1.]],
    [[0.90961589,  -0.26325835,  0.3213938,  20.82490128],
     [-0.04341204, -0.82959837, -0.5566704,   7.42774284],
     [0.41317591,   0.49240388, -0.76604444,  3.33947946],
     [0., 0., 0., 1.]]
]

defocusList = [15000.,20000.,25000.]
defocusAngle = [0., 10., 20.]
projSize = 128
samplingRate = 1.5

class TestSets(pwtests.BaseTest):
    """Run different tests related to the editor set protocol."""
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)

    def _createSetOfParticles(self, setPartSqliteName, partFn,
                             doCtf=False):
        # create a set of particles

        self.partSet = SetOfParticles(filename=setPartSqliteName)
        self.partSet.setAlignment(ALIGN_PROJ)
        self.partSet.setAcquisition(Acquisition(voltage=300,
                                           sphericalAberration=2,
                                           amplitudeContrast=0.1,
                                           magnification=60000))
        self.partSet.setSamplingRate(samplingRate)
        if doCtf:
            self.partSet.setHasCTF(True)
        aList = [np.array(m) for m in mList]

        for i, a in enumerate(aList):
            p = Particle()
            if doCtf:
                defocusU = defocusList[i]
                defocusV = defocusList[i]
                ctf = CTFModel(defocusU=defocusU,
                               defocusV=defocusV,
                               defocusAngle=defocusAngle[i])
                ctf.standardize()
                p.setCTF(ctf)

            p.setLocation(i + 1, partFn)
            p.setTransform(emobj.Transform(a))
            self.partSet.append(p)
        self.partSet.write()

    def importData(self, baseFn, objLabel, protType, importFrom):
        prot = self.newProtocol(protType,
                     objLabel=objLabel,
                     filesPath=baseFn,
                     maskPath=baseFn,
                     sqliteFile=baseFn,
                     haveDataBeenPhaseFlipped=False,
                     magnification=10000,
                     samplingRate=samplingRate,
                     importFrom=importFrom
                     )
        self.launchProtocol(prot)
        return prot

    def testOperation(self):
        """Make a trivial operation. Increment defocusU  by one"""
        setPartSqliteName = self.proj.getTmpPath("particles_operate.sqlite")
        setPartName = self.proj.getTmpPath('particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj   = self.importData(setPartSqliteName,
                                           "import projection\n operation",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_SCIPION)

        #launch operate set protocol
        protSetEditor = self.newProtocol(ProtSetEditor,
                                         objLabel="operate")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtSetEditor.CHOICE_FORMULA)
        protSetEditor.formula.set('item._ctfModel._defocusU.set(item._ctfModel._defocusU.get() + 1. )')
        self.launchProtocol(protSetEditor)
        for item1, item2 in zip(protSetEditor.outputParticles, protImportProj.outputParticles):
            self.assertAlmostEqual(item1._ctfModel._defocusU.get() ,
                                   item2._ctfModel._defocusU.get() + 1)

    def testRotation(self):
        """Rotate projections alignments so in the reconstruction
        vecor (0,0,1) y rotated to (1,0,0)"""
        setPartSqliteName = self.proj.getTmpPath("particles_rot_vec.sqlite")
        setPartName = self.proj.getTmpPath('particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj   = self.importData(setPartSqliteName,
                                           "import projection\n rot vector",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_SCIPION)

        #launch operate set protocol
        protSetEditor = self.newProtocol(ProtSetEditor, objLabel="rotate")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtSetEditor.CHOICE_ROTATE_VECTOR)
        protSetEditor.xs.set(0.)
        protSetEditor.ys.set(0.)
        protSetEditor.zs.set(1.)
        protSetEditor.xt.set(1.)
        protSetEditor.yt.set(0.)
        protSetEditor.zt.set(0.)
        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(  tOut,
                                protSetEditor.outputParticles, ):
            self.assertTrue(np.allclose(controlPartSet ,
                                   outPartSet.getTransform().getMatrix()))

    def testRotationIco(self):
        """Rotate projections alignments between icosahedral
        symmetries"""
        setPartSqliteName = self.proj.getTmpPath("particles_rot_ico.sqlite")
        setPartName = self.proj.getTmpPath('particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj   = self.importData(setPartSqliteName,
                                           "import projection\n ico sym",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_SCIPION)

        #launch operate set protocol
        protSetEditor = self.newProtocol(ProtSetEditor, objLabel="rotate")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtSetEditor.CHOICE_ROTATE_ICOSAHEDRAL)
        protSetEditor.originSymmetryGroup.set(SYM_I222 - SYM_I222)
        protSetEditor.targetSymmetryGroup.set(SYM_I222r - SYM_I222)

        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(tIcoOut,
                                protSetEditor.outputParticles, ):
            self.assertTrue(np.allclose(controlPartSet ,
                                   outPartSet.getTransform().getMatrix()))

if __name__ == '__main__':
    unittest.main()
