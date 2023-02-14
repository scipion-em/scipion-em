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
from pwem.protocols import ProtImportParticles, ProtProjectionEditor, ProtImportVolumes
from pwem.objects.data import Particle, SetOfParticles, Acquisition, CTFModel
from pwem.constants import ALIGN_PROJ
from pwem.constants import (SYM_I222, SYM_I222r, 
    SYM_DIHEDRAL_X, SYM_DIHEDRAL_Y, 
    SYM_TETRAHEDRAL_222, SYM_TETRAHEDRAL_Z3, SYM_TETRAHEDRAL_Z3R)

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

# operate move vector
tMoveOut = [
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
# operate rotate vector
tRotOut = [
    [[ 0.5,        -0.8660254,  0.,          0.        ],
     [ 0.8660254,   0.5,        0.,          0.        ],
     [ 0.,          0.,          1.,          0.        ],
     [ 0.,          0.,          0.,          1.        ]],
    [[-7.66044448e-01, -6.42787604e-01,  4.58118495e-09,  2.17487650e+01],
     [ 4.92403874e-01, -5.86824088e-01,  6.42787608e-01, -3.97983665e+00],
     [-4.13175910e-01,  4.92403880e-01,  7.66044440e-01,  3.33947946e+00],
     [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  1.00000000e+00]],
    [[-7.66044448e-01,  6.42787604e-01,  4.58118495e-09, -2.17487650e+01],
     [ 4.92403874e-01,  5.86824088e-01,  6.42787608e-01,  3.97983665e+00],
     [ 4.13175910e-01,  4.92403880e-01, -7.66044440e-01,  3.33947946e+00],
     [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  1.00000000e+00]]
]

# operate dihedral symmetry Y
tDiOutY = [
    [[-0.90145576, -0.15372423,  0.40465588,  0.        ],
     [-0.43275576,  0.29845287, -0.85067522,  0.        ],
     [ 0.00999868, -0.94196324, -0.33556711,  0.        ],
     [ 0.,          0.,          0.,          1.        ]],

    [[ 0.94052283,  0.2364644,  -0.24392908,  0.        ],
     [ 0.30988911, -0.89139634,  0.33072844,  0.        ],
     [-0.13923199, -0.38664861, -0.91165635,  0.        ],
     [ 0.,          0.,          0.,          1.        ]],

   [[ 0.11616951,  0.95112703,  0.2861154,   0.        ],
    [-0.25960861, -0.24897451,  0.93306755,  0.        ],
    [ 0.95870121, -0.18267202,  0.21799752,  0.        ],
    [ 0.,          0.,          0.,          1.        ]]

]
# operate dihedral symmetry no files downloaded
tDiOut = [
    [[ 0., -1.,  0.,      0.],
     [ 1.,  0.,  0.,  0.],
     [ 0.,  0.,  1.,  0.],
     [ 0.,  0.,  0.,  1.]],

    [[-0.90961589, -0.26325835, -0.3213938,  20.82490128],
     [ 0.04341204, -0.82959837,  0.5566704,   7.42774284],
     [-0.41317591,  0.49240388,  0.76604444,  3.33947946],
     [ 0.,          0.,          0.,          1.        ]],

   [[ -0.90961589,   0.26325835,  -0.3213938,  -20.82490128],
    [  0.04341204,   0.82959837,   0.5566704,   -7.42774284],
    [  0.41317591,   0.49240388,  -0.76604444,   3.33947946],
    [  0.,           0.,           0.,           1.        ]]
]

# operate ico symmetry
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

# enable FULLTEST to perform a detailed comparison
# very like does not make sense for automatic testing
downloadFileFromGithub = True
# recosntruct volume, this is useful if there is a human
# to check the results and relion is isntalled
# xmipp is not OK because the tetrahedral symmetry recosntruction
# is not implemented properlly
reconstructVolume = True

class TestProjectionEdit(pwtests.BaseTest):
    """Run different tests related to the editor set protocol."""
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)

    def _downloadFileFromGithub(self, host, dir, baseName):
        """ Download a file from github."""
        import os
        if os.path.exists(f"/tmp/{baseName}"):
            os.remove(f"/tmp/{baseName}")
        import shutil
        import requests
        url = f"{host}/{dir}/{baseName}"
        counter = 0
        print("retrieving", url)
        while counter < 5:
            try:
                req = requests.get(url, stream=True)
                break
            except requests.exceptions.Timeout:
                # continue in a retry loop
                countinue += 1
            except requests.exceptions.TooManyRedirects:
                raise SystemExit("Unknown url")
            except requests.exceptions.RequestException as e:
                # catastrophic error. bail.
                raise SystemExit(e)
        if counter >= 4:
            raise SystemExit("Can connect to server")

        assert req.status_code == 200

        # save retrieved file in /tmp
        with open(f"/tmp/{baseName}", "wb") as _fh:
            req.raw.decode_content = True
            shutil.copyfileobj(req.raw, _fh)

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
        """ import proejction data"""
        prot = self.newProtocol(protType,
                     objLabel=objLabel,
                     filesPath=baseFn,
                     maskPath=baseFn,
                     sqliteFile=baseFn,
                     mdFile=baseFn,
                     haveDataBeenPhaseFlipped=False,
                     magnification=10000,
                     samplingRate=samplingRate,
                     importFrom=importFrom
                     )
        self.launchProtocol(prot)
        return prot

    def _importVolume(self, baseFn, objLabel):
        args = {'filesPath': baseFn,
                'objLabel': objLabel,
                'samplingRate': 1,
                'setOrigCoord': False
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        self.launchProtocol(protImportVol)
        volume = protImportVol.outputVolume
        return volume

    def test_00_MoveVector(self):
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
        protSetEditor = self.newProtocol(ProtProjectionEditor, objLabel="move 1,0,0 -> 0,0,1")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_MOVE_VECTOR)
        protSetEditor.xs.set(0.)
        protSetEditor.ys.set(0.)
        protSetEditor.zs.set(1.)
        protSetEditor.xt.set(1.)
        protSetEditor.yt.set(0.)
        protSetEditor.zt.set(0.)
        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(  tMoveOut,
                                protSetEditor.outputParticles, ):
            self.assertTrue(np.allclose(controlPartSet ,
                                   outPartSet.getTransform().getMatrix()))

    def test_05_rotateAroundVector(self):   ####### TODO: fix this test
        """Rotate projections alignments around
        vecor (0,0,1)by 60 degrees"""
        setPartSqliteName = self.proj.getTmpPath("particles_rot_vec.sqlite")
        setPartName = self.proj.getTmpPath('particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj   = self.importData(setPartSqliteName,
                                           "import projection\n rot vector",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_SCIPION)

        #launch operate set protocol
        protSetEditor = self.newProtocol(ProtProjectionEditor, objLabel="rotate 60 degrees\n around 0,0,1")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_ROTATE_VECTOR)
        protSetEditor.x.set(0.)
        protSetEditor.y.set(0.)
        protSetEditor.z.set(1.)
        protSetEditor.angle.set(60.)
        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(  tRotOut,
                                protSetEditor.outputParticles, ):
            print(outPartSet.getTransform().getMatrix())
            self.assertTrue(np.allclose(controlPartSet ,
                                   outPartSet.getTransform().getMatrix()))

    def test_10_Dihedral(self):
        """Rotate projections alignments between dihedral x/y
        symmetries"""
        setPartSqliteName = self.proj.getTmpPath("particles_rot_ico.sqlite")
        setPartName = self.proj.getTmpPath('particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj   = self.importData(setPartSqliteName,
                                           "import projection\n d7x",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_SCIPION)

        #launch operate set protocol
        protSetEditor = self.newProtocol(ProtProjectionEditor, objLabel="D7x -> D7y")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_DIHEDRAL)
        protSetEditor.originSymmetryGroupI.set(SYM_DIHEDRAL_X - SYM_DIHEDRAL_X)
        protSetEditor.targetSymmetryGroupI.set(SYM_DIHEDRAL_Y - SYM_DIHEDRAL_X)

        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(tDiOut,  #### change to tDiOut
                                protSetEditor.outputParticles, ):
            print(controlPartSet, outPartSet.getTransform().getMatrix())    
            self.assertTrue(np.allclose(controlPartSet ,
                                   outPartSet.getTransform().getMatrix()))

    def test_11_DihedralY(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'D/dy'
        symFile = 'dy'
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrc')

        # import files
        protImportProj   = self.importData(f'/tmp/{symFile}.xmd',
                                           f"import projection\n {symFile}",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_XMIPP3)
        protImportVol    = self._importVolume(f'/tmp/{symFile}.mrc',
                                              f"import vol\n {symFile}")
        # reconstruct using d7y (since d7 is dx this should provide the wrong result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='D7',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=protImportProj.outputParticles)

            reconstruction1 = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(ProtProjectionEditor, objLabel="D7y -> D7x")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_DIHEDRAL)
        protSetEditor.originSymmetryGroupD.set(SYM_DIHEDRAL_Y - SYM_DIHEDRAL_X)
        protSetEditor.targetSymmetryGroupD.set(SYM_DIHEDRAL_X - SYM_DIHEDRAL_X)

        editProt = self.launchProtocol(protSetEditor)
        
        for controlPartSet, outPartSet in zip(tDiOutY,
                                protSetEditor.outputParticles, ):
            print(controlPartSet, outPartSet.getTransform().getMatrix())
            self.assertTrue(np.allclose(controlPartSet ,
                                   outPartSet.getTransform().getMatrix()))

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='D7',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=editProt.outputParticles)

            reconstruction2 = self.launchProtocol(recProt2)

    def test_12_DihedralX(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'D/dx'
        symFile = 'dx'
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrc')

        # import files
        protImportProj   = self.importData(f'/tmp/{symFile}.xmd',
                                           f"import projection\n {symFile}",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_XMIPP3)
        protImportVol    = self._importVolume(f'/tmp/{symFile}.mrc',
                                              f"import vol\n {symFile}")
        # reconstruct using d7x (since d7 is dx this should provide the right result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='D7',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=protImportProj.outputParticles)

            reconstruction1 = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(ProtProjectionEditor, objLabel="D7x -> D7y")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_DIHEDRAL)
        protSetEditor.originSymmetryGroupD.set(SYM_DIHEDRAL_X - SYM_DIHEDRAL_X)
        protSetEditor.targetSymmetryGroupD.set(SYM_DIHEDRAL_Y - SYM_DIHEDRAL_X)

        editProt = self.launchProtocol(protSetEditor)
        # edit projection direction, put back to original
        protSetEditor2 = self.newProtocol(ProtProjectionEditor, objLabel="D7y -> D7x")
        protSetEditor2.inputSet.set(editProt)
        protSetEditor2.inputSet.setExtended("outputParticles")
        protSetEditor2.operation.set(ProtProjectionEditor.CHOICE_DIHEDRAL)
        protSetEditor2.originSymmetryGroupD.set(SYM_DIHEDRAL_Y - SYM_DIHEDRAL_X)
        protSetEditor2.targetSymmetryGroupD.set(SYM_DIHEDRAL_X - SYM_DIHEDRAL_X)

        editProt2 = self.launchProtocol(protSetEditor2)
        for controlPartSet, outPartSet in zip(protImportProj.outputParticles,
                                              protSetEditor2.outputParticles, ):
            print(controlPartSet.getTransform().getMatrix(), 
                  outPartSet.getTransform().getMatrix())
            self.assertTrue(np.allclose(controlPartSet.getTransform().getMatrix() ,
                                   outPartSet.getTransform().getMatrix()))
        
        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='D7',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=editProt2.outputParticles)

            reconstruction2 = self.launchProtocol(recProt2)


    def test_15_Tetrahedral222_z3(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'T/t222'
        symFile = 't222'
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrc')

        # import files
        protImportProj   = self.importData(f'/tmp/{symFile}.xmd',
                                           f"import projection\n {symFile}",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_XMIPP3)
        protImportVol    = self._importVolume(f'/tmp/{symFile}.mrc',
                                              f"import vol\n {symFile}")
        # reconstruct using d7y (since d7 is dx this should provide the wrong result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='T',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=protImportProj.outputParticles,
                                            doCTF=False)

            reconstruction1 = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(ProtProjectionEditor, objLabel="T222 -> TZ3")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor.originSymmetryGroupT.set(SYM_TETRAHEDRAL_222 - SYM_TETRAHEDRAL_222)
        protSetEditor.targetSymmetryGroupT.set(SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)

        editProt = self.launchProtocol(protSetEditor)
        
        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='T',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=editProt.outputParticles,
                                            doCTF=False)

            reconstruction2 = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tDiOutY,
                                protSetEditor.outputParticles, ):
            print(controlPartSet, outPartSet.getTransform().getMatrix())
            self.assertTrue(np.allclose(controlPartSet ,
                                   outPartSet.getTransform().getMatrix()))

    def test_16_Tetrahedralz3r_z3(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'T/tz3r'
        symFile = 'tz3r'
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrc')

        # import files
        protImportProj   = self.importData(f'/tmp/{symFile}.xmd',
                                           f"import projection\n {symFile}",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_XMIPP3)
        protImportVol    = self._importVolume(f'/tmp/{symFile}.mrc',
                                              f"import vol\n {symFile}")
        # reconstruct using T (since only z3 is implemented this result should be wrong)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='T',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=protImportProj.outputParticles,
                                            doCTF=False)

            reconstruction1 = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(ProtProjectionEditor, objLabel="Tz3r -> TZ3")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor.originSymmetryGroupT.set(SYM_TETRAHEDRAL_Z3R - SYM_TETRAHEDRAL_222)
        protSetEditor.targetSymmetryGroupT.set(SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)

        editProt = self.launchProtocol(protSetEditor)
        
        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='T',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=editProt.outputParticles,
                                            doCTF=False)

            reconstruction2 = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tDiOutY,
                                protSetEditor.outputParticles, ):
            print(controlPartSet, outPartSet.getTransform().getMatrix())
            self.assertTrue(np.allclose(controlPartSet ,
                                   outPartSet.getTransform().getMatrix()))


    def test_17_Tetrahedralz3_z3r(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'T/tz3'
        symFile = 'tz3'
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host = 'https://raw.githubusercontent.com',
                                    dir = f'/I2PC/testDataSym/main/{symDir}',
                                    baseName=f'{symFile}.mrc')

        # import files
        protImportProj   = self.importData(f'/tmp/{symFile}.xmd',
                                           f"import projection\n {symFile}",
                                           ProtImportParticles,
                                           ProtImportParticles.IMPORT_FROM_XMIPP3)
        protImportVol    = self._importVolume(f'/tmp/{symFile}.mrc',
                                              f"import vol\n {symFile}")
        # reconstruct using T (since only z3 is implemented this result should be OK)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='T',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=protImportProj.outputParticles,
                                            doCTF=False)

            reconstruction1 = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(ProtProjectionEditor, objLabel="Tz3 -> TZ3r")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor.originSymmetryGroupT.set(SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)
        protSetEditor.targetSymmetryGroupT.set(SYM_TETRAHEDRAL_Z3R - SYM_TETRAHEDRAL_222)

        editProt = self.launchProtocol(protSetEditor)
        # undo edit
        protSetEditor2 = self.newProtocol(ProtProjectionEditor, objLabel="Tz3r -> TZ3")
        protSetEditor2.inputSet.set(protSetEditor)
        protSetEditor2.inputSet.setExtended("outputParticles")
        protSetEditor2.operation.set(ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor2.originSymmetryGroupT.set(SYM_TETRAHEDRAL_Z3R - SYM_TETRAHEDRAL_222)
        protSetEditor2.targetSymmetryGroupT.set(SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)

        editProt2 = self.launchProtocol(protSetEditor2)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(ProtRelionReconstruct,
                                            symmetryGroup='T',
                                            objLabel='Fourier reconstruction',
                                            inputParticles=editProt2.outputParticles,
                                            doCTF=False)

            reconstruction2 = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(protImportProj.outputParticles,
                                              protSetEditor2.outputParticles, ):
            print(controlPartSet.getTransform().getMatrix(), 
                  outPartSet.getTransform().getMatrix())
            self.assertTrue(np.allclose(controlPartSet.getTransform().getMatrix() ,
                                   outPartSet.getTransform().getMatrix()))
if __name__ == '__main__':
    unittest.main()
