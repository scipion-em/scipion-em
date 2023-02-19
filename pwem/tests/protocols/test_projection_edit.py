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
from os.path import abspath
import inspect

import pyworkflow.tests as pwtests
import pwem.objects as emobj
from pwem.protocols import (ProtImportParticles, ProtProjectionEditor,
                            ProtImportVolumes)
from pwem.objects.data import Particle, SetOfParticles, Acquisition, CTFModel
from pwem.constants import ALIGN_PROJ
from pwem.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                            SYM_DIHEDRAL_X, SYM_DIHEDRAL_Y,
                            SYM_TETRAHEDRAL_222, SYM_TETRAHEDRAL_Z3,
                            SYM_TETRAHEDRAL_Z3R)

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
    [[0.000000e-00,  0.000000e+00,  1.000000e+00,  0.000000e+00],
     [0.000000e+00,  1.000000e+00,  0.000000e+00,  0.000000e+00],
     [-1.000000e+00,  0.000000e+00,  0.000000e-00,  0.000000e+00],
     [0.000000e+00,  0.000000e+00,  0.000000e+00,  1.000000e+00]],
    [[-0.41317591,   0.49240388,   0.76604444,   3.33947946],
     [0.90961589,   0.26325835,   0.3213938,  -20.82490128],
     [-0.04341204,   0.82959837,  -0.5566704,   -7.42774284],
     [0.,           0.,           0.,           1.]],
    [[0.41317591,  0.49240388, -0.76604444,  3.33947946],
     [0.90961589, -0.26325835,  0.3213938,  20.82490128],
     [-0.04341204, -0.82959837, -0.5566704,   7.42774284],
     [0.,          0.,          0.,          1.]]
]
tMoveOut2 = [
    [[-0.76808666,  0.51809096,  0.37633049, 0.],
     [-0.16183324,  0.41156446, -0.89689726, 0.],
     [-0.61955862, -0.74979761, -0.23227283, 0.],
     [0., 0., 0., 1.]],
    [[-0.41809412, 0.85616425, -0.30361173,  0.],
     [0.64719338, 0.5152753, 0.56181145,  0.],
     [0.6374465, 0.03839456, -0.7695374,   0.],
     [0., 0., 0., 1.]],
    [[0.86562616, 0.49977971, 0.03019261, 0.],
     [0.04570132, -0.13891773, 0.98924883, 0.],
     [0.49860078, -0.85493983, -0.14309141, 0.],
     [0.,         0.,         0.,          1.]]
]
# operate rotate vector
tRotOut = [
    [[0.5,        -0.8660254,  0.,          0.],
     [0.8660254,   0.5,        0.,          0.],
     [0.,          0.,          1.,          0.],
     [0.,          0.,          0.,          1.]],
    [[-7.66044448e-01, -6.42787604e-01,  4.58118495e-09, 2.17487650e+01],
     [4.92403874e-01, -5.86824088e-01,  6.42787608e-01, -3.97983665e+00],
     [-4.13175910e-01,  4.92403880e-01,  7.66044440e-01, 3.33947946e+00],
     [0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]],
    [[-7.66044448e-01,  6.42787604e-01,  4.58118495e-09, -2.17487650e+01],
     [4.92403874e-01,  5.86824088e-01,  6.42787608e-01, 3.97983665e+00],
     [4.13175910e-01,  4.92403880e-01, -7.66044440e-01, 3.33947946e+00],
     [0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 1.00000000e+00]]
]
tRotOut2 = [
    [[0.61955862, 0.74979761, 0.23227283, 0.],
     [0.76808666, -0.51809096, -0.37633049, 0.],
     [-0.16183324, 0.41156446, -0.89689726, 0.],
     [0., 0., 0., 1.]],
    [[-0.6374465, -0.03839456, 0.7695374, 0.],
     [0.41809412, -0.85616425, 0.30361173, 0.],
     [0.64719338, 0.5152753, 0.56181145, 0.],
     [0., 0., 0., 1.]],
    [[-0.49860078, 0.85493983, 0.14309141, 0.],
     [-0.86562616, -0.49977971, -0.03019261, 0.],
     [0.04570132, -0.13891773, 0.98924883, 0.],
     [0., 0., 0., 1.]]
]

# operate dihedral symmetry Y
tDiOutY = [
    [[-0.90145576, -0.15372423,  0.40465588,  0.],
     [-0.43275576,  0.29845287, -0.85067522,  0.],
     [0.00999868, -0.94196324, -0.33556711,  0.],
     [0.,          0.,          0.,          1.]],

    [[0.94052283,  0.2364644,  -0.24392908,  0.],
     [0.30988911, -0.89139634,  0.33072844,  0.],
     [-0.13923199, -0.38664861, -0.91165635,  0.],
     [0.,          0.,          0.,          1.]],

    [[0.11616951,  0.95112703,  0.2861154,   0.],
     [-0.25960861, -0.24897451,  0.93306755,  0.],
     [0.95870121, -0.18267202,  0.21799752,  0.],
     [0.,          0.,          0.,          1.]]

]
# operate dihedral symmetry no files downloaded
tDiOut = [
    [[0., -1.,  0., 0.],
     [1.,  0.,  0., 0.],
     [0.,  0.,  1., 0.],
     [0.,  0.,  0., 1.]],

    [[-0.90961589, -0.26325835, -0.3213938,  20.82490128],
     [0.04341204, -0.82959837,  0.5566704,   7.42774284],
     [-0.41317591,  0.49240388,  0.76604444,  3.33947946],
     [0.,          0.,          0.,          1.]],

    [[-0.90961589,   0.26325835,  -0.3213938,  -20.82490128],
     [0.04341204,   0.82959837,   0.5566704,   -7.42774284],
     [0.41317591,   0.49240388,  -0.76604444,   3.33947946],
     [0.,           0.,           0.,           1.]]
]

# operate tetrahedral symmetry
tTetraOut = [
    [[+0.70710678, -0.70710678, 0., 0.],
     [+0.40824829,  0.40824829,  0.81649658,  0.],
     [-0.57735027, -0.57735027,  0.57735027,  0.],
     [+0., 0., 0., 1.]],
    [[-0.61249862, -0.7727664,  0.16636568, 19.97763624],
     [+0.05171531,  0.17083874,  0.98394087, -2.74269347],
     [-0.78877815,  0.61126608, -0.06467464,  9.6629024],
     [+0.,  0.,  0.,  1.]],
    [[-0.61249862,   0.77276640,   0.16636568, -19.97763624],
     [+0.72642874,   0.63325343,  -0.26700446,   8.19604059],
     [-0.31168371,  -0.04268705,  -0.94922657,  -5.80680367],
     [+0., 0., 0., 1.]]
]
# operate tetrahedral symmetry t222-> tz3
tTetraOut_222_z3 = [
    [[-0.03729279, -0.09362153,  0.99490917,  0.],
     [-0.99003983,  0.13871696, -0.02405692,  0.],
     [-0.13575853, -0.98589686, -0.09786219,  0.],
     [+0.,          0.,          0.,          1.]],
    [[+0.82401380, -0.31138528, -0.47332913,  0.],
     [-0.06255460, -0.88032370,  0.47023091,  0.],
     [-0.56310583, -0.35786785, -0.74487746,  0.],
     [+0.,          0.,          0.,          1.]],
    [[+0.88956655,  0.08324056,  0.44915740,  0.],
     [-0.42656693, -0.20042519,  0.88196961,  0.],
     [+0.16343810, -0.97616635, -0.14278387,  0.],
     [+0.00000000,  0.00000000,  0.00000000,  1.]]
]

tTetraOut_222_z3r = [
    [[+0.03729279,  0.09362153, -0.99490917,  0.],
     [-0.99003983,  0.13871696, -0.02405692,  0.],
     [-0.13575853, -0.98589686, -0.09786219,  0.],
     [+0.,          0.,          0.,          1.]],
    [[-0.82401380,  0.31138528,  0.47332913,  0.],
     [-0.06255460, -0.88032370,  0.47023091,  0.],
     [-0.56310583, -0.35786785, -0.74487746,  0.],
     [+0.,          0.,          0.,          1.]],
    [[-0.88956655, -0.08324056, -0.44915740,  0.],
     [-0.42656693, -0.20042519,  0.88196961,  0.],
     [+0.16343810, -0.97616635, -0.14278387,  0.],
     [+0.,          0.,          0.,          1.]]
]

tTetraOut_z3r_z3 = [
    [[0.26401850, -0.20164842, +0.94320313, 0.],
     [0.89488254, +0.41603801, -0.16154753, 0.],
     [0.35983255, -0.88670755, -0.29029338, 0.],
     [0.,         +0.,          +0.,         1.]],
    [[0.24286299, +0.18374571, -0.95249939,  0.],
     [0.97005715, -0.04862108, +0.23796033,  0.],
     [0.00258736, +0.98177060, +0.19005210,  0.],
     [0.,         +0.,         +0.,          1.]],
    [[0.80062774, -0.42224097, +0.42509739, 0.],
     [0.35666936, +0.90595285, +0.22811487, 0.],
     [0.48143763, +0.03101587, -0.87593140, 0.],
     [0.,         +0.,         +0.,         1.]]
]

# operate ico symmetry
tIco222_222r_Out = [
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

tIcoOut_222_222r = [
    [[0.39134375,  0.61792733, -0.68192073,  0.],
     [0.83169429,  0.07966828,  0.54948846,  0.],
     [0.39387139, -0.78218845, -0.48274895,  0.],
     [0.00000000,  0.00000000,  0.00000000,  1.]],
    [[-0.99242925,  0.08980610, -0.08377972,  0.],
     [-0.12047057, -0.57908928,  0.80631411,  0.],
     [+0.02389598,  0.81030271,  0.58552413,  0.],
     [+0.00000000,  0.00000000,  0.00000000,  1.]],
    [[-0.86075607,  0.50709402,  0.04421125,  0.],
     [-0.26093301, -0.51414980,  0.81704587,  0.],
     [+0.43705028,  0.69174102,  0.57487513,  0.],
     [+0.00000000,  0.00000000,  0.00000000,  1.]]
]

tIcoOut_222r_n25 = [
    [[-0.64908948,  0.01532560, -0.76055768, 0.],
     [-0.63113568, -0.56899860,  0.52717012, 0.],
     [-0.42467706,  0.82219567,  0.37900354, 0.],
     [+0.00000000,  0.00000000,  0.00000000, 1.]],
    [[+0.23103367, -0.96948462,  0.08199404, 0.],
     [-0.37158628, -0.16580791, -0.91347215, 0.],
     [+0.89919246,  0.18057497, -0.39855439, 0.],
     [+0.00000000,  0.00000000,  0.00000000, 1.]],
    [[+0.54151319,  0.27229642, -0.79537295, 0.],
     [+0.77322535, -0.53266751,  0.3440754 , 0.],
     [-0.32997882, -0.80132389, -0.49899298, 0.],
     [+0.00000000,  0.00000000,  0.00000000, 1.]],
     ]

tIcoOut_n25_n25r = [
    [[-0.17452361,  0.16585493,  0.97058418, 0.],
     [+0.32982927,  0.93861316, -0.10108409, 0.],
     [-0.92776838,  0.30248552, -0.21851397, 0.],
     [+0.00000000,  0.00000000,  0.00000000, 1.]],
    [[+0.94832549, -0.17894467, -0.26202591, 0.],
     [-0.19912627,  0.30730334, -0.93054467, 0.],
     [+0.24703745,  0.93463547,  0.25579101, 0.],
     [+0.00000000,  0.00000000,  0.00000000, 1.]],
    [[-0.12318540, -0.90116932,  0.41559500, 0.],
     [-0.83278152, -0.13388985, -0.53716707, 0.],
     [+0.53972243, -0.41227097, -0.73398388, 0.],
     [+0.00000000,  0.00000000,  0.00000000, 1.]],
    ]

defocusList = [15000., 20000., 25000.]
defocusAngle = [0., 10., 20.]
projSize = 128
samplingRate = 1.5

# enable FULLTEST to perform a detailed comparison
# very like does not make sense for automatic testing
downloadFileFromGithub = True
# recosntruct volume, this is useful if there is a human
# to check the results and relion is installed
# xmipp is not OK because the tetrahedral symmetry reconstruction
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
        if os.path.exists(self.proj.getTmpPath(f"{baseName}")):
            os.remove(self.proj.getTmpPath(f"{baseName}"))
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
                counter += 1
            except requests.exceptions.TooManyRedirects:
                raise SystemExit("Unknown url")
            except requests.exceptions.RequestException as e:
                # catastrophic error. bail.
                raise SystemExit(e)
        if counter >= 4:
            raise SystemExit("Can connect to server")

        assert req.status_code == 200

        # save retrieved file in /tmp
        with open(self.proj.getTmpPath(f"{baseName}"), "wb") as _fh:
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
        vecor (0,0,1) becames (1,0,0). Rotation along plane
        np.cross([0,0,1], [1,0,0])"""
        funcName = inspect.stack()[0][3]
        setPartSqliteName = self.proj.getTmpPath(
            f"{funcName}__particles_rot_vec.sqlite")
        setPartName = self.proj.getTmpPath(f'{funcName}__particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj = self.importData(
            setPartSqliteName,
            "import projection\n rot vector",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_SCIPION)

        # launch operate set protocol
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="move 1,0,0 -> 0,0,1")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(
            ProtProjectionEditor.CHOICE_MOVE_VECTOR)
        protSetEditor.xs.set(0.)
        protSetEditor.ys.set(0.)
        protSetEditor.zs.set(1.)
        protSetEditor.xt.set(1.)
        protSetEditor.yt.set(0.)
        protSetEditor.zt.set(0.)
        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(tMoveOut,
                                              protSetEditor.outputParticles):
            self.assertTrue(np.allclose(controlPartSet,
                                        outPartSet.getTransform().getMatrix()))

    def test_01_MoveVector(self):
        """Rotate projections alignments so the reconstruction
        is rotated by 90 degrees
        """
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'C/C7'
        symFile = 'c7'
        # Wed Feb 8 14:27:58 2023 +0100
        # set hash = 'main' for last version
        hash = 'd8301c0274429056c2b3e98f3e0bd479bcbb1f55'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            abspath(self.proj.getTmpPath(f'{symFile}.xmd')),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)
        _ = self._importVolume(abspath(self.proj.getTmpPath(f'{symFile}.mrc')),
                               f"import vol\n {symFile}")
        # reconstruct using C7 symmetry
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                numberOfMpis=4,
                symmetryGroup='C1',
                objLabel='Fourier reconstruction',
                inputParticles=protImportProj.outputParticles)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="move 1,0,0 -> 0,0,1")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(
            ProtProjectionEditor.CHOICE_MOVE_VECTOR)
        protSetEditor.xs.set(0.)
        protSetEditor.ys.set(0.)
        protSetEditor.zs.set(1.)
        protSetEditor.xt.set(1.)
        protSetEditor.yt.set(0.)
        protSetEditor.zt.set(0.)
        editProt = self.launchProtocol(protSetEditor)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='C1',
                numberOfMpis=4,
                objLabel='Fourier reconstruction',
                inputParticles=editProt.outputParticles)

            _ = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tMoveOut2,
                                              protSetEditor.outputParticles):
            # print(controlPartSet, outPartSet.getTransform().getMatrix())
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_05_rotateAroundVector(self):
        """Rotate projections alignments around
        vector (0,0,1) by 60 degrees"""
        funcName = inspect.stack()[0][3]
        setPartSqliteName = self.proj.getTmpPath(
            f"{funcName}__particles_rot_vec.sqlite")
        setPartName = self.proj.getTmpPath(f'{funcName}__particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj = self.importData(
            setPartSqliteName,
            "import projection\n rot vector",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_SCIPION)

        # launch operate set protocol
        protSetEditor = self.newProtocol(
            ProtProjectionEditor,
            objLabel="rotate 60 degrees\n around 0,0,1")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_ROTATE_VECTOR)
        protSetEditor.x.set(0.)
        protSetEditor.y.set(0.)
        protSetEditor.z.set(1.)
        protSetEditor.angle.set(60.)
        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(tRotOut,
                                              protSetEditor.outputParticles):
            # print(outPartSet.getTransform().getMatrix())
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_06_rotateAroundVector(self):
        """Rotate projections alignments around vector x by 90
        """
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'C/C7'
        symFile = 'c7'
        # Wed Feb 8 14:27:58 2023 +0100
        hash = 'd8301c0274429056c2b3e98f3e0bd479bcbb1f55'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            abspath(self.proj.getTmpPath(f'{symFile}.xmd')),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)
        _ = self._importVolume(abspath(self.proj.getTmpPath(f'{symFile}.mrc')),
                               f"import vol\n {symFile}")
        # reconstruct using C7 symmetry
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                numberOfMpis=4,
                symmetryGroup='C1',
                objLabel='Fourier reconstruction',
                inputParticles=protImportProj.outputParticles)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="rotate 90 degree around 1,0,0")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(
            ProtProjectionEditor.CHOICE_ROTATE_VECTOR)
        protSetEditor.x.set(1.)
        protSetEditor.y.set(0.)
        protSetEditor.z.set(0.)
        protSetEditor.angle.set(90.)
        editProt = self.launchProtocol(protSetEditor)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='C1',
                numberOfMpis=4,
                objLabel='Fourier reconstruction',
                inputParticles=editProt.outputParticles)

            _ = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tRotOut2,
                                              protSetEditor.outputParticles):
            # print(controlPartSet, outPartSet.getTransform().getMatrix())
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_10_Dihedral(self):
        """Rotate projections alignments between dihedral x/y
        symmetries"""
        funcName = inspect.stack()[0][3]
        setPartSqliteName = self.proj.getTmpPath(
            f"{funcName}__particles_rot_vec.sqlite")
        setPartName = self.proj.getTmpPath(f'{funcName}__particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj = self.importData(
            setPartSqliteName,
            "import projection\n d7x",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_SCIPION)

        # launch operate set protocol
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="D7x -> D7y")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_DIHEDRAL)
        protSetEditor.originSymmetryGroupD.set(SYM_DIHEDRAL_X - SYM_DIHEDRAL_X)
        protSetEditor.targetSymmetryGroupD.set(SYM_DIHEDRAL_Y - SYM_DIHEDRAL_X)

        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(tDiOut,
                                              protSetEditor.outputParticles):
            # print(controlPartSet, outPartSet.getTransform().getMatrix())
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_11_DihedralY(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'D/dy'
        symFile = 'dy'
        # Wed Feb 8 14:27:58 2023 +0100
        hash = 'd8301c0274429056c2b3e98f3e0bd479bcbb1f55'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            self.proj.getTmpPath(f'{symFile}.xmd'),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)
        _ = self._importVolume(
                abspath(self.proj.getTmpPath(f'{symFile}.mrc')),
                f"import vol\n {symFile}")
        # reconstruct using d7y (since d7 is dx
        # this should provide a wrong result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='D7',
                objLabel='Fourier reconstruction',
                inputParticles=protImportProj.outputParticles)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="D7y -> D7x")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_DIHEDRAL)
        protSetEditor.originSymmetryGroupD.set(
            SYM_DIHEDRAL_Y - SYM_DIHEDRAL_X)
        protSetEditor.targetSymmetryGroupD.set(
            SYM_DIHEDRAL_X - SYM_DIHEDRAL_X)

        editProt = self.launchProtocol(protSetEditor)

        for controlPartSet, outPartSet in zip(tDiOutY,
                                              protSetEditor.outputParticles):
            # print(controlPartSet, outPartSet.getTransform().getMatrix())
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='D7',
                objLabel='Fourier reconstruction',
                inputParticles=editProt.outputParticles)

            _ = self.launchProtocol(recProt2)

    def test_12_DihedralX(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'D/dx'
        symFile = 'dx'
        # Wed Feb 8 14:27:58 2023 +0100
        hash = 'd8301c0274429056c2b3e98f3e0bd479bcbb1f55'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            self.proj.getTmpPath(f'{symFile}.xmd'),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)
        _ = self._importVolume(
                abspath(self.proj.getTmpPath(f'{symFile}.mrc')),
                f"import vol\n {symFile}")
        # reconstruct using d7x (since d7 is dx
        #  this should provide a good result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='D7',
                objLabel='Fourier reconstruction',
                inputParticles=protImportProj.outputParticles)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="D7x -> D7y")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_DIHEDRAL)
        protSetEditor.originSymmetryGroupD.set(
            SYM_DIHEDRAL_X - SYM_DIHEDRAL_X)
        protSetEditor.targetSymmetryGroupD.set(
            SYM_DIHEDRAL_Y - SYM_DIHEDRAL_X)

        editProt = self.launchProtocol(protSetEditor)
        # edit projection direction, put back to original
        protSetEditor2 = self.newProtocol(
            ProtProjectionEditor, objLabel="D7y -> D7x")
        protSetEditor2.inputSet.set(editProt)
        protSetEditor2.inputSet.setExtended("outputParticles")
        protSetEditor2.operation.set(ProtProjectionEditor.CHOICE_DIHEDRAL)
        protSetEditor2.originSymmetryGroupD.set(
            SYM_DIHEDRAL_Y - SYM_DIHEDRAL_X)
        protSetEditor2.targetSymmetryGroupD.set(
            SYM_DIHEDRAL_X - SYM_DIHEDRAL_X)

        editProt2 = self.launchProtocol(protSetEditor2)
        for controlPartSet, outPartSet in zip(protImportProj.outputParticles,
                                              protSetEditor2.outputParticles):
            # print(controlPartSet.getTransform().getMatrix(),
            #       outPartSet.getTransform().getMatrix())
            self.assertTrue(
                np.allclose(controlPartSet.getTransform().getMatrix(),
                            outPartSet.getTransform().getMatrix()))

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='D7',
                objLabel='Fourier reconstruction',
                inputParticles=editProt2.outputParticles)

            _ = self.launchProtocol(recProt2)

    def test_14_Tetrahedral(self):
        """Rotate projections alignments between dihedral x/y
        symmetries"""
        funcName = inspect.stack()[0][3]
        setPartSqliteName = self.proj.getTmpPath(
            f"{funcName}__particles_rot_vec.sqlite")
        setPartName = self.proj.getTmpPath(f'{funcName}__particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj = self.importData(
            setPartSqliteName,
            "import projection\n t222",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_SCIPION)

        # launch operate set protocol
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="T222 -> TZ3")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor.originSymmetryGroupT.set(
            SYM_TETRAHEDRAL_222 - SYM_TETRAHEDRAL_222)
        protSetEditor.targetSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)

        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(tTetraOut,
                                              protSetEditor.outputParticles):
            #print(controlPartSet)
            #print(repr(outPartSet.getTransform().getMatrix()))
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_15_Tetrahedral222_z3(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'T/t222'
        symFile = 't222'
        # Wed Feb 8 14:27:58 2023 +0100
        hash = 'd8301c0274429056c2b3e98f3e0bd479bcbb1f55'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            self.proj.getTmpPath(f'{symFile}.xmd'),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)

        _ = self._importVolume(
            self.proj.getTmpPath(f'{symFile}.mrc'),
            f"import vol\n {symFile}")
        # reconstruct using d7y (since d7 is dx
        #  this should provide a wrong result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='T',
                objLabel='Fourier reconstruction',
                inputParticles=protImportProj.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="T222 -> TZ3")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor.originSymmetryGroupT.set(
            SYM_TETRAHEDRAL_222 - SYM_TETRAHEDRAL_222)
        protSetEditor.targetSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)

        editProt = self.launchProtocol(protSetEditor)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='T',
                objLabel='Fourier reconstruction',
                inputParticles=editProt.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tTetraOut_222_z3,
                                              protSetEditor.outputParticles):
            # print(controlPartSet)
            # print(repr(outPartSet.getTransform().getMatrix()))
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_16_Tetrahedral222_z3r(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'T/t222'
        symFile = 't222'
        # Wed Feb 8 14:27:58 2023 +0100
        hash = 'd8301c0274429056c2b3e98f3e0bd479bcbb1f55'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            self.proj.getTmpPath(f'{symFile}.xmd'),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)
        _ = self._importVolume(
                self.proj.getTmpPath(f'{symFile}.mrc'),
                f"import vol\n {symFile}")
        # reconstruct using T (since T is tz3 this
        # should provide a wrong result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='T',
                objLabel='Fourier reconstruction',
                inputParticles=protImportProj.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="T222 -> TZ3R")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(
            ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor.originSymmetryGroupT.set(
            SYM_TETRAHEDRAL_222 - SYM_TETRAHEDRAL_222)
        protSetEditor.targetSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3R - SYM_TETRAHEDRAL_222)

        _ = self.launchProtocol(protSetEditor)
        # now let us move back to something that we know how to reconstruct
        protSetEditor2 = self.newProtocol(
            ProtProjectionEditor, objLabel="TZ3R -> TZ3")
        protSetEditor2.inputSet.set(protSetEditor)
        protSetEditor2.inputSet.setExtended("outputParticles")
        protSetEditor2.operation.set(ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor2.originSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3R - SYM_TETRAHEDRAL_222)
        protSetEditor2.targetSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)

        editProt2 = self.launchProtocol(protSetEditor2)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='T',
                objLabel='Fourier reconstruction',
                inputParticles=editProt2.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tTetraOut_222_z3r,  #
                                              protSetEditor.outputParticles):
            # print(repr(controlPartSet))
            # print(repr(outPartSet.getTransform().getMatrix()))
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_17_Tetrahedralz3r_z3(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'T/tz3r'
        symFile = 'tz3r'
        hash = 'd8301c0274429056c2b3e98f3e0bd479bcbb1f55'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            self.proj.getTmpPath(f'{symFile}.xmd'),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)
        _ = self._importVolume(
                self.proj.getTmpPath(f'{symFile}.mrc'),
                f"import vol\n {symFile}")
        # reconstruct using T (since only z3 is
        # implemented this result should be wrong)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='T',
                objLabel='Fourier reconstruction',
                inputParticles=protImportProj.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="Tz3r -> TZ3")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor.originSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3R - SYM_TETRAHEDRAL_222)
        protSetEditor.targetSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)

        editProt = self.launchProtocol(protSetEditor)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='T',
                objLabel='Fourier reconstruction',
                inputParticles=editProt.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tTetraOut_z3r_z3,  #
                                              protSetEditor.outputParticles):
            # print(controlPartSet, outPartSet.getTransform().getMatrix())
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_18_Tetrahedralz3_z3r(self):

        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'T/tz3'
        symFile = 'tz3'
        hash = 'd8301c0274429056c2b3e98f3e0bd479bcbb1f55'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj =\
            self.importData(
                self.proj.getTmpPath(f'{symFile}.xmd'),
                f"import projection\n {symFile}",
                ProtImportParticles,
                ProtImportParticles.IMPORT_FROM_XMIPP3)
        _ = self._importVolume(
                self.proj.getTmpPath(f'{symFile}.mrc'),
                f"import vol\n {symFile}")
        # reconstruct using T (z3 is implemented this result should be OK)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='T',
                objLabel='Fourier reconstruction',
                inputParticles=protImportProj.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="Tz3 -> TZ3r")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor.originSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)
        protSetEditor.targetSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3R - SYM_TETRAHEDRAL_222)

        _ = self.launchProtocol(protSetEditor)
        # undo edit
        protSetEditor2 = self.newProtocol(
            ProtProjectionEditor, objLabel="Tz3r -> TZ3")
        protSetEditor2.inputSet.set(protSetEditor)
        protSetEditor2.inputSet.setExtended("outputParticles")
        protSetEditor2.operation.set(
            ProtProjectionEditor.CHOICE_TETRAHEDRAL)
        protSetEditor2.originSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3R - SYM_TETRAHEDRAL_222)
        protSetEditor2.targetSymmetryGroupT.set(
            SYM_TETRAHEDRAL_Z3 - SYM_TETRAHEDRAL_222)

        editProt2 = self.launchProtocol(protSetEditor2)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='T',
                objLabel='Fourier reconstruction',
                inputParticles=editProt2.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(protImportProj.outputParticles,
                                              protSetEditor2.outputParticles):
            # print(controlPartSet.getTransform().getMatrix(),
            #       outPartSet.getTransform().getMatrix())
            self.assertTrue(
                np.allclose(controlPartSet.getTransform().getMatrix(),
                            outPartSet.getTransform().getMatrix()))

    def test_20_Ico222_222r(self):
        """Rotate projections alignments between icosahedral
        symmetries"""
        setPartSqliteName = self.proj.getTmpPath("particles_rot_ico.sqlite")
        setPartName = self.proj.getTmpPath('particles.stk')

        self._createSetOfParticles(setPartSqliteName, setPartName,
                                   doCtf=True)
        protImportProj = self.importData(
            setPartSqliteName,
            "import projection\n i222 sym",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_SCIPION)

        # launch operate set protocol
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="rotate")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(
            ProtProjectionEditor.CHOICE_ICOSAHEDRAL)
        protSetEditor.originSymmetryGroupI.set(SYM_I222 - SYM_I222)
        protSetEditor.targetSymmetryGroupI.set(SYM_I222r - SYM_I222)

        self.launchProtocol(protSetEditor)
        for controlPartSet, outPartSet in zip(tIco222_222r_Out,
                                              protSetEditor.outputParticles):
            self.assertTrue(
                np.allclose(
                    controlPartSet, outPartSet.getTransform().getMatrix()))

    def test_21_Icosahedral222_222r(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'I/i222'
        symFile = 'i222'
        hash = 'd8301c0274429056c2b3e98f3e0bd479bcbb1f55'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            self.proj.getTmpPath(f'{symFile}.xmd'),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)

        _ = self._importVolume(
            self.proj.getTmpPath(f'{symFile}.mrc'),
            f"import vol\n {symFile}")
        # reconstruct using d7y (since d7 is dx
        #  this should provide a wrong result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='I1',
                objLabel='Fourier reconstruction i222',
                inputParticles=protImportProj.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="I222 -> I222r")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_ICOSAHEDRAL)
        protSetEditor.originSymmetryGroupI.set(
            SYM_I222 - SYM_I222)
        protSetEditor.targetSymmetryGroupI.set(
            SYM_I222r - SYM_I222)

        editProt = self.launchProtocol(protSetEditor)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='I2',
                objLabel='Fourier reconstruction i222r',
                inputParticles=editProt.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tIcoOut_222_222r,
                                              protSetEditor.outputParticles):
            # print(controlPartSet)
            # print(repr(outPartSet.getTransform().getMatrix()))
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_22_Icosahedral222r_In25(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'I/i222r'
        symFile = 'i222r'
        hash = '190b0ad726b8afd17762fa145e1504728c4d7e8c'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'I2PC/testDataSym/{hash}/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            self.proj.getTmpPath(f'{symFile}.xmd'),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)

        _ = self._importVolume(
            self.proj.getTmpPath(f'{symFile}.mrc'),
            f"import vol\n {symFile}")
        # reconstruct using i3
        #  this should provide a good result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='I2',
                objLabel='Fourier reconstruction i222r',
                inputParticles=protImportProj.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="I222r -> In25")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_ICOSAHEDRAL)
        protSetEditor.originSymmetryGroupI.set(
            SYM_I222r - SYM_I222)
        protSetEditor.targetSymmetryGroupI.set(
            SYM_In25 - SYM_I222)

        editProt = self.launchProtocol(protSetEditor)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='I3',
                objLabel='Fourier reconstruction In25',
                inputParticles=editProt.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tIcoOut_222r_n25,
                                              protSetEditor.outputParticles):
            # print(controlPartSet)
            # print(repr(outPartSet.getTransform().getMatrix()))
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))

    def test_23_IcosahedralIn25_In25r(self):
        if not downloadFileFromGithub:
            self.assertTrue(True)
            return
        symDir = 'I/in25'
        symFile = 'in25'
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'/I2PC/testDataSym/main/{symDir}',
                                     baseName=f'{symFile}.xmd')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'/I2PC/testDataSym/main/{symDir}',
                                     baseName=f'{symFile}.mrcs')
        self._downloadFileFromGithub(host='https://raw.githubusercontent.com',
                                     dir=f'/I2PC/testDataSym/main/{symDir}',
                                     baseName=f'{symFile}.mrc')

        # import files
        protImportProj = self.importData(
            self.proj.getTmpPath(f'{symFile}.xmd'),
            f"import projection\n {symFile}",
            ProtImportParticles,
            ProtImportParticles.IMPORT_FROM_XMIPP3)

        _ = self._importVolume(
            self.proj.getTmpPath(f'{symFile}.mrc'),
            f"import vol\n {symFile}")
        # reconstruct using i3
        #  this should provide a good result)
        if reconstructVolume:
            from relion.protocols import ProtRelionReconstruct
            recProt1 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='I3',
                objLabel='Fourier reconstruction in25',
                inputParticles=protImportProj.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt1)

        # edit projection direction
        protSetEditor = self.newProtocol(
            ProtProjectionEditor, objLabel="In25 -> In25r")
        protSetEditor.inputSet.set(protImportProj)
        protSetEditor.inputSet.setExtended("outputParticles")
        protSetEditor.operation.set(ProtProjectionEditor.CHOICE_ICOSAHEDRAL)
        protSetEditor.originSymmetryGroupI.set(
            SYM_In25 - SYM_I222)
        protSetEditor.targetSymmetryGroupI.set(
            SYM_In25r - SYM_I222)

        editProt = self.launchProtocol(protSetEditor)

        # reconstruct again
        if reconstructVolume:
            recProt2 = self.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup='I4',
                objLabel='Fourier reconstruction In25r',
                inputParticles=editProt.outputParticles,
                doCTF=False)

            _ = self.launchProtocol(recProt2)

        for controlPartSet, outPartSet in zip(tIcoOut_n25_n25r,
                                              protSetEditor.outputParticles):
            #print(controlPartSet)
            #print(repr(outPartSet.getTransform().getMatrix()))
            self.assertTrue(
                np.allclose(controlPartSet,
                            outPartSet.getTransform().getMatrix()))


if __name__ == '__main__':
    unittest.main()
