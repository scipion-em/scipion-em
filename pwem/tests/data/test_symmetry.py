# **************************************************************************
# *
# * Authors:     roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.,  See the
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
#
import numpy as np
import pyworkflow.tests as pwtests

import pwem.constants as emcts
import pwem.convert as emconv
import os

try:
    from itertools import izip
except ImportError:
    izip = zip


class TestSymmetry(pwtests.unittest.TestCase):

    def assertArrayAlmostEqual(self, a1, a2, decimal=3):
        try:
            np.testing.assert_array_almost_equal(a1, a2, decimal=decimal)
            res = True
        except AssertionError as err:
            res = False
            print(err)
        self.assertTrue(res)

    def test_01_SymmetryCyclicSymmetryMatrices(self):
        n = 7
        matrices = emconv.getSymmetryMatrices(emcts.SYM_CYCLIC, n=n)
        refMatrices = []
        angle = 2 * np.pi / n
        for i in range(n):
            c = np.cos(angle * i)
            s = np.sin(angle * i)
            refMatrices.append([[c, -s, 0, 0],
                                [s, c, 0, 0],
                                [0, 0, 1.0, 0],
                                [0, 0, 0.0, 1.0]])
        for i, (m, r) in enumerate(izip(matrices, refMatrices)):
            print(f"Symmetry matrix {i}:\n ", m)
            self.assertArrayAlmostEqual(r, m)

    def test_02_SymmetryCyclicUnitCell(self):
        n = 7
        angle = 2 * np.pi / n
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_CYCLIC,
                                                       n=n, 
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0), 
                                                       offset=None)
        v = []
        v.append([np.cos(angle / 2.), np.sin(angle / 2.), 0])
        v.append([np.cos(-angle / 2.), np.sin(-angle / 2.), 0])
        for v, w in zip(v, vectorsEdge):
            self.assertArrayAlmostEqual(v, w, decimal=3)


    def test_10_SymmetryDihedralXSymmetryMatrices(self):
        n = 7
        angle = 2 * np.pi / n
        matrices = emconv.getSymmetryMatrices(emcts.SYM_DIHEDRAL_X, n=7)
        refMatrices = []
        # this are the matrices for the 7-fold dihedral symmetry
        for i in range(n):
            c = np.cos(angle * i)
            s = np.sin(angle * i)
            refMatrices.append([[c, -s, 0, 0],
                                [s, c, 0, 0],
                                [0, 0, 1.0, 0],
                                [0, 0, 0.0, 1.0]])

        # reflection x axis
        for i in range(n):
            mat = np.copy(refMatrices[i])
            
            mat[1][1] *= -1
            mat[1][0] *= -1
            mat[2][2] *= -1
            refMatrices.append(mat)

        for i, (m, r) in enumerate(izip(matrices, refMatrices)):
            print(f"Symmetry matrix {i}:\n ", m)
            self.assertArrayAlmostEqual(r, m)

    def test_11_SymmetryDihedralYSymmetryMatrices(self):
        n = 7
        angle = 2 * np.pi / n
        matrices = emconv.getSymmetryMatrices(emcts.SYM_DIHEDRAL_Y, n=7)
        refMatrices = []
        # this are the matrices for the 7-fold dihedral symmetry
        for i in range(n):
            c = np.cos(angle * i)
            s = np.sin(angle * i)
            refMatrices.append([[c, -s, 0, 0],
                                [s, c, 0, 0],
                                [0, 0, 1.0, 0],
                                [0, 0, 0.0, 1.0]])

        # reflection y axis
        for i in range(n):
            mat = np.copy(refMatrices[i])
            
            mat[0][0] *= -1
            mat[0][1] *= -1
            mat[2][2] *= -1
            refMatrices.append(mat)

        for i, (m, r) in enumerate(izip(matrices, refMatrices)):
            print(f"Symmetry matrix {i}:\n ", m)
            self.assertArrayAlmostEqual(r, m)

    def test_13_SymmetrydihedralXUnitCell(self):
        n = 7
        angle = np.pi / n  + np.pi/2.
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_DIHEDRAL_X,
                                                       n=n, 
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0), 
                                                       offset=None)
        v = []
        v.append([np.cos(angle), np.sin(angle), 0])
        v.append([-np.cos(angle), np.sin(angle), 0])
        v.append([0, 0, 1])
        for v, w in zip(v, vectorsEdge):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_14_SymmetrydihedralYUnitCell(self):
        n = 7
        angle = np.pi / n  + np.pi/2.
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_DIHEDRAL_Y,
                                                       n=n, 
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0), 
                                                       offset=None)
        v = []
        v.append([np.cos(angle), np.sin(angle), 0])
        v.append([-np.cos(angle), np.sin(angle), 0])
        v.append([0, 0, 1])
        for v, w in zip(v, vectorsEdge):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_21_SymmetryTetrahedral222SymmetryMatrices(self):
        n = 7 # dummy n since T symmetry does not require order but general function does
        matrices = emconv.getSymmetryMatrices(emcts.SYM_TETRAHEDRAL_222, n=n)
        refMatrices = [
                       [[1.0, 0.0, 0,   0],  #0
                        [0.0, 1.0, 0,   0], 
                        [0,   0,   1.,  0],
                        [0,   0,   0,   1.]],
                       [[1.0, 0.0, 0,   0], #1
                        [0.0, -1.0, 0,  0],
                        [0,   0,   -1,  0],
                        [0,   0,   0,   1]],
                        [[-1.0, 0.0, 0,  0], #2
                            [0.0, 1.0, 0,   0],
                            [0,   0,   -1,  0],
                            [0,   0,   0,   1]],
                        [[-1.0, 0.0, 0,  0], #3
                            [0.0, -1.0, 0,  0],
                            [0,   0,   1,   0],
                            [0,   0,   0,   1]],
                        [[0.0, 0.0, 1.0,   0], #4
                            [1.0, 0.0, 0,   0],
                            [0,   1.0,   0,  0],
                            [0,   0,   0,   1]],
                            [[0.0, 1.0, 0.0,  0], #5
                            [0,   0,   1.0, 0],
                            [1.0, 0.0, 0,  0],
                            [0,   0,   0,   1]],
                        [[0.0, 0.0, -1.0,  0], #6
                            [1.0, 0.0, 0,   0],
                            [0,   -1.0,   0,  0],
                            [0,   0,   0,   1]],
                        [[0.0, 1.0, 0.0,   0], #7
                            [ 0,   0.0, -1,  0],
                            [-1.0, 0.0, 0,  0],
                            [0,   0,   0,   1]],
                        [[0.0, 0.0, 1.0,  0], #8
                            [ -1,   0.0, 0,  0],
                            [0,   -1.0, 0,  0],
                            [0,   0,   0,   1]],
                        [[0.0, -1.0, 0.0,   0], #9
                            [0.0, 0.0, -1.0,  0],
                            [1.0, 0.0, 0,   0],
                            [0,   0,   0,   1]],
                        [[0.0, 0.0, -1.0,   0], #10
                            [-1.0, 0.0, 0,  0],
                            [0,   1.0, 0.0,  0],
                            [0,   0,   0,   1]],
                        [[0.0, -1.0, 0.0,   0], #11
                            [0.0, 0.0, 1.0,  0],
                            [-1.0, 0.0, 0,  0],
                            [0,   0,   0,   1]],
        ]
        for i, (m, r) in enumerate(izip(matrices, refMatrices)):
             self.assertArrayAlmostEqual(r, m)

        #t222 phantom
        v = []
        v1 = np.array([ 0.5773502691896258,  0.5773502691896258,  0.5773502691896258, 1.])
        v.append(np.array([-0.5773502691896258, -0.5773502691896258,  0.5773502691896258, 1.]))
        v.append(np.array([ 0.5773502691896258, -0.5773502691896258, -0.5773502691896258, 1.]))
        v.append(np.array([-0.5773502691896258,  0.5773502691896258, -0.5773502691896258, 1.]))

        # some matrix should transform v1 into vX
        for j, w in enumerate(v):
            loopEnd = False
            for i, m in enumerate(matrices):
                vv = np.dot(m, w)
                res = np.allclose(v1, vv, atol=.001)
                if res:
                    print("found match", m, w, vv)
                    loopEnd = True
                    break
            self.assertTrue(loopEnd, f"No matrix transform vector {w} into vector {v}" )

    def test_22_SymmetryTetrahedralZ3RSymmetryMatrices(self):
        n = 7
        matrices = emconv.getSymmetryMatrices(emcts.SYM_TETRAHEDRAL_Z3R, n=n)
        print(np.round(matrices, 7))
        refMatrices = [
            [[1.,         0.,        -0.,         0.,       ],
            [ 0.,         1.,        -0.,         0.,       ],
            [ 0.,        -0.,         1.,         0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.,         0.5773503,  0.8164966,  0.,       ],
            [ 0.5773503, -0.6666667,  0.4714045,  0.,       ],
            [ 0.8164966,  0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.,       -0.5773503, -0.8164966,  0.,       ],
            [-0.5773503, -0.6666667,  0.4714045,  0.,       ],
            [-0.8164966,  0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-1.,        0.,        -0.,         0.,       ],
            [ 0.,         0.3333333, -0.942809,   0.,       ],
            [ 0.,        -0.942809,  -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.5,      -0.8660254,  0.,         0.,       ],
            [ 0.8660254, -0.5,       -0.,         0.,       ],
            [ 0.,         0.,         1.,         0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.5,       0.8660254,  0.,         0.,       ],
            [-0.8660254, -0.5,         0.,         0.,       ],
            [ 0.,         0.,         1.,         0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.5,       0.2886751, -0.8164966,  0.,       ],
            [-0.2886751,  0.8333333,  0.4714045,  0.,       ],
            [ 0.8164966,  0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.5,      -0.2886751,  0.8164966,  0.,       ],
            [ 0.2886751,  0.8333333,  0.4714045,  0.,       ],
            [-0.8164966,  0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.5,       -0.2886751,  0.8164966,  0.,       ],
            [-0.8660254, -0.1666667,  0.4714045,  0.,       ],
            [-0.,        -0.942809,  -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.5,       -0.8660254,  0.,         0.,       ],
            [-0.2886751, -0.1666667, -0.942809,   0.,       ],
            [ 0.8164966,  0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.5,        0.8660254, -0.,         0.,       ],
            [ 0.2886751, -0.1666667, -0.942809,   0.,       ],
            [-0.8164966,  0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.5,        0.2886751, -0.8164966,  0.,       ],
            [ 0.8660254, -0.1666667,  0.4714045,  0.,       ],
            [-0.,        -0.942809,  -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]]
        ]


        for i, (m, r) in enumerate(izip(matrices, refMatrices)):
            # print(f"Symmetry matrix {i}:\n ", m, r, type(m), type(r))
            self.assertArrayAlmostEqual(m, r)

        #tz3r phantom
        v = []
        v1 = np.array([0, -0.9428090415820634, -0.333333, 1.])
        v.append(np.array([0.816496580927726,  0.4714045207910317,  -0.333333, 1.]))
        v.append(np.array([-0.816496580927726, 0.4714045207910317,  -0.333333, 1.]))
        v.append(np.array([0, 0, 1, 1.]))

        # some matrix should transform v1 into vX
        for j, w in enumerate(v):
            loopEnd = False
            for i, m in enumerate(matrices):
                vv = np.dot(m, w)
                res = np.allclose(v1, vv, atol=.001)
                if res:
                    # print("found match", m, w, vv)
                    loopEnd = True
                    break
            self.assertTrue(loopEnd, f"No matrix transform vector {w} into vector {v}" )

    def test_23_SymmetryTetrahedralZ3SymmetryMatrices(self):
        n = 7
        matrices = emconv.getSymmetryMatrices(emcts.SYM_TETRAHEDRAL_Z3, n=n)
        #print(np.round(matrices, 7))
        refMatrices = [
            [[1.,         0.,         0.,         0.,       ], #0
            [ 0.,         1.,         0.,         0.,       ],
            [ 0.,         0.,         1.,         0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.,         0.5773503, -0.8164966,  0.,       ], #1
            [ 0.5773503, -0.6666667, -0.4714045,  0.,       ],
            [-0.8164966, -0.4714045,  -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.,       -0.5773503,  0.8164966,  0.,       ], #2
            [-0.5773503, -0.6666667, -0.4714045,  0.,       ],
            [ 0.8164966, -0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-1.,        0.,         0.,         0.,       ], #3
            [ 0.,         0.3333333,  0.942809,   0.,       ],
            [-0.,         0.942809,  -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.5,       -0.8660254, -0.,         0.,       ], #4
            [ 0.8660254, -0.5,        0.,         0.,       ],
            [ 0.,        -0.,         1.,         0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.5,       0.8660254, -0.,         0.,       ], #5
            [-0.8660254, -0.5,       -0.,         0.,       ],
            [-0.,         0.,         1.,         0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.5,       0.2886751,  0.8164966,  0.,       ], #6
            [-0.2886751,  0.8333333, -0.4714045,  0.,       ],
            [-0.8164966, -0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[-0.5,       -0.2886751, -0.8164966,  0.,       ], #7
            [ 0.2886751,  0.8333333, -0.4714045,  0.,       ],
            [ 0.8164966, -0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.5,       -0.2886751, -0.8164966,  0.,       ], #8
            [-0.8660254, -0.1666667, -0.4714045,  0.,       ],
            [ 0.,         0.942809,  -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.5,       -0.8660254, -0.,         0.,       ], #9
            [-0.2886751, -0.1666667,  0.942809,   0.,       ],
            [-0.8164966, -0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.5,        0.8660254,  0.,         0.,       ], #10
            [ 0.2886751, -0.1666667,  0.942809,   0.,       ],
            [ 0.8164966, -0.4714045, -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]],

            [[0.5,        0.2886751,  0.8164966,  0.,       ], #11
            [ 0.8660254, -0.1666667, -0.4714045,  0.,       ],
            [ 0.,         0.942809,  -0.3333333,  0.,       ],
            [ 0.,         0.,         0.,         1.,       ]]
        ]

        for i, (m, r) in enumerate(izip(matrices, refMatrices)):
            # print(f"Symmetry matrix {i}:\n ", m, r, type(m), type(r))
            self.assertArrayAlmostEqual(m, r)

        #tz3 phantom
        v = []
        v1 = np.array([0, 0.9428090415820634, -0.333333, 1.])
        v.append(np.array([0.816496580927726, -0.4714045207910317,  -0.333333, 1.]))
        v.append(np.array([-0.816496580927726, -0.4714045207910317,  -0.333333, 1.]))
        v.append(np.array([0, 0, 1, 1.]))

        # some matrix should transform v1 into vX
        for j, w in enumerate(v):
            loopEnd = False
            for i, m in enumerate(matrices):
                vv = np.dot(m, w)
                res = np.allclose(v1, vv, atol=.001)
                if res:
                    # print("found match", m, w, vv)
                    loopEnd = True
                    break
            self.assertTrue(loopEnd, f"No matrix transform vector {w} into vector {v}" )


    def test_24_SymmetryTetahedral222UnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_TETRAHEDRAL_222,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [ [-0.7071067811865475,  0.7071067811865475,  0.,        ],
                           [ 0.,         -0.7071067811865475,  0.7071067811865475],
                           [0.7071067811865475,   0.,          0.7071067811865475]]
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_25_SymmetryTetahedralZ3RUnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_TETRAHEDRAL_Z3R,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        
        #print("vectorsPlane", vectorsPlane)
        refereceVectors= [ [-1.0, 0.0, 0.0],
                            [0.5, 0.8660254037844387, 0.],
                            [0.5, 0.28867513459481275, -0.816496580927726]
        ]                            
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_26_SymmetryTetahedralZ3UnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_TETRAHEDRAL_Z3,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [  [-1.0, 1.1775693440128312e-16, 0.0],
                            [0.5, 0.28867513459481275, 0.816496580927726],
                            [0.5000000000000001, 0.8660254037844387, 0.0]
        ]                            
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_31_SymmetryOctahedral(self):
        matrices = emconv.getSymmetryMatrices(emcts.SYM_OCTAHEDRAL)
        refMatrices = [
            [[1., 0., 0., 0.],  # 0
             [0., 1., 0., 0.],
             [0., 0., 1., 0.],
             [0., 0., 0., 1.]],
            [[0., -1., 0., 0.],  # 1
             [1., 0., 0., 0.],
             [0., 0., 1., 0.],
             [0., 0., 0., 1.]],
            [[-1., 0., 0., 0.],  # 2
             [0., -1., 0., 0.],
             [0., 0., 1., 0.],
             [0., 0., 0., 1.]],
            [[0., 1., 0., 0.],  # 3
             [-1., 0., 0., 0.],
             [0., 0., 1., 0.],
             [0., 0., 0., 1.]],

            [[1., 0., 0., 0.],  # 4
             [0., 0., -1., 0.],
             [0., 1., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., -1., 0., 0.],  # 5
             [0., 0., -1., 0.],
             [1., 0., 0., 0.],
             [0., 0., 0., 1.]],
            [[-1., 0., 0., 0.],  # 6
             [0., 0., -1., 0.],
             [0., -1., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., 1., 0., 0.],  # 7
             [0., 0., -1., 0.],
             [-1., 0., 0., 0.],
             [0., 0., 0., 1.]],

            [[1., 0., 0., 0.],  # 8
             [0., -1., 0., 0.],
             [0., 0., -1., 0.],
             [0., 0., 0., 1.]],
            [[0., -1., 0., 0.],  # 9
             [-1., 0., 0., 0.],
             [0., 0., -1., 0.],
             [0., 0., 0., 1.]],
            [[-1., 0., 0., 0.],  # 10
             [0., 1., 0., 0.],
             [0., 0., -1., 0.],
             [0., 0., 0., 1.]],
            [[0., 1., 0., 0.],  # 11
             [1., 0., 0., 0.],
             [0., 0., -1., 0.],
             [0., 0., 0., 1.]],

            [[1., 0., 0., 0.],  # 12
             [0., 0., 1., 0.],
             [0., -1., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., -1., 0., 0.],  # 13
             [0., 0., 1., 0.],
             [-1., 0., 0., 0.],
             [0., 0., 0., 1.]],
            [[-1., 0., 0., 0.],  # 14
             [0., 0., 1., 0.],
             [0., 1., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., 1., 0., 0.],  # 15
             [0., 0., 1., 0.],
             [1., 0., 0., 0.],
             [0., 0., 0., 1.]],

            [[0., 0., 1., 0.],  # 16
             [0., 1., 0., 0.],
             [-1., 0., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., 0., 1., 0.],  # 17
             [1., 0., 0., 0.],
             [0., 1., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., 0., 1., 0.],  # 18
             [0., -1., 0., 0.],
             [1., 0., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., 0., 1., 0.],  # 19
             [-1., 0., 0., 0.],
             [0., -1., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., 0., -1., 0.],  # 20
             [0., 1., 0., 0.],
             [1., 0., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., 0., -1., 0.],  # 21
             [1., 0., 0., 0.],
             [0., -1., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., 0., -1., 0.],  # 22
             [0., -1., 0., 0.],
             [-1., 0., 0., 0.],
             [0., 0., 0., 1.]],
            [[0., 0., -1., 0.],  # 23
             [-1., 0., 0., 0.],
             [0., 1., 0., 0.],
             [0., 0., 0., 1.]],
        ]

        for m1, m2 in zip(matrices[:len(refMatrices)], refMatrices):
            self.assertArrayAlmostEqual(m1, m2)

    def test_32_SymmetryOctahedralUnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_OCTAHEDRAL,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [ [0.7071067811865476, 0.7071067811865476, -0.0],
                           [0.0, -0.7071067811865475, 0.7071067811865475],
                           [-0.7071067811865476, 0.7071067811865476, 0.0]

        ]                            
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)


    def test_41_SymmetryIcosahedral222(self):
        matrices = emconv.getSymmetryMatrices(emcts.SYM_I222)
        refMatrices = [
            [[1., 0., 0., 0.],  # 1
             [0., 1., 0., 0.],
             [0., 0., 1., 0.],
             [0., 0., 0., 1.]],

            [[-0.80901699, -0.5, 0.30901699, 0.],  # 2
             [-0.5, 0.30901699, -0.80901699, 0.],
             [0.30901699, -0.80901699, -0.5, 0.],
             [0., 0., 0., 1.]],

            [[0., 1., 0., 0.],  # 3
             [0., 0., -1., 0.],
             [-1., 0., 0., 0.],
             [0., 0., 0., 1.]],

            [[0.80901699, -0.5, -0.30901699, 0.],  # 4
             [-0.5, -0.30901699, -0.80901699, 0.],
             [0.30901699, 0.80901699, -0.5, 0.],
             [0., 0., 0., 1.]],

            [[0.5, 0.30901699, -0.80901699, 0.],  # 5
             [-0.30901699, -0.80901699, -0.5, 0.],
             [-0.80901699, 0.5, -0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.30901699, -0.80901699, -0.5, 0.],  # 6
             [-0.80901699, 0.5, -0.30901699, 0.],
             [0.5, 0.30901699, -0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.80901699, 0.5, -0.30901699, 0.],  # 7
             [0.5, 0.30901699, -0.80901699, 0.],
             [-0.30901699, -0.80901699, -0.5, 0.],
             [0., 0., 0., 1.]],

            [[-0.80901699, -0.5, -0.30901699, 0.],  # 8
             [0.5, -0.30901699, -0.80901699, 0.],
             [0.30901699, -0.80901699, 0.5, 0.],
             [0., 0., 0., 1.]],

            [[-0.30901699, 0.80901699, -0.5, 0.],  # 9
             [-0.80901699, -0.5, -0.30901699, 0.],
             [-0.5, 0.30901699, 0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[0.5, -0.30901699, -0.80901699, 0.],  # 10
             [-0.30901699, 0.80901699, -0.5, 0.],
             [0.80901699, 0.5, 0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0., 0., -1., 0.],  # 11
             [-1., 0., 0., 0.],
             [0., 1., 0., 0.],
             [0., 0., 0., 1.]],

            [[-0.5, -0.30901699, -0.80901699, 0.],  # 12
             [0.30901699, 0.80901699, -0.5, 0.],
             [0.80901699, -0.5, -0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.5, 0.30901699, -0.80901699, 0.],  # 13
             [0.30901699, -0.80901699, -0.5, 0.],
             [-0.80901699, -0.5, 0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.30901699, -0.80901699, -0.5, 0.],  # 14
             [0.80901699, -0.5, 0.30901699, 0.],
             [-0.5, -0.30901699, 0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.80901699, 0.5, -0.30901699, 0.],  # 15
             [-0.5, -0.30901699, 0.80901699, 0.],
             [0.30901699, 0.80901699, 0.5, 0.],
             [0., 0., 0., 1.]],

            [[0.5, 0.30901699, -0.80901699, 0.],  # 16
             [0.30901699, 0.80901699, 0.5, 0.],
             [0.80901699, -0.5, 0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.5, 0.30901699, -0.80901699, 0.],  # 17
             [-0.30901699, 0.80901699, 0.5, 0.],
             [0.80901699, 0.5, -0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0., 0., -1., 0.],
             [1., 0., 0., 0.],
             [0., -1., 0., 0.],
             [0., 0., 0., 1.]],

            [[-0.5, -0.30901699, -0.80901699, 0.],
             [-0.30901699, -0.80901699, 0.5, 0.],
             [-0.80901699, 0.5, 0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0., -1., 0., 0.],
             [0., 0., 1., 0.],
             [-1., 0., 0., 0.],
             [0., 0., 0., 1.]],

            [[-0.80901699, 0.5, 0.30901699, 0.],
             [0.5, 0.30901699, 0.80901699, 0.],
             [0.30901699, 0.80901699, -0.5, 0.],
             [0., 0., 0., 1.]],

            [[0.80901699, 0.5, -0.30901699, 0.],
             [0.5, -0.30901699, 0.80901699, 0.],
             [0.30901699, -0.80901699, -0.5, 0.],
             [0., 0., 0., 1.]],

            [[-0.30901699, 0.80901699, -0.5, 0.],
             [0.80901699, 0.5, 0.30901699, 0.],
             [0.5, -0.30901699, -0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[0.5, -0.30901699, -0.80901699, 0.],
             [0.30901699, -0.80901699, 0.5, 0.],
             [-0.80901699, -0.5, -0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.80901699, -0.5, -0.30901699, 0.],
             [-0.5, 0.30901699, 0.80901699, 0.],
             [-0.30901699, 0.80901699, -0.5, 0.],
             [0., 0., 0., 1.]],

            [[-0.30901699, -0.80901699, 0.5, 0.],
             [-0.80901699, 0.5, 0.30901699, 0.],
             [-0.5, -0.30901699, -0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.30901699, 0.80901699, 0.5, 0.],
             [0.80901699, 0.5, -0.30901699, 0.],
             [-0.5, 0.30901699, -0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[1., 0., 0., 0.],
             [0., -1., 0., 0.],
             [0., 0., -1., 0.],
             [0., 0., 0., 1.]],

            [[0.30901699, 0.80901699, -0.5, 0.],
             [0.80901699, -0.5, -0.30901699, 0.],
             [-0.5, -0.30901699, -0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[0.30901699, -0.80901699, -0.5, 0.],
             [-0.80901699, -0.5, 0.30901699, 0.],
             [-0.5, 0.30901699, -0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[-1., 0., 0., 0.],
             [0., 1., 0., 0.],
             [0., 0., -1., 0.],
             [0., 0., 0., 1.]],

            [[0.80901699, 0.5, -0.30901699, 0.],
             [-0.5, 0.30901699, -0.80901699, 0.],
             [-0.30901699, 0.80901699, 0.5, 0.],
             [0., 0., 0., 1.]],

            [[0., -1., 0., 0.],
             [0., 0., -1., 0.],
             [1., 0., 0., 0.],
             [0., 0., 0., 1.]],

            [[-0.80901699, 0.5, 0.30901699, 0.],
             [-0.5, -0.30901699, -0.80901699, 0.],
             [-0.30901699, -0.80901699, 0.5, 0.],
             [0., 0., 0., 1.]],

            [[-0.5, -0.30901699, 0.80901699, 0.],
             [-0.30901699, -0.80901699, -0.5, 0.],
             [0.80901699, -0.5, 0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0.30901699, 0.80901699, 0.5, 0.],
             [-0.80901699, 0.5, -0.30901699, 0.],
             [-0.5, -0.30901699, 0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[0.80901699, -0.5, 0.30901699, 0.],
             [0.5, 0.30901699, -0.80901699, 0.],
             [0.30901699, 0.80901699, 0.5, 0.],
             [0., 0., 0., 1.]],

            [[0.80901699, 0.5, 0.30901699, 0.],
             [0.5, -0.30901699, -0.80901699, 0.],
             [-0.30901699, 0.80901699, -0.5, 0.],
             [0., 0., 0., 1.]],

            [[0.30901699, -0.80901699, 0.5, 0.],
             [-0.80901699, -0.5, -0.30901699, 0.],
             [0.5, -0.30901699, -0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.5, 0.30901699, 0.80901699, 0.],
             [-0.30901699, 0.80901699, -0.5, 0.],
             [-0.80901699, -0.5, -0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0., 0., 1., 0.],
             [-1., 0., 0., 0.],
             [0., -1., 0., 0.],
             [0., 0., 0., 1.]],

            [[0.5, 0.30901699, 0.80901699, 0.],
             [0.30901699, 0.80901699, -0.5, 0.],
             [-0.80901699, 0.5, 0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0.5, -0.30901699, 0.80901699, 0.],
             [0.30901699, -0.80901699, -0.5, 0.],
             [0.80901699, 0.5, -0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0.30901699, 0.80901699, 0.5, 0.],
             [0.80901699, -0.5, 0.30901699, 0.],
             [0.5, 0.30901699, -0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[0.80901699, -0.5, 0.30901699, 0.],
             [-0.5, -0.30901699, 0.80901699, 0.],
             [-0.30901699, -0.80901699, -0.5, 0.],
             [0., 0., 0., 1.]],

            [[-0.5, -0.30901699, 0.80901699, 0.],
             [0.30901699, 0.80901699, 0.5, 0.],
             [-0.80901699, 0.5, -0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0.5, -0.30901699, 0.80901699, 0.],
             [-0.30901699, 0.80901699, 0.5, 0.],
             [-0.80901699, -0.5, 0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0., 0., 1., 0.],
             [1., 0., 0., 0.],
             [0., 1., 0., 0.],
             [0., 0., 0., 1.]],

            [[0.5, 0.30901699, 0.80901699, 0.],
             [-0.30901699, -0.80901699, 0.5, 0.],
             [0.80901699, -0.5, -0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0., 1., 0., 0.],
             [0., 0., 1., 0.],
             [1., 0., 0., 0.],
             [0., 0., 0., 1.]],

            [[0.80901699, -0.5, -0.30901699, 0.],
             [0.5, 0.30901699, 0.80901699, 0.],
             [-0.30901699, -0.80901699, 0.5, 0.],
             [0., 0., 0., 1.]],

            [[-0.80901699, -0.5, 0.30901699, 0.],
             [0.5, -0.30901699, 0.80901699, 0.],
             [-0.30901699, 0.80901699, 0.5, 0.],
             [0., 0., 0., 1.]],

            [[0.30901699, -0.80901699, 0.5, 0.],
             [0.80901699, 0.5, 0.30901699, 0.],
             [-0.5, 0.30901699, 0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.5, 0.30901699, 0.80901699, 0.],
             [0.30901699, -0.80901699, 0.5, 0.],
             [0.80901699, 0.5, 0.30901699, 0.],
             [0., 0., 0., 1.]],

            [[0.80901699, 0.5, 0.30901699, 0.],
             [-0.5, 0.30901699, 0.80901699, 0.],
             [0.30901699, -0.80901699, 0.5, 0.],
             [0., 0., 0., 1.]],

            [[0.30901699, 0.80901699, -0.5, 0.],
             [-0.80901699, 0.5, 0.30901699, 0.],
             [0.5, 0.30901699, 0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[0.30901699, -0.80901699, -0.5, 0.],
             [0.80901699, 0.5, -0.30901699, 0.],
             [0.5, -0.30901699, 0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[-1., 0., 0., 0.],
             [0., -1., 0., 0.],
             [0., 0., 1., 0.],
             [0., 0., 0., 1.]],

            [[-0.30901699, -0.80901699, 0.5, 0.],
             [0.80901699, -0.5, -0.30901699, 0.],
             [0.5, 0.30901699, 0.80901699, 0.],
             [0., 0., 0., 1.]],

            [[-0.30901699, 0.80901699, 0.5, 0.],
             [-0.80901699, -0.5, 0.30901699, 0.],
             [0.5, -0.30901699, 0.80901699, 0.],
             [0., 0., 0., 1.]]
        ]

        for m1, m2 in izip(matrices[:len(refMatrices)], refMatrices):
            self.assertArrayAlmostEqual(m1, m2)

    def test_42_SymmetryIcosahedral222UnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_I222,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [[-0.30901699437494745, -0.8090169943749475, -0.5],
                          [-0.3090169943749475, -0.8090169943749475, 0.5],
                          [0.0, 1.0, 1.0558848421106697e-16],
                          ]
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_43_SymmetryIcosahedral222rUnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_I222r,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [ [-0.5, -0.8090169943749476, 0.30901699437494734],
                           [0.5, -0.8090169943749475, 0.3090169943749474],
                           [0.0, 1.0, -0.0]
                          ]
                         
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_44_SymmetryIcosahedraln25UnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_In25,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [ 
                           [-0.5877852522924731, -0.8090169943749475, -0.0],
                           [0.26286555605956674, -0.8090169943749477, 0.5257311121191334],
                           [0.0, 1.0, 0.0]
                          ]
                         
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)
    
    def test_45_SymmetryIcosahedraln25RUnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_In25r,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [ 
                           [0.16245984811645317, -0.5000000000000001, 0.8506508083520399],
                           [0.1624598481164532, 0.5000000000000001, 0.8506508083520399],
                           [-0.5257311121191337, 6.206335383118184e-17, -0.8506508083520399]
                          ]
                         
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)
    
    def test_46_SymmetryIcosahedral2n3UnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_I2n3,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [ 
                           [0.3090169943749475, -0.17841104488654505, -0.9341723589627157],
                           [0.3090169943749475, 0.7557613140761706, -0.577350269189626],
                           [0.0, -0.3568220897730899, 0.9341723589627158]
                          ]
                         
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_47_SymmetryIcosahedral2n3rUnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_I2n3r,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [ 
                           [0.3090169943749475, -0.7557613140761706, -0.577350269189626],
                           [0.3090169943749475, 0.17841104488654505, -0.9341723589627157],
                           [0.0, 0.3568220897730899, 0.9341723589627158]    
                          ]
                         
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_48_SymmetryIcosahedral2n5UnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_I2n5,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [ 
                           [0.3090169943749475, 0.42532540417602, -0.8506508083520399],
                           [0.3090169943749475, 0.9510565162951535, 0.],
                           [0.0, -0.8506508083520401, 0.5257311121191336]
                          ]
                         
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_49_SymmetryIcosahedral2n5rUnitCell(self):
        vectorsEdge, vectorsPlane = emconv.getUnitCell(sym=emcts.SYM_I2n5r,
                                                       circumscribed_radius=1, 
                                                       center=(0, 0, 0))
        print("vectorsPlane", vectorsPlane)
        refereceVectors= [ 
                           [0.3090169943749475, -0.9510565162951535, 0.],
                           [0.3090169943749475, -0.42532540417602, -0.8506508083520399],
                           [0.0, 0.8506508083520401, 0.5257311121191336]
                          ]
                         
        for v, w in zip(refereceVectors, vectorsPlane):
            self.assertArrayAlmostEqual(v, w, decimal=3)

    def test_50_moveParticlesInsideUnitCell(self):
        """ Test moveParticlesInsideUnitCell function. """
        import pwem.objects as emobj
        transforms =[
            [[-0.8185494886822994, -0.3483236755177939, 0.456779324895207, -2.156149758375051], #0
             [0.5364525666015211, -0.1791764432657179, 0.8246905152633245, 0.6413349977154054],
             [-0.20541513664914549, 0.9200904408004801, 0.33352391575866236, 1.367163551807022],
             [0.0, 0.0, 0.0, 1.0]],
            [[-0.858752733295496, -0.2365739023950405, 0.4545069105779325, -8.985651938689228], #1
             [-0.43079706089217373, -0.14689274461793103, -0.8904136195641194, -4.780540644113263],
              [0.2774123922685857, -0.9604453707963785, 0.024229204106350065, -7.124360976049757],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.13871745845956557, 0.829766457494788, -0.5405969781039005, -2.5924361521394617],   #2 
             [0.7771818232984424, -0.24711551312977775, -0.5787247503823518, -0.8228381223887112],
              [-0.613796285629872, -0.5004213716331039, -0.6106001724203945, 3.075105480752516],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.5608446409022507, -0.16808852524103843, 0.8106784420801619, -3.334470959628195], #3
             [-0.1193335532577222, -0.9525367503974368, -0.28005935479676053, -1.8257001199859955],
              [0.8192757727640382, -0.25381092731516547, 0.5141665307439837, 4.262975503064537],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.12030913630850786, -0.800596293924665, -0.5870019470792173, -5.928069556967794], #4
             [-0.7650422218167392, -0.3020459445010162, 0.5687518318635961, -20.73792550983718],
              [-0.632642166282337, 0.5175073154696272, -0.5761510807739681, -14.431882587531991],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.3812631664773085, -0.0024653097546033565, -0.9244632605655619, 0.8811534929936544], #5
             [0.24984464412018964, 0.9625096854638254, -0.10560662475803186, -2.559716603556499],
              [0.8900651951919882, -0.27123611049442115, -0.3663535460067314, -1.485644144257723],
               [0.0, 0.0, 0.0, 1.0]],
            [[0.7813806932698447, -0.28648187016239474, -0.5544117154727108, -4.581844214196944], #6
             [0.16569704928607887, -0.7612721837777564, 0.6269044186031334, -5.881947754563014],
              [-0.6016549676043771, -0.5817153945655064, -0.5473741861674454, -2.0958075033442602],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.994625504490887, 0.01788179047098669, 0.10198209345664482, 9.183780304028408], #7
             [-0.066294630386591, 0.6465989370701519, -0.759943969356725, -3.109049922358411],
              [-0.07953067205898232, -0.7626205190979988, -0.6419383272967468, 5.139566746340785],
               [0.0, 0.0, 0.0, 1.0]],
            [[0.11819744348090917, 0.7176250505080871, -0.6863261988572509, -1.0799299361841568], #8
             [0.6544950244545856, -0.5760914836805305, -0.4896476951798675, 4.094478985594572],
              [-0.7467701301731136, -0.39132197652831935, -0.5377745655636763, -2.349809885115758],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.6223378990980357, 0.2839941514768155, 0.7294126824186749, -2.0435766188882525], #9
             [-0.000998353913662789, -0.9321478333605614, 0.36207653893986863, -7.478677564449396],
              [0.7827480709881733, 0.22460574055035804, 0.5803944509349818, 7.233812082336363],
               [0.0, 0.0, 0.0, 1.0]]
        ]
        transformsOUT = [
            [[-0.7252858385122594, -0.6864764274765147, -0.0520631056761925, -2.626126986227842],
             [0.6873563224915256, -0.7263190893242406, 0.0013661676806916044, -0.11910413088886435],
              [-0.038752269410898904, -0.03479504278323646, 0.9986428623953699, -0.13674719534996074],
               [0.0, 0.0, 0.0, 1.0]],
            [[0.9665753642787465, -0.24343802349768162, -0.0804362721970048, 5.184632937244055],
             [0.25633637920483676, 0.9235683557176511, 0.2851546124704075, 10.930718580335737],
              [0.0048709204020018545, -0.29624216619606375, 0.9551004413682893, -2.8268359616178897],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.88778551534138, 0.46016517026001924, 0.00921383902659667, 0.31973180161037973],
             [-0.43483575262317287, -0.845142958161346, 0.3108878391171456, 3.8891486077267636],
              [0.1508467665900711, 0.27199521383630415, 0.9504019447893948, -1.2752867043868823],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.7764517718598236, 0.6296635178958486, -0.025426368336134816, -3.5663777623651547],
             [-0.630135329899526, -0.7762307915602822, 0.019880247637118934, -4.4615897805661815],
              [-0.007218863356232397, 0.03145810650249034, 0.9994790020541313, -0.001983558348627401],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.8571618070910606, 0.5150471011340995, -0.00034652236317347316, -20.212726792116044],
             [-0.4835895462789405, -0.8045770043391588, 0.34466939930657253, -15.281671963652643],
              [0.1772421710375828, 0.2956050197509946, 0.9387240729329882, 5.603478898568271],
               [0.0, 0.0, 0.0, 1.0]],
            [[0.7129723447225714, -0.6999276429994896, 0.04209192590007993, -0.19434446250775977],
             [0.6678054177014275, 0.6961027800370845, 0.26358460446935095, -2.970520267727017],
              [-0.21379045757334622, -0.15981931732372437, 0.9637175032449903, 0.8209498529268272],
               [0.0, 0.0, 0.0, 1.0]],
            [[0.9841809478811525, -0.17615704639119623, 0.01888271258284585, -4.476507899377257],
             [0.16244065560214985, 0.9397814984423979, 0.3007054515516885, 6.052385492190926],
              [-0.07071700810484034, -0.29288125612966626, 0.9535301120429331, -1.820033450896662],
               [0.0, 0.0, 0.0, 1.0]],
            [[0.5263699056237151, -0.7677131245303491, 0.3654466867787004, -0.4884024601224816],
             [0.8239481047761119, 0.5666536472675733, 0.00362831561603133, -10.960366974104154],
              [-0.20986720346336674, 0.299199268819636, 0.9308252008020953, 0.23447230621585655],
               [0.0, 0.0, 0.0, 1.0]],
            [[0.9679222571562066, 0.2502987236397183, -0.02184154403796548, 3.61455892039102],
             [-0.23963421466206813, 0.9458054668859154, 0.21915168713019173, -3.12633686249817],
              [0.07551123932851517, -0.2068878144141011, 0.9754463004092169, 0.7833230240052959],
               [0.0, 0.0, 0.0, 1.0]],
            [[-0.8260686557653656, -0.5600247587102615, -0.06310978999739547, -10.223115595493407],
             [0.5522437205160172, -0.8267135233326576, 0.10757148083730465, -2.7932180680893124],
              [-0.11241640944554798, 0.05400944334481576, 0.9921922852536971, -0.3474206356852143],
               [0.0, 0.0, 0.0, 1.0]]
        ]
        # TODO: move database to memory. 
        testFileIN = ':memory:'
        testFileOUT = ':memory:'
        imgSetIN  = emobj.SetOfParticles(filename=testFileIN)
        imgSetOUT = emobj.SetOfParticles(filename=testFileOUT)
        img    = emobj.Particle()
        t      = emobj.Transform()
        
        for i in range(1, 11):
            t.setMatrix(np.array(transforms[i-1]))
            img.setTransform(t)

            img.setLocation(i, 'mystack.stk')
            img.setSamplingRate(1.5)

            #img.setTransform(t)
            imgSetIN.append(img)
            img.setObjId(None)
            t.setObjId(None)

        emconv.moveParticlesInsideUnitCell(setIN=imgSetIN, 
                                           setOUT=imgSetOUT,
                                           sym=emcts.SYM_I222r)
        for i, img in enumerate(imgSetOUT):
            print(img.getTransform().getMatrix(), transformsOUT[i])
            self.assertTrue(np.allclose(img.getTransform().getMatrix(), transformsOUT[i]))

