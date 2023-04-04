# **************************************************************************
# *
# * Authors:     Roberto Marabini
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
"""
This module returns the matrices related with the different
point symmetries. Code based on chimera file sym.py
"""

from math import sin, cos, pi, acos, sqrt, asin
import numpy as np
from numpy.linalg import inv as matrix_inverse
from numpy import dot as matrix_multiply
import operator

import pwem.constants as cts

DEBUG = False  # set to True for debuging
if DEBUG:
    from pwem.viewers.viewer_chimera import Chimera
    # expansion factor is use to create chimera bild files
    # this is the factor to multiply the coordinates
    expansionFactor = 100


def _column(matrix, i):
    """Returns column i of matrix as a list."""
    return [row[i] for row in matrix]


def _applyMatrix(tf, points):
    """
    Args:multiply point by a matrice list
    """
    tf = np.array(tf)
    r = matrix_multiply(points, np.transpose(tf[:, :3]))
    np.add(r, tf[:, 3], r)
    return r


def __length(v):
    """Returns the length of a vector."""
    d = sqrt(sum([e*e for e in v]))
    return d


def _normalizeVector(v):
    """Returns a unit vector in the same direction as v."""
    d = __length(v)
    if d == 0:
        d = 1
    return tuple([e/d for e in v])


def _reflectMatrix(axis):
    """Returns a matrix that reflects points in the plane
    perpendicular to axis."""
    axis = _normalizeVector(axis)
    return ((1 - 2*axis[0]*axis[0], -2*axis[0]*axis[1], -2*axis[0]*axis[2], 0),
            (-2*axis[1]*axis[0], 1 - 2*axis[1]*axis[1], -2*axis[1]*axis[2], 0),
            (-2*axis[2]*axis[0], -2*axis[2]*axis[1], 1 - 2*axis[2]*axis[2], 0))


def _rotationTransform(axis, angle, center=(0, 0, 0)):
    """ returns matrix that rotates "angle" angles arround vector axis
        angle is in degrees
    """
    axis = _normalizeVector(axis)

    arad = angle*pi/180.0
    sa = sin(arad)
    ca = cos(arad)
    if axis == (0, 0, 1) and center == (0, 0, 0):
        return ((ca, -sa, 0, 0),
                (sa, ca, 0, 0),
                (0, 0, 1, 0))

    k = 1 - ca
    ax, ay, az = axis
    tf = ((1 + k*(ax*ax-1), -az*sa+k*ax*ay, ay*sa+k*ax*az, 0),
          (az*sa+k*ax*ay, 1 + k*(ay*ay-1), -ax*sa+k*ay*az, 0),
          (-ay*sa+k*ax*az, ax*sa+k*ay*az, 1 + k*(az*az-1), 0))

    if center == (0, 0, 0):
        return tf

    cx, cy, cz = center
    c_tf = ((1, 0, 0, cx), (0, 1, 0, cy), (0, 0, 1, cz))
    inv_c_tf = ((1, 0, 0, -cx), (0, 1, 0, -cy), (0, 0, 1, -cz))
    rtf = _multiplyMatrices(c_tf, tf, inv_c_tf)
    return rtf


def _translationMatrix(shift):
    """ returns matriz that shift point by vector shift"""
    tf = np.array(((1.0, 0, 0, shift[0]),
                   (0, 1.0, 0, shift[1]),
                   (0, 0, 1.0, shift[2])))
    return tf


def _identityMatrix():
    """ returns identity matrix"""
    return (1.0, 0, 0, 0), (0, 1.0, 0, 0), (0, 0, 1.0, 0)


def _invertMatrix(tf):
    """ returns inverse 3x3 matrix"""
    tf = np.array(tf)
    r = tf[:, :3]
    t = tf[:, 3]
    tfinv = np.zeros((3, 4), np.float)
    rinv = tfinv[:, :3]
    tinv = tfinv[:, 3]
    rinv[:, :] = matrix_inverse(r)
    tinv[:] = matrix_multiply(rinv, -t)
    return tfinv


def _multiplyMatrices(*mlist):
    """ returns matrix product of matrices in mlist"""
    if len(mlist) == 2:
        m1, m2 = mlist
        p = [[0, 0, 0, 0],
             [0, 0, 0, 0],
             [0, 0, 0, 0]]
        for r in range(3):
            for c in range(4):
                p[r][c] = m1[r][0]*m2[0][c] +\
                          m1[r][1]*m2[1][c] +\
                          m1[r][2]*m2[2][c]
            p[r][3] += m1[r][3]
        p = tuple(map(tuple, p))
    else:
        p = _multiplyMatrices(*mlist[1:])
        p = _multiplyMatrices(mlist[0], p)
    return p


def _matrixProducts(mlist1, mlist2):
    """ returns list of matrix products of matrices in mlist1 and mlist2"""
    plist = []
    for m1 in mlist1:
        for m2 in mlist2:
            m1xm2 = _multiplyMatrices(m1, m2)
            plist.append(m1xm2)
    return plist


def _coordinateTransformList(tflist, ctf):
    """ returns list of coordinate transforms of tflist by ctf"""
    ctfinv = _invertMatrix(ctf)
    return [_multiplyMatrices(ctfinv, tf, ctf) for tf in tflist]


def _recenterSymmetries(tflist, center):

    if center == (0, 0, 0):
        return tflist
    ctf = _translationMatrix([-x for x in center])
    return _coordinateTransformList(tflist, ctf)


def _transposeMatrix(tf):
    """ returns transpose of 3x3 matrix"""
    return ((tf[0][0], tf[1][0], tf[2][0], tf[0][3]),
            (tf[0][1], tf[1][1], tf[2][1], tf[1][3]),
            (tf[0][2], tf[1][2], tf[2][2], tf[2][3]))

# ==================== End of utility functions ====================


def moveParticlesInsideUnitCell(setIN, setOUT, sym=cts.SYM_CYCLIC, n=1):
    """ Move particles projection direction inside the unit cell associated
    to symmetry sym.
    This function (1) gets the symmetry matrices and (2) applies them to
         each proejction direction until is inside the unit cell.
    """
    matrixSet = getSymmetryMatrices(sym=sym, n=n)
    # get unit cell vectors
    _, (u, v, w) = getUnitCell(sym=sym, n=n, generalize=False)
    counter = 0  # counter for the number of particles
    totalNumberOfParticles = len(setIN)
    totalNumberOFMatrices = len(matrixSet)
    # for each particle in set
    for par in setIN:
        counter += 1
        if counter % 1000 == 0:
            # print something to keep the user engaged
            print(f"{counter}/{totalNumberOfParticles}")
        # get projection direction
        colunm = par.getTransform().getMatrix()[0:3, 2]
        if (np.dot(colunm, u) > 0) and \
           (np.dot(colunm, v) > 0) and \
           (np.dot(colunm, w) > 0):
            pass  # particle is inside unit cell so nothing needs to be done
        else:
            # loop over symmetry matrices
            for i, matrix in enumerate(matrixSet):
                matrix3 = np.array(matrix)
                matrix3 = np.delete(matrix3, -1, 1)
                columPrime = matrix3.dot(colunm)[:3]
                if (np.dot(columPrime, u) > 0) and \
                   (np.dot(columPrime, v) > 0) and \
                   (np.dot(columPrime, w) > 0):
                    matrix = np.dot(matrix, par.getTransform().getMatrix())
                    t = par.getTransform()
                    t.setMatrix(matrix)
                    par.setTransform(t)
                    setOUT.append(par)
                    break
                if i == totalNumberOFMatrices-1:
                    print("Error: something went wrong in "
                          "moveParticlesInsideUnitCell")
                    print("       no matrix move particle direction "
                          "inside the unit cell")
                    print("       particle id: ", par.getObjId())


def getSymmetryMatrices(sym=cts.SYM_CYCLIC,
                        n=1,
                        circumscribed_radius=1,
                        center=(0, 0, 0),
                        offset=None):
    """ interface between scipion and chimera code
        chimera code uses tuples of tuples as matrices
        but scipion uses np.arrays (lists of lists)
        so let us convert them here.
        Direct use of the classes return 3x3 arrays
        which may be OK in many cases, but not in all
    """
    if sym == cts.SYM_CYCLIC:
        c = Cyclic(n=n,
                   center=center,
                   circumscribed_radius=circumscribed_radius,
                   offset=offset)
        matrices = c.symmetryMatrices()
    elif sym == cts.SYM_DIHEDRAL_X or sym == cts.SYM_DIHEDRAL_Y:
        d = Dihedral(n=n,
                     center=center,
                     circumscribed_radius=circumscribed_radius,
                     offset=offset, sym=sym)
        matrices = d.symmetryMatrices()
    elif sym == cts.SYM_OCTAHEDRAL:
        o = Octahedral(center=center,
                       circumscribed_radius=circumscribed_radius,
                       offset=offset)
        matrices = o.symmetryMatrices()
    elif (sym == cts.SYM_TETRAHEDRAL_222 or
          sym == cts.SYM_TETRAHEDRAL_Z3 or
          sym == cts.SYM_TETRAHEDRAL_Z3R):

        t = Tetrahedral(center=center,
                        circumscribed_radius=circumscribed_radius,
                        sym=sym)
        matrices = t.symmetryMatrices()
    elif (sym == cts.SYM_I222 or sym == cts.SYM_I222r or
          sym == cts.SYM_In25 or sym == cts.SYM_In25r or
          sym == cts.SYM_I2n3 or sym == cts.SYM_I2n3r or
          sym == cts.SYM_I2n5 or sym == cts.SYM_I2n5r):
        matrices = __icosahedralSymmetryMatrices(sym, center)

        # '222'         2-fold symmetry along x, y, and z axes.
        # '222r'        '222' with 90 degree rotation around z.
        # '2n5'         2-fold symmetry along x and 5-fold along z.
        # '2n5r'        '2n5' with 180 degree rotation about y.
        # 'n25'         2-fold symmetry along y and 5-fold along z.
        # 'n25r'        'n25' with 180 degree rotation about x.
        # '2n3'         2-fold symmetry along x and 3-fold along z.
        # '2n3r'        '2n3' with 180 degree rotation about y.

    # convert from 4x3 to 4x4 matrix, Scipion standard
    extraRow = (0., 0., 0., 1.)
    for i in range(len(matrices)):
        matrices[i] += (extraRow,)
    # convert from sets to lists Scipion standard
    return np.array(matrices)


def getUnitCell(sym=cts.SYM_CYCLIC,
                n=1,
                circumscribed_radius=1,
                center=(0, 0, 0),
                offset=None,
                generalize=False):
    """ return vectors normal to the unit cell faces
    """

    if sym == cts.SYM_CYCLIC:
        c = Cyclic(n=n,
                   center=center,
                   circumscribed_radius=circumscribed_radius,
                   offset=offset)
        vectorsEdge, vectorsPlane = c.unitCellPlanes()
        if DEBUG:
            bildFileName = f'/tmp/C{n}.bild'
            vLabels = ['v1', 'v2', 'v3']
            pLabels = ['p1', 'p2']
    elif sym == cts.SYM_DIHEDRAL_X or sym == cts.SYM_DIHEDRAL_Y:
        print("getUnitCell:", cts.SCIPION_SYM_NAME[sym])
        d = Dihedral(n=n,
                     center=center,
                     circumscribed_radius=circumscribed_radius,
                     offset=offset, sym=sym)
        vectorsEdge, vectorsPlane = d.unitCellPlanes()
        if DEBUG:
            vLabels = ['v1', 'v2', 'eigenvector']
            pLabels = ['p1', 'p2', 'p3']
    elif sym == cts.SYM_OCTAHEDRAL:
        o = Octahedral(center=center,
                       circumscribed_radius=circumscribed_radius)
        vectorsEdge, vectorsPlane = o.unitCellPlanes()
        if DEBUG:
            vLabels = ['v1', 'v2', 'v3']
            pLabels = ['p1', 'p2', 'p3']
    elif (sym == cts.SYM_TETRAHEDRAL_Z3 or
          sym == cts.SYM_TETRAHEDRAL or
          sym == cts.SYM_TETRAHEDRAL_Z3R):
        t = Tetrahedral(center=center,
                        circumscribed_radius=circumscribed_radius,
                        sym=sym)
        vectorsEdge, vectorsPlane = t.unitCellPlanes()
        if DEBUG:
            vLabels = ['v1', 'v2', 'v3']
            pLabels = ['p1', 'p2', 'p3']
    elif (sym == cts.SYM_I222 or sym == cts.SYM_I222r or
          sym == cts.SYM_In25 or sym == cts.SYM_In25r or
          sym == cts.SYM_I2n3 or sym == cts.SYM_I2n3r or
          sym == cts.SYM_I2n5 or sym == cts.SYM_I2n5r):
        i_sym = cts.SCIPION_SYM_NAME[sym][1:]
        i = Icosahedron(center=center,
                        circumscribed_radius=circumscribed_radius,
                        orientation=i_sym)
        vectorsEdge, vectorsPlane = i.unitCellPlanes()
        if DEBUG:
            vLabels = ['v1', 'v2', 'v3']
            pLabels = ['p1', 'p2', 'p3']

    if DEBUG:
        bildFileName = cts.SCIPION_SYM_NAME[sym] + '.bild'
        Chimera.createCoordinateAxisFile(dim=expansionFactor,
                                         bildFileName=bildFileName,
                                         sampling=1)
        colors = ['cyan', 'green', 'magenta', 'navy blue']
        with open(bildFileName, 'a') as f:
            for i, vector in enumerate(vectorsEdge):
                f.write(f'.comment {vLabels[i]}\n')
                f.write(f'.color {colors[i]}\n')
                x = vector[0]*expansionFactor
                y = vector[1]*expansionFactor
                z = vector[2]*expansionFactor
                f.write(f'.arrow 0 0 0  {x} {y} {z} 0.20 0.40 0.75\n')

            for i, vector in enumerate(vectorsPlane):
                f.write(f'.comment {pLabels[i]}\n')
                f.write(f'.color {colors[i]}\n')
                x = vector[0]
                y = vector[1]
                z = vector[2]
                r = expansionFactor/2
                f.write(f'.cone 0 0 0  {x/(expansionFactor*100)} \
                                       {y/(expansionFactor*100)} \
                                       {z/(expansionFactor*100)} {r}\n')
                f.write(f'.arrow 0 0 0  {x*expansionFactor} \
                                        {y*expansionFactor} \
                                        {z*expansionFactor} 0.20 0.40 0.75\n')
        f.close()
    if generalize:
        for i, v in enumerate(vectorsEdge):
            vectorsEdge[i] = np.append(v, 1.)
        for i, v in enumerate(vectorsPlane):
            vectorsPlane[i] = np.append(v, 1.)
    return vectorsEdge, vectorsPlane


def __icosahedralSymmetryMatrices(orientation=cts.SYM_I222, center=(0, 0, 0)):
    if orientation == cts.SYM_I222:
        sym = '222'
    elif orientation == cts.SYM_I222r:
        sym = '222r'
    elif orientation == cts.SYM_In25:
        sym = 'n25'
    elif orientation == cts.SYM_In25r:
        sym = 'n25r'
    elif orientation == cts.SYM_I2n5:
        sym = '2n5'
    elif orientation == cts.SYM_I2n5r:
        sym = '2n5r'
    elif orientation == cts.SYM_I2n3:
        sym = '2n3'
    elif orientation == cts.SYM_I2n3r:
        sym = '2n3r'

    i = Icosahedron(orientation=sym, center=center)
    return list(i.icosahedralSymmetryMatrices())


def __icosahedralUnitCellPlanes(orientation=cts.SYM_I222, center=(0, 0, 0)):
    if orientation == cts.SYM_I222:
        sym = '222'
    elif orientation == cts.SYM_I222r:
        sym = '222r'
    elif orientation == cts.SYM_In25:
        sym = 'n25'
    elif orientation == cts.SYM_In25r:
        sym = 'n25r'
    elif orientation == cts.SYM_I2n5:
        sym = '2n5'
    elif orientation == cts.SYM_I2n5r:
        sym = '2n5r'
    elif orientation == cts.SYM_I2n3:
        sym = '2n3'
    elif orientation == cts.SYM_I2n3r:
        sym = '2n3r'

    i = Icosahedron(orientation=sym, center=center)
    return list(i.unitCellPlanes())


class Cyclic(object):
    """cyclic class.
        Allow to compute symmetry matrices, vertices,
        symmetry axis and unit cell planes
    """

    def __init__(self, n, circumscribed_radius=1, center=(
            0, 0, 0), offset=None):
        """
        :Parameters:
            n: int order of symmetry
            circumscribed_radius: float radius of the circumscribed sphere
            center: tuple of 3 floats, center of the circumscribed sphere
            offset: float, angle in degrees to rotate the symmetry axis
                    it modifies the unit cell but not the symmetry matrices
            """
        self.n = n
        self.circumscribed_radius = circumscribed_radius
        self.center = center
        self.matrices = None
        self.offset = offset

    def symmetryMatrices(self):
        """ get Matrices for cyclic symmetry of order n
        """
        if self.matrices is not None:
            return self.matrices
        self.matrices = []
        for k in range(self.n):
            a = 2*pi * np.float32(k) / self.n
            c = cos(a)
            s = sin(a)
            tf = ((c, -s, 0, 0),
                  (s, c, 0, 0),
                  (0, 0, 1, 0))
            self.matrices.append(tf)
        self.matrices = _recenterSymmetries(self.matrices, self.center)
        return self.matrices

    def unitCellPlanes(self):
        """ get planes that define a unit cell for cyclic symmetry of order n
        """
        if self.offset is None:
            self.offset = pi/(self.n)
        _ = self.symmetryMatrices()
        # these three vectors are enges of the unit cell
        angle = -2*pi/self.n
        a1 = self.offset
        a2 = angle + self.offset
        c1 = cos(a1)
        s1 = sin(a1)
        c2 = cos(a2)
        s2 = sin(a2)
        v1 = np.array([c1, s1, 0.])
        v2 = np.array([c2, s2, 0.])
        v3 = np.array([0., 0., 1.])  # keep this dim 3 instead of dim 4

        # cross product of v1/v2 with eigenvetor
        # these two vectors are normal to the planes that define the unit cell
        plane1 = np.cross(v1, v3)
        plane2 = np.cross(v3, v2)
        if DEBUG:
            print("v1", v1)
            print("v2", v2)
            print("v3", v3)
            print("plane1", plane1)
            print("plane2", plane2)
        c = self.circumscribed_radius
        return [v1*c, v2*c, v3*c], [plane1*c, plane2*c]


class Dihedral(Cyclic):
    """dihedral class.
        Allow to compute symmetry matrices,  vertices,
        symmetry axis and unit cell planes
    """
    def __init__(self, n, circumscribed_radius=1, center=(
            0, 0, 0), offset=None, sym=cts.SYM_DIHEDRAL_X):
        """
        :Parameters:
            n: int order of symmetry
            circumscribed_radius: float radius of the circumscribed sphere
            center: tuple of 3 floats, center of the circumscribed sphere
            offset: float, angle in degrees to rotate the symmetry axis
                    it modifies the unit cell but not the symmetry matrices
            """
        super().__init__(n=n,
                         circumscribed_radius=circumscribed_radius,
                         center=center,
                         offset=offset)
        self.sym = sym

    def symmetryMatrices(self):
        clist = super().symmetryMatrices()
        if self.sym == cts.SYM_DIHEDRAL_X:
            # matrix that reflect in the x axis
            reflect = ((1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0))
        else:
            print("Be carefull untested code")
            # matrix that reflect in the y axis
            reflect = ((-1, 0, 0, 0), (0, 1, 0, 0), (0, 0, -1, 0))

        tflist = _matrixProducts([_identityMatrix(), reflect], clist)
        tflist = _recenterSymmetries(tflist, self.center)
        return tflist

    def unitCellPlanes(self):
        """ get planes that define a unit cell for dihedral symmetry of order n
        """
        if self.offset is None:
            self.offset = pi/(self.n) + pi/2
        vectorsEdge, vectorsPlane = super().unitCellPlanes()
        vectorsPlane.append(np.array([0., 0., self.circumscribed_radius]))
        if DEBUG:
            print("v1", vectorsEdge[0])
            print("v2", vectorsEdge[1])
            print("eigenvector", vectorsEdge[2])
            print("plane1", vectorsPlane[0])
            print("plane2", vectorsPlane[1])
            print("plane3", vectorsPlane[2])

        return vectorsEdge, vectorsPlane

    def coordinateSystemTransform(self, origSym, targetSym):
        print("========================",
              cts.SCIPION_SYM_NAME[origSym],
              "->",
              cts.SCIPION_SYM_NAME[targetSym])
        if origSym == targetSym:
            return np.identity(4)
        elif origSym == cts.SYM_DIHEDRAL_X and targetSym == cts.SYM_DIHEDRAL_Y:
            return ((0., -1., 0., 0.), (1, 0., 0., 0.), (0., 0., 1., 0.))
        elif origSym == cts.SYM_DIHEDRAL_Y and targetSym == cts.SYM_DIHEDRAL_X:
            return ((0., 1., 0., 0.), (-1, 0., 0., 0.), (0., 0., 1., 0.))
        else:
            raise Exception("unknown dihedral symmetry pair: %s %s" % (
                cts.SCIPION_SYM_NAME[origSym],
                cts.SCIPION_SYM_NAME[targetSym]))


class Tetrahedral(object):
    """Tetrahedral class.
        Allow to compute symmetry matrices,  vertices,
        symmetry axis and unit cell planes
    """

    def __init__(self, circumscribed_radius=1, center=(
            0, 0, 0), sym=cts.SYM_TETRAHEDRAL_Z3):
        """
        :Parameters:
            n: int order of symmetry
            circumscribed_radius: float radius of the circumscribed sphere
            center: tuple of 3 floats, center of the circumscribed sphere
            """
        self.circumscribed_radius = circumscribed_radius
        self.center = center
        self.matrices = None
        self.sym = sym

    def symmetryMatrices(self):
        """
        identity
        4 * rotation by 120 clockwise (seen from a vertex):
                            (234), (143), (412), (321)
        4 * rotation by 120 counterclockwise (ditto)
        3 * rotation by 180
         """
        # simmetries are rotations around theses axis by these angles
        self.aa = (((0, 0, 1), 0), ((1, 0, 0), 180),  # 0 1
                   ((0, 1, 0), 180), ((0, 0, 1), 180),  # 2 3
                   ((1, 1, 1), 120), ((1, 1, 1), 240),  # 4 5
                   ((-1, -1, 1), 120), ((-1, -1, 1), 240),  # 6 7
                   ((-1, 1, -1), 120), ((-1, 1, -1), 240),  # 8 9
                   ((1, -1, -1), 120),  ((1, -1, -1), 240))  # 10 11

        syms = [_rotationTransform(axis, angle) for axis, angle in self.aa]

        if self.sym == cts.SYM_TETRAHEDRAL_Z3R:
            # convention, 3-fold on z, 3-fold in yz plane sign(y) = sign(z)
            self.tf = _multiplyMatrices(
                     _rotationTransform((0, 0, 1), -45.0),
                     _rotationTransform((1, 0, 0),
                                        -acos(1 / sqrt(3)) * 180 / pi))
            syms = _coordinateTransformList(syms, self.tf)
        elif self.sym == cts.SYM_TETRAHEDRAL_Z3:
            # convention, 3-fold on z, 3-fold in yz plane sign(y) != sign(z)
            tfaux = _multiplyMatrices(
                     _rotationTransform((0, 0, 1), -45.0),
                     _rotationTransform((1, 0, 0),
                                        -acos(1 / sqrt(3)) * 180 / pi)
            )
            self.tf = _multiplyMatrices(
                     tfaux,
                     _reflectMatrix((0, 0, 1))
            )
            syms = _coordinateTransformList(syms, self.tf)

        syms = _recenterSymmetries(syms, self.center)
        return syms

    def unitCellPlanes(self):
        """ get planes that define a unit cell for tetrahedral symmetry of
            order n
        """
        print("Unit cell sym", cts.SCIPION_SYM_NAME[self.sym])
        _ = self.symmetryMatrices()
        # four corners tetahedron
        v0 = _normalizeVector(np.array(self.aa[4][0])) *\
            self.circumscribed_radius
        v1 = _normalizeVector(np.array(self.aa[6][0])) *\
            self.circumscribed_radius
        v2 = _normalizeVector(np.array(self.aa[8][0])) *\
            self.circumscribed_radius
        v3 = _normalizeVector(np.array(self.aa[10][0])) *\
            self.circumscribed_radius

        if self.sym == cts.SYM_TETRAHEDRAL_Z3 or\
           self.sym == cts.SYM_TETRAHEDRAL_Z3R:
            self.tf = _invertMatrix(self.tf)
            if self.sym == cts.SYM_TETRAHEDRAL_Z3:
                _3fold_1 = np.dot(self.tf, np.append(v0, 1))
                _3fold_2 = np.dot(self.tf, np.append(v1, 1))
                _3fold_3 = np.dot(self.tf, np.append(v3, 1))
            else:
                _3fold_1 = np.dot(self.tf, np.append(v1, 1))
                _3fold_2 = np.dot(self.tf, np.append(v0, 1))
                _3fold_3 = np.dot(self.tf, np.append(v3, 1))
        else:
            _3fold_1 = v1
            _3fold_2 = v0
            _3fold_3 = v2

        # colors = ['cyan', 'green', 'magenta', 'navy blue']
        vp = np.add(_3fold_1, _3fold_2)
        vp = _normalizeVector(np.add(vp, _3fold_3))

        plane1 = _normalizeVector(np.cross(_3fold_1, _3fold_2))  # cyan
        plane2 = _normalizeVector(np.cross(_3fold_2, vp))  # green
        plane3 = _normalizeVector(np.cross(vp, _3fold_1))  # magenta
        if DEBUG:
            print("v0", _3fold_1)
            print("v1", _3fold_2)
            print("v2", vp)

            print("plane1", plane1)
            print("plane2", plane2)
            print("plane3", plane3)
        return [_3fold_1, vp, _3fold_2], [plane1, plane2, plane3]

    def coordinateSystemTransform(self, origSym, targetSym):
        if origSym == targetSym:
            return np.identity(4)

        elif (origSym == cts.SYM_TETRAHEDRAL_222 and
              targetSym == cts.SYM_TETRAHEDRAL_Z3 or
              origSym == cts.SYM_TETRAHEDRAL_222 and
              targetSym == cts.SYM_TETRAHEDRAL_Z3R):
            print("========================222 -> ")
            matrix = _multiplyMatrices(
                     _rotationTransform((1, 0, 0),
                                        -acos(1 / sqrt(3)) * 180 / pi),
                     _rotationTransform((0, 0, 1),  45.0)
            )
            # matrix = _multiplyMatrices(
            #          matrix,
            #          _reflectMatrix((0, 0, 1))
            # )
            # matrix = _invertMatrix(matrix)
            if targetSym == cts.SYM_TETRAHEDRAL_Z3:
                return matrix
            else:
                matrix2 = _reflectMatrix((1, 0, 0))
                return _multiplyMatrices(matrix2, matrix)

        elif (origSym == cts.SYM_TETRAHEDRAL_Z3R and
              targetSym == cts.SYM_TETRAHEDRAL_Z3) or\
             (origSym == cts.SYM_TETRAHEDRAL_Z3 and
              targetSym == cts.SYM_TETRAHEDRAL_Z3R):
            matrix = _reflectMatrix((1, 0, 0))
            return matrix

        elif (origSym == cts.SYM_TETRAHEDRAL_Z3 and
              targetSym == cts.SYM_TETRAHEDRAL_222) or\
             (origSym == cts.SYM_TETRAHEDRAL_Z3R and
              targetSym == cts.SYM_TETRAHEDRAL_222):
            matrix = self.coordinateSystemTransform(targetSym, origSym)
            return _invertMatrix(matrix)

        else:
            raise Exception("unknown Tetrahedral symmetry pair: %s %s" % (
                cts.SCIPION_SYM_NAME[origSym],
                cts.SCIPION_SYM_NAME[targetSym]))


class Octahedral(object):
    """Octahedral class.
        Allow to compute symmetry matrices,
        symmetry axis and unit cell planes of an octahedral symmetry.
        4-folds along x, y, z axes.
    """
    def __init__(self, circumscribed_radius=1, center=(
            0, 0, 0), offset=None):
        """
        :Parameters:
            circumscribed_radius: float radius of the circumscribed sphere
            center: tuple of 3 floats, center of the circumscribed sphere
            offset: float, angle in degrees to rotate the symmetry axis
                    it modifies the unit cell but not the symmetry matrices
            """
        self.circumscribed_radius = circumscribed_radius
        self.center = center
        self.matrices = None
        self.offset = offset

    def symmetryMatrices(self):
        """ get Matrices for cyclic symmetry of order n
        """
        self.c4 = (((0, 0, 1), 0),
                   ((0, 0, 1), 90),
                   ((0, 0, 1), 180),
                   ((0, 0, 1), 270))
        self.cube = (((1, 0, 0), 0),
                     ((1, 0, 0), 90),
                     ((1, 0, 0), 180),
                     ((1, 0, 0), 270),
                     ((0, 1, 0), 90),
                     ((0, 1, 0), 270))
        c4syms = [_rotationTransform(axis, angle)
                  for axis, angle in self.c4]
        cubesyms = [_rotationTransform(axis, angle)
                    for axis, angle in self.cube]
        syms = _matrixProducts(cubesyms, c4syms)
        self.matrices = _recenterSymmetries(syms, self.center)
        return self.matrices

    def unitCellPlanes(self):
        """ get planes that define a unit cell for an octahedral symmetry
        """
        _ = self.symmetryMatrices()
        # four corners tetahedron
        _4fold_1 = _normalizeVector(
            np.array([0., 0., 1.])) * self.circumscribed_radius
        _3fold_1 = _normalizeVector(
            np.array([1., 1., 1.])) * self.circumscribed_radius
        _3fold_2 = _normalizeVector(
            np.array([-1., 1., 1.])) * self.circumscribed_radius

        plane1 = _normalizeVector(
            np.cross(_3fold_2, _4fold_1))  # cyan
        plane2 = _normalizeVector(
            np.cross(_3fold_1, _3fold_2))  # green
        plane3 = _normalizeVector(
            np.cross(_4fold_1, _3fold_1))  # magenta
        if DEBUG:
            print("v0", _4fold_1)
            print("v1", _3fold_2)
            print("v2", _3fold_1)

            print("plane1", plane1)
            print("plane2", plane2)
            print("plane3", plane3)
        return [_4fold_1, _3fold_2, _3fold_1], [plane1, plane2, plane3]


# TODO: Why is this a global variable instead of
# a class variable?
icos_matrices = {}  # Maps orientation name to 60 matrices.


class Icosahedron(object):
    """Icosahedron class.
        Allow to compute symmetry matrices, icosahedron vertices,
         symmetry axis and unit cell planes
    """
    def __init__(self, circumscribed_radius=1, orientation='222', center=(
            0, 0, 0)):
        """point = np.array([0, 1, PHI]"""
        self.circumscribed_radius = circumscribed_radius
        self.orientation = orientation
        self.center = center
        # Triangle edge length of unit icosahedron.
        self.e = sqrt(2 - 2 / sqrt(5))
        self.vertices = None  # icosahedron vertices
        self.triangles = None  # icosahedron faces
        self._3foldAxis = []
        self._2foldAxis = []
        self.edges = None  # icosahedron edges
        # ---------------------------------------------------------------------------------------------
        # Coordinates systems.
        # '222'         2-fold symmetry along x, y, and z axes.
        # '222r'        '222' with 90 degree rotation around z.
        # '2n5'         2-fold symmetry along x and 5-fold along z.
        # '2n5r'        '2n5' with 180 degree rotation about y.
        # 'n25'         2-fold symmetry along y and 5-fold along z.
        # 'n25r'        'n25' with 180 degree rotation about x.
        # '2n3'         2-fold symmetry along x and 3-fold along z.
        # '2n3r'        '2n3' with 180 degree rotation about y.
        #
        self.coordinate_system_names = ('222', '222r', '2n5',
                                        '2n5r', 'n25', 'n25r', '2n3', '2n3r')
        self.icosahedronEdgeLength = \
            self.icosahedronEdgeLength(circumscribed_radius)
        self.angle23, self.angle25, self.angle35 = self.icosahedronAngles()
        self.icosahedronGeometry()

    # ---------------------------------------------------------------------------------------------
    # 60 icosahedral symmetry matrices.
    #
    def icosahedralSymmetryMatrices(self):
        t = self.icosahedralMatrixTable()
        tflist = _recenterSymmetries(t.get(self.orientation, None),
                                     self.center)
        return tflist

    # -------------------------------------------------------------------------
    # Edge length of icosahedron with a certain radio of the circumscribed
    # sphere:
    # According to Radio of the circumscribed sphere = 1 (vertices),
    #  the icosahedronEdgeLength should be 1.0515

    def icosahedronEdgeLength(self, circumscribed_radius):
        return 4 * circumscribed_radius / sqrt(10 + 2 * sqrt(5))

    # -------------------------------------------------------------------------
    # Matrices for mapping between different icosahedron coordinate frames.
    #
    def coordinateSystemTransform(self, from_cs, to_cs):

        self.cst = {}
        if self.cst:
            return self.cst[(from_cs, to_cs)]

        transform = self.cst

        s25 = self.e / 2  # Sin/Cos for angle between 2-fold and 5-fold axis
        c25 = sqrt(1 - s25 * s25)
        # Sin/Cos for angle between 3-fold and 5-fold axis
        s35 = self.e / sqrt(3)
        c35 = sqrt(1 - s35 * s35)

        transform[('2n5', '222')] = ((1, 0, 0, 0),
                                     (0, c25, -s25, 0),
                                     (0, s25, c25, 0))
        transform[('2n5', '2n3')] = ((1, 0, 0, 0),
                                     (0, c35, s35, 0),
                                     (0, -s35, c35, 0))

        # Axes permutations.
        # 90 degree rotation about z
        transform[('222', '222r')] = ((0, 1, 0, 0),
                                      (-1, 0, 0, 0),
                                      (0, 0, 1, 0))
        # 180 degree rotation about y
        transform[('2n3', '2n3r')] = \
            transform[('2n5', '2n5r')] = ((-1, 0, 0, 0),
                                          (0, 1, 0, 0),
                                          (0, 0, -1, 0))
        # 180 degree rotation about x
        transform[('n25', 'n25r')] = ((1, 0, 0, 0),
                                      (0, -1, 0, 0),
                                      (0, 0, -1, 0))
        # x <-> y and z -> -z
        transform[('n25', '2n5')] = ((0, 1, 0, 0),
                                     (1, 0, 0, 0),
                                     (0, 0, -1, 0))

        # Extend to all pairs of transforms.
        tlist = []
        while len(transform) > len(tlist):

            tlist = list(transform.keys())

            # Add inverse transforms
            for key in tlist:
                f = key[0]
                t = key[1]
                if not (t, f) in transform:
                    transform[(t, f)] = _transposeMatrix(transform[(f, t)])

            # Use transitivity
            for key1 in tlist:
                f1 = key1[0]
                t1 = key1[1]
                for key2 in tlist:
                    f2 = key2[0]
                    t2 = key2[1]
                    if f2 == t1 and f1 != t2 and not (f1, t2) in transform:
                        transform[(f1, t2)] = _multiplyMatrices(transform[(
                            f2, t2)], transform[(f1, t1)])

        i = _identityMatrix()
        for s in self.coordinate_system_names:
            transform[(s, s)] = i

        return transform[(from_cs, to_cs)]

    # ---------------------------------------------------------------------------------------
    # Compute icosahedral transformation matrices for different
    # coordinate systems.
    #

    def icosahedralMatrixTable(self):
        global icos_matrices
        if icos_matrices:
            return icos_matrices

        c = cos(2 * pi / 5)  # .309016994
        c2 = cos(4 * pi / 5)  # -.809016994

        icos_matrices['222'] = (

            ((1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0)),

            ((c2, -0.5, c, 0.0),
             (-0.5, c, c2, 0.0),
             (c, c2, -0.5, 0.0)),

            ((0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0)),

            ((-c2, -0.5, -c, 0.0),
             (-0.5, -c, c2, 0.0),
             (c, -c2, -0.5, 0.0)),

            ((0.5, c, c2, 0.0),
             (-c, c2, -0.5, 0.0),
             (c2, 0.5, -c, 0.0)),

            ((-c, c2, -0.5, 0.0),
             (c2, 0.5, -c, 0.0),
             (0.5, c, c2, 0.0)),

            ((c2, 0.5, -c, 0.0),
             (0.5, c, c2, 0.0),
             (-c, c2, -0.5, 0.0)),

            ((c2, -0.5, -c, 0.0),
             (0.5, -c, c2, 0.0),
             (c, c2, 0.5, 0.0)),

            ((-c, -c2, -0.5, 0.0),
             (c2, -0.5, -c, 0.0),
             (-0.5, c, -c2, 0.0)),

            ((0.5, -c, c2, 0.0),
             (-c, -c2, -0.5, 0.0),
             (-c2, 0.5, c, 0.0)),

            ((0.0, 0.0, -1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0)),

            ((-0.5, -c, c2, 0.0),
             (c, -c2, -0.5, 0.0),
             (-c2, -0.5, -c, 0.0)),

            ((-0.5, c, c2, 0.0),
             (c, c2, -0.5, 0.0),
             (c2, -0.5, c, 0.0)),

            ((-c, c2, -0.5, 0.0),
             (-c2, -0.5, c, 0.0),
             (-0.5, -c, -c2, 0.0)),

            ((c2, 0.5, -c, 0.0),
             (-0.5, -c, -c2, 0.0),
             (c, -c2, 0.5, 0.0)),

            ((0.5, c, c2, 0.0),
             (c, -c2, 0.5, 0.0),
             (-c2, -0.5, c, 0.0)),

            ((-0.5, c, c2, 0.0),
             (-c, -c2, 0.5, 0.0),
             (-c2, 0.5, -c, 0.0)),

            ((0.0, 0.0, -1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0)),

            ((-0.5, -c, c2, 0.0),
             (-c, c2, 0.5, 0.0),
             (c2, 0.5, c, 0.0)),

            ((0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0)),

            ((c2, 0.5, c, 0.0),
             (0.5, c, -c2, 0.0),
             (c, -c2, -0.5, 0.0)),

            ((-c2, 0.5, -c, 0.0),
             (0.5, -c, -c2, 0.0),
             (c, c2, -0.5, 0.0)),

            ((-c, -c2, -0.5, 0.0),
             (-c2, 0.5, c, 0.0),
             (0.5, -c, c2, 0.0)),

            ((0.5, -c, c2, 0.0),
             (c, c2, 0.5, 0.0),
             (c2, -0.5, -c, 0.0)),

            ((c2, -0.5, -c, 0.0),
             (-0.5, c, -c2, 0.0),
             (-c, -c2, -0.5, 0.0)),

            ((-c, c2, 0.5, 0.0),
             (c2, 0.5, c, 0.0),
             (-0.5, -c, c2, 0.0)),

            ((-c, -c2, 0.5, 0.0),
             (-c2, 0.5, -c, 0.0),
             (-0.5, c, c2, 0.0)),

            ((1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0)),

            ((c, -c2, -0.5, 0.0),
             (-c2, -0.5, -c, 0.0),
             (-0.5, -c, c2, 0.0)),

            ((c, c2, -0.5, 0.0),
             (c2, -0.5, c, 0.0),
             (-0.5, c, c2, 0.0)),

            ((-1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0)),

            ((-c2, 0.5, -c, 0.0),
             (-0.5, c, c2, 0.0),
             (-c, -c2, 0.5, 0.0)),

            ((0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, -1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0)),

            ((c2, 0.5, c, 0.0),
             (-0.5, -c, c2, 0.0),
             (-c, c2, 0.5, 0.0)),

            ((-0.5, -c, -c2, 0.0),
             (-c, c2, -0.5, 0.0),
             (-c2, -0.5, c, 0.0)),

            ((c, -c2, 0.5, 0.0),
             (c2, 0.5, -c, 0.0),
             (-0.5, -c, -c2, 0.0)),

            ((-c2, -0.5, c, 0.0),
             (0.5, c, c2, 0.0),
             (c, -c2, 0.5, 0.0)),

            ((-c2, 0.5, c, 0.0),
             (0.5, -c, c2, 0.0),
             (-c, -c2, -0.5, 0.0)),

            ((c, c2, 0.5, 0.0),
             (c2, -0.5, -c, 0.0),
             (0.5, -c, c2, 0.0)),

            ((-0.5, c, -c2, 0.0),
             (-c, -c2, -0.5, 0.0),
             (c2, -0.5, -c, 0.0)),

            ((0.0, 0.0, 1.0, 0.0),
             (-1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0)),

            ((0.5, c, -c2, 0.0),
             (c, -c2, -0.5, 0.0),
             (c2, 0.5, c, 0.0)),

            ((0.5, -c, -c2, 0.0),
             (c, c2, -0.5, 0.0),
             (-c2, 0.5, -c, 0.0)),

            ((c, -c2, 0.5, 0.0),
             (-c2, -0.5, c, 0.0),
             (0.5, c, c2, 0.0)),

            ((-c2, -0.5, c, 0.0),
             (-0.5, -c, -c2, 0.0),
             (-c, c2, -0.5, 0.0)),

            ((-0.5, -c, -c2, 0.0),
             (c, -c2, 0.5, 0.0),
             (c2, 0.5, -c, 0.0)),

            ((0.5, -c, -c2, 0.0),
             (-c, -c2, 0.5, 0.0),
             (c2, -0.5, c, 0.0)),

            ((0.0, 0.0, 1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0),
             (0.0, 1.0, 0.0, 0.0)),

            ((0.5, c, -c2, 0.0),
             (-c, c2, 0.5, 0.0),
             (-c2, -0.5, -c, 0.0)),

            ((0.0, 1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0),
             (1.0, 0.0, 0.0, 0.0)),

            ((-c2, -0.5, -c, 0.0),
             (0.5, c, -c2, 0.0),
             (-c, c2, 0.5, 0.0)),

            ((c2, -0.5, c, 0.0),
             (0.5, -c, -c2, 0.0),
             (-c, -c2, 0.5, 0.0)),

            ((c, c2, 0.5, 0.0),
             (-c2, 0.5, c, 0.0),
             (-0.5, c, -c2, 0.0)),

            ((-0.5, c, -c2, 0.0),
             (c, c2, 0.5, 0.0),
             (-c2, 0.5, c, 0.0)),

            ((-c2, 0.5, c, 0.0),
             (-0.5, c, -c2, 0.0),
             (c, c2, 0.5, 0.0)),

            ((c, -c2, -0.5, 0.0),
             (c2, 0.5, c, 0.0),
             (0.5, c, -c2, 0.0)),

            ((c, c2, -0.5, 0.0),
             (-c2, 0.5, -c, 0.0),
             (0.5, -c, -c2, 0.0)),

            ((-1.0, 0.0, 0.0, 0.0),
             (0.0, -1.0, 0.0, 0.0),
             (0.0, 0.0, 1.0, 0.0)),

            ((-c, c2, 0.5, 0.0),
             (-c2, -0.5, -c, 0.0),
             (0.5, c, -c2, 0.0)),

            ((-c, -c2, 0.5, 0.0),
             (c2, -0.5, c, 0.0),
             (0.5, -c, -c2, 0.0)),

        )

        for cs in self.coordinate_system_names:
            if cs != '222':
                t = self.coordinateSystemTransform(cs, '222')
                tinv = self.coordinateSystemTransform('222', cs)
                icos_matrices[cs] = [_multiplyMatrices(tinv, m, t)
                                     for m in icos_matrices['222']]
        return icos_matrices

    # -----------------------------------------------------------------------------

    def icosahedronAngles(self):
        # Angles between 2-fold, 3-fold and 5-fold
        # symmetry axes of an icosahedron.
        angle25 = asin(self.e/2)
        angle35 = asin(self.e/sqrt(3))
        angle23 = asin(self.e/(2*sqrt(3)))
        return angle23, angle25, angle35

    def icosahedronGeometry(self):

        a = 2 * self.angle25  # Angle spanned by edge from center
        # 5-fold symmetry axis along z
        c5 = cos(2 * pi / 5)
        s5 = sin(2 * pi / 5)
        tf5 = ((c5, -s5, 0, 0),
               (s5, c5, 0, 0),
               (0, 0, 1, 0))

        # 2-fold symmetry axis along x
        tf2 = ((1, 0, 0, 0),
               (0, -1, 0, 0),
               (0, 0, -1, 0))

        p = tuple(map(operator.add, self.center,
                      (0, 0, self.circumscribed_radius)))
        p50 = tuple(map(operator.add, self.center,
                        (0, self.circumscribed_radius * sin(a),
                         self.circumscribed_radius * cos(a))))
        p51 = _applyMatrix(tf5, p50)
        p52 = _applyMatrix(tf5, p51)
        p53 = _applyMatrix(tf5, p52)
        p54 = _applyMatrix(tf5, p53)
        self.vertices = [p, p50, p51, p52, p53, p54]
        self.vertices.extend([_applyMatrix(tf2, q) for q in self.vertices])
        if self.orientation != '2n5':
            self.tf = self.coordinateSystemTransform('2n5', self.orientation)
            self.vertices = [_applyMatrix(self.tf, p) for p in self.vertices]

            #
            # Vertex numbering
            #
            #  Top   1          Bottom
            #      2   5        9   10
            #        0            6
            #      3   4        8   11
            #                     7
        # 20 triangles composing icosahedron.
        #
        self.triangles = ((0, 1, 2), (0, 2, 3),
                          (0, 3, 4), (0, 4, 5),
                          (0, 5, 1), (6, 7, 8),
                          (6, 8, 9), (6, 9, 10),
                          (6, 10, 11), (6, 11, 7),
                          (1, 9, 2), (2, 9, 8),
                          (2, 8, 3), (3, 8, 7),
                          (3, 7, 4), (4, 7, 11),
                          (4, 11, 5), (5, 11, 10),
                          (5, 10, 1), (1, 10, 9))

        for triangle in self.triangles:
            self._3foldAxis.append([(item1 + item2 + item3) / 3.
                                    for item1,  item2,  item3 in zip(
                                        self.vertices[triangle[0]],
                                        self.vertices[triangle[1]],
                                        self.vertices[triangle[2]])])

        self.edges = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
                      (1, 2), (2, 3), (3, 4), (4, 5), (5, 1),
                      (6, 7), (6, 8), (6, 9), (6, 10), (6, 11),
                      (1, 9), (2, 8), (2, 9), (3, 7), (3, 8),
                      (7, 8), (8, 9), (9, 10), (10, 11), (11, 7),
                      (4, 7), (4, 11), (5, 10), (5, 11), (1, 10))

        for edge in self.edges:
            self._2foldAxis.append([(item1 + item2) / 2.
                                    for item1, item2 in zip(
                                        self.vertices[edge[0]],
                                        self.vertices[edge[1]])
                                    ])

    def getVertices(self):
        return self.vertices

    def getTriangles(self):
        return self.triangles

    def getEdges(self):
        return self.edges

    def get3foldAxis(self):
        return self._3foldAxis

    def get2foldAxis(self):
        return self._2foldAxis

    def unitCellPlanes(self):
        """ get planes that define a unit cell for an octahedral symmetry
        """
        if self.orientation == '222':
            # Unit cell in triangle 11 -> (2,8,9)
            _3f = 11; _5f1 = 2; _5f2 = 8  # noqa
        elif self.orientation == '222r':
            # Unit cell in triangle 0 -> (0,1,2)
            _3f = 0; _5f1 = 1; _5f2 = 0  # noqa
        elif self.orientation == 'n25':
            # Unit cell in triangle 9 -> (6,11,7)
            _3f = 9; _5f1 = 6; _5f2 = 7  # noqa
        elif self.orientation == 'n25r':
            # Unit cell in triangle 14 -> (3,7,4)
            _3f = 14; _5f1 = 3; _5f2 = 4  # noqa
        elif self.orientation == '2n3':  # TODO: Untested
            _3f = 11; _5f1 = 8; _5f2 = 2  # noqa
        elif self.orientation == '2n3r':  # TODO: Untested
            _3f = 11; _5f1 = 2; _5f2 = 8  # noqa
        elif self.orientation == '2n5':  # TODO: Untested
            _3f = 11; _5f1 = 8; _5f2 = 2  # noqa
        elif self.orientation == '2n5r':  # TODO: Untested
            _3f = 11; _5f1 = 2; _5f2 = 8  # noqa
        _3fold = self._3foldAxis[_3f]
        _5fold_1 = self.getVertices()[_5f1]
        _5fold_2 = self.getVertices()[_5f2]

        plane1 = _normalizeVector(np.cross(_5fold_1, _3fold))  # cyan
        plane2 = _normalizeVector(np.cross(_3fold, _5fold_2))  # magenta
        plane3 = _normalizeVector(np.cross(_5fold_2, _5fold_1))  # green
        if DEBUG:
            print("v0", _3fold)
            print("v1", _5fold_1)
            print("v2", _5fold_2)

            print("plane1", plane1)
            print("plane2", plane2)
            print("plane3", plane3)
        return [_3fold, _5fold_1, _5fold_2], [plane1, plane2, plane3]
