# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This modules contains constants related to EM
"""
# ------------------ Constants values --------------------------------------

NO_INDEX = 0  # This is the index value for single images
URL = 'https://github.com/scipion-em/scipion-em'
# Sampling rate input modes
SAMPLING_FROM_IMAGE = 0
SAMPLING_FROM_SCANNER = 1

RELATION_CTF = 'relation_ctf'

UNIT_PIXEL = 'px'
UNIT_PIXEL_FOURIER = '1/px'
UNIT_ANGSTROM = 'A'
UNIT_ANGSTROM_FOURIER = '1/A'
# No needed: Å should be valid as long as first line
# in file has utf-8 directive --> # -*- coding: utf-8 -*-
UNIT_ANGSTROM_SYMBOL = u"\u212B"

# Fourier Filter options
FILTER_LOW_PASS = 0
FILTER_HIGH_PASS = 1
FILTER_BAND_PASS = 2
FILTER_GAUSSIAN = 3
FILTER_LOW_PASS_NO_DECAY = 4
FILTER_NO_DECAY = 5

# Transform
ALIGN_NONE = 'None'
ALIGN_2D = '2D'            # 2D image alignment
ALIGN_3D = '3D'            # 3D map alignment
ALIGN_PROJ = 'Projection'  # relate projections with 3d map

ALIGNMENTS = [ALIGN_NONE, ALIGN_2D, ALIGN_3D, ALIGN_PROJ]

# Constants related with colormaps for viewers
# Color maps
COLOR_JET = 0
COLOR_TERRAIN = 1
COLOR_GIST_EARTH = 2
COLOR_GIST_NCAR = 3
COLOR_GNU_PLOT = 4
COLOR_GNU_PLOT2 = 5
COLOR_OTHER = 6

COLOR_CHOICES = ['jet', 'terrain', 'gist_earth',
                 'gist_ncar', 'gnuplot', 'gnuplot2', 'other']

# Axis code
AX_X = 0
AX_Y = 1
AX_Z = 2

# SYMMETRY, conventions described at:
# https://scipion-em.github.io/docs/docs/developer/symmetries
# Some notes:
# Icosahedral
#   xmipp and relion define I1, I2, I3 and I4 that corerspond to
#   I222, I222r, In25 and In25z
#   cryosparc has I1 and I2
#   EMAN always puts the highest symmetry axis on Z.
#   and a 2 fold axis in X (so it uses I2n5 or I2n5r, not sure)
# DIEDRAL: first axis Z, second X (DX) except cryosparc that uses y (DY)
# Tetraedral, most of the programs use TZ3
# CN cyclic symmetry Cn around z axis
SYM_CYCLIC = 0
# DN cyclic symmetry plus and extra 2fold symmetry axis around X
SYM_DIHEDRAL = 1
SYM_DIHEDRAL_X = SYM_DIHEDRAL
# DN cyclic symmetry plus and extra 2fold symmetry axis around Y (cryoSparc)
SYM_DIHEDRAL_Y = 2
# T 222 Tetrahedral symmetry with two-fold symmetry axis along the X, Y and Z
SYM_TETRAHEDRAL = 3
SYM_TETRAHEDRAL_222 = SYM_TETRAHEDRAL
# T  three-fold symmetry axis along Z, another three-fold axis
# in the YZ plane such that the sign(y) != sign(z)
SYM_TETRAHEDRAL_Z3 = 4
# SYM_TETRAHEDRAL222R  is a degenerated symmetry, that is
# for each symmetry axis (x,y,z)  defined in T222, there is a
# symmetry axis (-x,-y,-z) in T222R, therefore the produced
# symmetry matrices are equivalent
# SYM_TETRAHEDRAL222R = 14

# T  three-fold symmetry axis along Z, another three-fold axis
# in the YZ plane such that  sign(y) = sign(z)
SYM_TETRAHEDRAL_Z3R = 5

SYM_OCTAHEDRAL = 6  # O

# icosahedric IXXX
# (no crowther 222 and standard in heyman et al 2005 article).
# 2-fold axes on x,y,z axes. With the positive z-axis pointing at the viewer,
# the front-most 5-fold vertices are in yz plane, and the front-most 3-fold
# axes are in the xz plane.
SYM_I222 = 7

# (crowther) 2-fold axes on x,y,z axes. With the positive z-axis pointing at
# the viewer, the front-most 5-fold vertices are in xz plane,
# and the front-most 3-fold axes are in the yz plane.
# (222 rotated 90 degrees around Z)
SYM_I222r = 8

# '2-fold symmetry along y and 5-fold along z
SYM_In25 = 9

# 'n25' with 180 degree rotation about x
SYM_In25r = 10
SYM_I2n3 = 11  # Two-fold symmetry along X and 3-fold along Z
SYM_I2n3r = 12  # Idem but rotated 180 degree about Y
SYM_I2n5 = 13  # Two-fold symmetry along Y and 5-fold along Z
SYM_I2n5r = 14  # Idem but rotated 180 degree about X

# Symmetry dictionary
SCIPION_SYM_NAME = dict()
SCIPION_SYM_NAME[SYM_CYCLIC] = 'Cn'
SCIPION_SYM_NAME[SYM_DIHEDRAL_X] = 'Dxn'
SCIPION_SYM_NAME[SYM_DIHEDRAL_Y] = 'Dyn'
SCIPION_SYM_NAME[SYM_TETRAHEDRAL] = 'T222'
SCIPION_SYM_NAME[SYM_TETRAHEDRAL_Z3] = 'Tz3'
SCIPION_SYM_NAME[SYM_TETRAHEDRAL_Z3R] = 'Tz3r'
SCIPION_SYM_NAME[SYM_OCTAHEDRAL] = 'O'
SCIPION_SYM_NAME[SYM_I222] = 'I222'
SCIPION_SYM_NAME[SYM_I222r] = 'I222r'
SCIPION_SYM_NAME[SYM_In25] = 'In25'
SCIPION_SYM_NAME[SYM_In25r] = 'In25r'
SCIPION_SYM_NAME[SYM_I2n3] = 'I2n3'
SCIPION_SYM_NAME[SYM_I2n3r] = 'I2n3r'
SCIPION_SYM_NAME[SYM_I2n5] = 'I2n5'
SCIPION_SYM_NAME[SYM_I2n5r] = 'I2n5r'

# Entry points names
CONVERT_ENTRY_POINT = 'emconvert'
CHIMERA_ENTRY_POINT = 'emchimera'
EM_PROGRAM_ENTRY_POINT = 'emprogram'

# maxit
MAXIT_HOME = 'MAXIT_HOME'
MAXIT = 'maxit'

# Explanation messages for the wizards
RING_MASK_WIZ_MSG = \
    'The values of the outer and inner radius can be controlled via '\
    'both the sliders placed below or '\
    'the mousewheel (when the mouse cursor is over the image). '\
    'The arrow (=>) located at the left of '\
    'the sliders is used to indicate which one is active. '\
    'To change from one to another, simply '\
    'click and drag the corresponding slider or use up and down keys.'
CIRCLE_MASK_WIZ_MSG = \
    'The values of the mask radius can be controlled via '\
    'both the slider or the mousewheel (when ' \
    'the mouse cursor is over the image).'
FREQ_BANDPASS_WIZ_MSG =\
    'The values of the low and high frequencies can '\
    'be controlled via both the sliders placed below or '\
    'the mousewheel (when the mouse cursor is over the image). '\
    'The arrow (=>) located at the left of '\
    'the sliders is used to indicate which one is active. '\
    'To change from one to another, simply '\
    'click and drag the corresponding slider or use up and down keys. '\
    'Values are shown in both digital '\
    'frequency and angstroms.'

EM_ROOT_VAR = 'EM_ROOT'

DEFAULT_MAX_PREVIEW_FILE_SIZE = 500  # Unit is MB: 500 MB

RESIDUES3TO1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

RESIDUES1TO3 = {v: k for k, v in RESIDUES3TO1.items()}
