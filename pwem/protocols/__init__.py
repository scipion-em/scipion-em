# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
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
In this module are protocol base classes related to EM.
There should be sub-classes in the different packages from
each EM-software package.
"""
from .protocol import *
from .protocol_micrographs import *
from .protocol_movies import *
from .protocol_align_movies import (ProtAlignMovies, ProtAverageFrames,
                                    createAlignmentPlot)
from .protocol_particles import *
from .protocol_particles_picking import *
from .protocol_2d import *
from .protocol_3d import *
from .protocol_sets import *
from .protocol_pdf_report import *
from .protocol_tiltpairs import *
from .protocol_ctf_assign import ProtCTFAssign
from .protocol_alignment_assign import ProtAlignmentAssign
from .protocol_alignment_invertHand import ProtAlignmentInvertHand
from .protocol_batch import *
from .protocol_extract_coordinates import ProtExtractCoords
from .protocol_create_stream_data import ProtCreateStreamData
from .parallel import ProtTestParallel
from .protocol_create_stream_data import SET_OF_RANDOM_MICROGRAPHS

from .protocol_import import *

from .protocol_export import *
from .protocol_tests import *

from .protocol_origin_sampling_volume import ProtOrigSampling
from .protocol_set_editor import ProtSetEditor
from .protocol_set_filter import ProtSetFilter
from .protocol_movie_eraser import ProtMovieEraser
from .protocol_classes_selector import ProtClassesSelector
