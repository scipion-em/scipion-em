# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

# Expose many basic views
from pyworkflow.utils import weakImport
from .views import (DataView, ObjectView, MicrographsView, CtfView,
                    ClassesView, Classes3DView, CoordinatesObjectView,
                    ImageView, TableView)
from .viewer_base import EmProtocolViewer
from .plotter import EmPlotter
from .viewers_data import DataViewer

from .viewer_localres import LocalResolutionViewer, ChimeraAttributeViewer
from .viewer_vmd import Vmd, VmdView, VmdViewer
from .viewer_fsc import FscViewer
from .viewer_pdf import PDFReportViewer
from .viewer_chimera import (Chimera, ChimeraView, ChimeraViewer, ChimeraAngDist)
from .viewer_sequence import SequenceViewer
from .viewer_volumes import viewerProtImportVolumes
from .viewer_angular_dist import AngularDistributionViewer
from .viewers_data import BasicMDViewer

with weakImport("metadataviewer"):
    from .mdviewer import MDViewer

from .showj import *

if emConfig.CHIMERA_OLD_BINARY_PATH:
    from .viewer_chimera import ChimeraOldViewer
