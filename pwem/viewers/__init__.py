# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# *
# **************************************************************************

# Expose many basic views
from pyworkflow.gui.browser import FileTreeProvider, STANDARD_IMAGE_EXTENSIONS
from .views import (DataView, ObjectView, MicrographsView, CtfView,
                    ClassesView, Classes3DView, CoordinatesObjectView,
                    ImageView, TableView)
from .viewer_base import EmProtocolViewer
from .plotter import EmPlotter
from .viewers_data import DataViewer

from .viewer_localres import LocalResolutionViewer
from .viewer_vmd import Vmd, VmdView, VmdViewer
from .viewer_fsc import FscViewer
from .viewer_pdf import PDFReportViewer
from .viewer_chimera import (Chimera, ChimeraView, ChimeraClientView,
                             ChimeraDataView, ChimeraViewer,
                             ChimeraClient)
from .viewer_sequence import SequenceViewer
from .viewer_volumes import viewerProtImportVolumes
from .showj import *

# register file handlers to preview info in the Filebrowser....
from .filehandlers import *
register = FileTreeProvider.registerFileHandler
register(MdFileHandler(), '.xmd', '.star', '.pos', '.ctfparam', '.doc')
register(ParticleFileHandler(),
         '.xmp', '.tif', '.tiff', '.spi', '.mrc', '.map', '.raw',
         '.inf', '.dm3', '.em', '.pif', '.psd', '.spe', '.ser', '.img',
         '.hed', *STANDARD_IMAGE_EXTENSIONS)
register(VolFileHandler(), '.vol')
register(StackHandler(), '.stk', '.mrcs', '.st', '.pif', '.dm4')
register(ChimeraHandler(), '.bild')
