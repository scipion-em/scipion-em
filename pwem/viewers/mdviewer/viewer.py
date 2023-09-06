# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from pwem.objects import EMSet
from pyworkflow.viewer import Viewer, View
from metadataviewer.model import ObjectManager, ImageRenderer

from .readers import STKImageReader, MRCImageReader
from .sqlite_dao import SqliteFile
from ...protocols import ProtUserSubSet

# Register the readers only once
ImageRenderer.registerImageReader(MRCImageReader)
ImageRenderer.registerImageReader(STKImageReader)


class MDView(View):

    def __init__(self, emSet: EMSet, protocol, projectWindow):
        self._emSet = emSet
        self.protocol = protocol
        self.projectWindow = projectWindow

    def show(self):
        objectManager = ObjectManager()
        objectManager.registerDAO(SqliteFile)
        SqliteFile.userSubsetCreationCallback = self.userSubsetCreationCallback
        objectManager.open(self._emSet.getFileName())

    def userSubsetCreationCallback(self, subsetName, selectionFile, outputType):
        project = self.protocol.getProject()
        protocol = project.newProtocol(ProtUserSubSet)
        protocol.setObjLabel(subsetName)
        protocol.sqliteFile.set(selectionFile)
        protocol.inputObject.set(self._emSet)
        protocol.outputClassName.set(outputType)
        # protocol.other.set()
        self.projectWindow.getViewWidget().executeProtocol(protocol)


class MDViewer(Viewer):
    _name = 'Scipion metadata viewer'
    _targets = [EMSet]

    def _visualize(self, obj, **kwargs):
        return [MDView(obj, self.protocol, self.formWindow)]



