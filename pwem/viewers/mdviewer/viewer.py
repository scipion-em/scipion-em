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
import os.path
import subprocess

from pwem.objects import EMSet
from pwem.viewers.viewers_data import RegistryViewerConfig
from pyworkflow import PYTHON
from pyworkflow.viewer import Viewer, View
from pwem.viewers.mdviewer.readers import SCIPION_PORT, SCIPION_OBJECT_ID


class MDView(View):

    def __init__(self, emSet: EMSet, protocol, port):
        self._emSet = emSet
        self.protocol = protocol
        self.port = port

    def show(self):
        env = os.environ
        env[SCIPION_PORT] = str(self.port)
        env[SCIPION_OBJECT_ID] = str(self._emSet.getObjId())
        visibleLabels = self.getVisibleLabels()

        subprocess.Popen(
            [PYTHON, "-m", "metadataviewer", "--extensionpath", os.path.join(os.path.dirname(__file__), "readers.py"),
            self._emSet.getFileName(), "--visiblelabels", visibleLabels])

    def getVisibleLabels(self):
        from pwem.viewers import VISIBLE
        config = RegistryViewerConfig.getConfig(type(self._emSet))
        if config is not None and VISIBLE in config:
            return config[VISIBLE]
        return ''


class MDViewer(Viewer):
    _name = 'Scipion metadata viewer'
    _targets = [EMSet]

    def _visualize(self, obj, **kwargs):
        return [MDView(obj, self.protocol, self._project.port)]




