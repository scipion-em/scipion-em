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
from pwem.viewers.mdviewer.sqlite_dao import SCIPION_PORT, SCIPION_OBJECT_ID


class MDView(View):

    def __init__(self, emSet: EMSet, protocol=None, port=None):
        self._emSet = emSet
        self.protocol = protocol
        self.port = port

    def show(self):
        env = os.environ
        if self.port:
            env[SCIPION_PORT] = str(self.port)
            env[SCIPION_OBJECT_ID] = str(self._emSet.getObjId())
            visibleLabels, orderLabels = self.getVisibleAndOrderLabels()
            orderLabels = orderLabels.replace(' _filename', ' stack', 1)
            fn = self._emSet.getFileName()
        else:
            visibleLabels = orderLabels = ""
            fn = self._emSet

        subprocess.Popen(
            [PYTHON, "-m", "metadataviewer", "--extensionpath", os.path.join(os.path.dirname(__file__), "readers.py"),
            fn, "--visiblelabels", visibleLabels, "--orderlabels", orderLabels])

    def getVisibleAndOrderLabels(self):
        from pwem.viewers import VISIBLE, ORDER
        config = RegistryViewerConfig.getConfig(type(self._emSet))
        visible = config[VISIBLE] if config and VISIBLE in config else ''
        order = config[ORDER] if config and ORDER in config else ''
        return visible, order


class MDViewer(Viewer):
    _name = 'Scipion'
    _targets = [EMSet]

    def _visualize(self, obj, **kwargs):
        return [MDView(obj, self.protocol, self._project.port)]




