# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
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
from pyworkflow.gui import showInfo
from pwem.viewers import EmPlotter
from pyworkflow.viewer import Viewer
import pwem.objects as emobj

class AngularDistributionViewer(Viewer):
    """ Visualize the output of protocol reconstruct swarm """
    _label = 'Score transformation viewer'
    _targets = [emobj.SetOfParticles, emobj.SetOfVolumes]

    @staticmethod
    def plotAngularDistribution(emSet:emobj.EMSet, histogram=False):

        plotter = EmPlotter(x=1, y=1, windowTitle='Angular distribution')

        plotter.plotAngularDistributionFromSet(emSet, "Angular distribution", histogram=histogram)

        return plotter

    def _visualize(self, outputSet: emobj.EMSet, **kwargs):
        # Keep input object in case we need to launch
        # a new protocol and set dependencies

        views =[]

        # Weird plot. Not sure if used at all:
        # views.append(self.plotAngularDistribution(outputSet, histogram=True))

        firstItem = outputSet.getFirstItem()

        if not firstItem.hasTransform():
            showInfo("Missing alignment information", "This set does not have alignment information.")
            return

        views.append(self.plotAngularDistribution(outputSet))

        return views




