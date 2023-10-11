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
from matplotlib import cm

from pwem.viewers import EmPlotter
from pyworkflow.viewer import ProtocolViewer
import pwem.objects as emobj
from  pyworkflow.protocol.params import LabelParam

class AngularDistributionViewer(ProtocolViewer):
    """ Visualize particles, subtomograms with orientation information """
    _label = "Angular distribution viewer"
    _targets = [emobj.SetOfParticles, emobj.SetOfVolumes]


    def _defineParams(self, form):

        form.addSection(label='Visualization')

        form.addParam(self.getParamName(self.doShowHeatMap), LabelParam,
                      label="Show heatmap")

        form.addParam(self.getParamName(self.doShow2DPolar), LabelParam,
                      label="Show 2D polar plot")

        form.addParam(self.getParamName(self.doShow3DPlot), LabelParam,
                      label="Show 3D plot")

        from pwem.wizards import ColorScaleWizardBase
        ColorScaleWizardBase.defineColorScaleParams(form, defaultHighest=0,
                                                    defaultLowest=1)

    def getParamName(self, method):
        return method.__name__ + "P"

    def _getVisualizeDict(self):

        return {self.getParamName(self.doShowHeatMap): self.doShowHeatMap,
                self.getParamName(self.doShow2DPolar): self.doShow2DPolar,
                self.getParamName(self.doShow3DPlot): self.doShow3DPlot,
                }

    def getColorMap(self):
        cmap = cm.get_cmap(self.colorMap.get())
        if cmap is None:
            cmap = cm.jet
        return cmap

    @staticmethod
    def plotAngularDistribution(emSet:emobj.EMSet, colormap, type=1):

        plotter = EmPlotter(x=1, y=1, windowTitle='Angular distribution')

        plotter.plotAngularDistributionFromSet(emSet, "Angular distribution", type=type,
                                               colormap=colormap)

        return [plotter]

    def doShowHeatMap(self, e):
        return self.plotAngularDistribution(self.protocol, self.getColorMap(), type=1)

    def doShow2DPolar(self,e):
        return self.plotAngularDistribution(self.protocol, self.getColorMap(), type=2)

    def doShow3DPlot(self,e):
        return self.plotAngularDistribution(self.protocol, self.getColorMap(), type=3)
