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

import pyworkflow.viewer as pwviewer

import pwem.protocols as emprot

from .views import DataView

from pyworkflow.protocol.params import LabelParam
import matplotlib.pyplot as plt
import os
from .views import ObjectView
from pwem.viewers import EmProtocolViewer

# class ViewerClassesConsensus(pwviewer.Viewer):
#     _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]
#     _targets = [emprot.ProtClassesConsensus]
#
#     def _visualize(self, obj, **kwargs):
#         labels = ('class1.id class1._representative._filename class2.id '
#                   'class2._representative._filename jaccard intersection union')
#         return [DataView(obj.outputConsensus.getFileName(),
#                          viewParams={'order': labels, 'mode': 'metadata',
#                                      'visible': labels,
#                                      'render': 'class1._representative._filename class2._representative._filename'
#                                      })
#                 ]
#
#     def visualize(self, obj, **kwargs):
#         self._visualize(obj, **kwargs)[0].show()


class ViewerGoodClassesExtractor(EmProtocolViewer):
    _label = 'viewer Good Classes Extractor'
    _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]
    _targets = [emprot.ProtGoodClassesExtractor]

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('visualizeParticles', LabelParam,
                      label="Visualize accepted particles",
                      help="Visualize particles that come from good 2D classes.")
        form.addParam('visualizeDiscardedParticles', LabelParam,
                      label="Visualize discarded particles",
                      help="Visualize particles that come from bad 2D classes.")
        form.addParam('visualizeDistribution', LabelParam,
                      label="Visualize particles distribution",
                      help="Visualize plot that shows the particles distribution in good and bad classes.")
        form.addParam('visualizeDistributionVsTime', LabelParam,
                      label="Visualize particles distribution vs Time",
                      help="Visualize plot particles distribution vs time.")

    def _getVisualizeDict(self):
        return {
                 'visualizeParticles': self._visualizeParticlesF,
                 'visualizeDiscardedParticles': self._visualizeDiscardedParticlesF,
                 'visualizeDistribution': self._visualizeDistribution,
                 'visualizeDistributionVsTime': self._visualizeDistributionTime
                }

    def _visualizeParticlesF(self, e=None):
        return self._visualizeParticles("outputParticles")

    def _visualizeDiscardedParticlesF(self, e=None):
        return self._visualizeParticles("outputParticlesDiscarded")

    def _visualizeDistribution(self, e=None):
        if os.path.exists(self.protocol.getDistributionPlot()):
            image = plt.imread(self.protocol.getDistributionPlot())
            plt.figure()
            fig = plt.imshow(image)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.show()

    def _visualizeDistributionTime(self, e=None):
        if os.path.exists(self.protocol.getDistributionTimePlot()):
            image = plt.imread(self.protocol.getDistributionTimePlot())
            plt.figure()
            fig = plt.imshow(image)
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            plt.show()

    def _visualizeParticles(self, objName):
        views = []

        labels = 'id _filename _samplingRate _acquisition._dosePerFrame ' \
                 '_acquisition._doseInitial _MEAN_DOSE_PER_ANGSTROM2 _STD_DOSE_PER_ANGSTROM2 ' \
                 '_DIFF_TO_DOSE_PER_ANGSTROM2 '

        if self.protocol.hasAttribute(objName):
            setParticles = getattr(self.protocol, objName)
            views.append(ObjectView(
                self._project, setParticles.getObjId(), setParticles.getFileName()))
                #viewParams={'mode': MODE_MD, 'order': labels, 'visible': labels}))
        else:
            self.infoMessage('%s does not have %s%s'
                             % (self.protocol.getObjLabel(), objName,
                                getStringIfActive(self.protocol)),
                             title='Info message').show()
        return views

def getStringIfActive(prot):
    return ', yet.' if prot.isActive() else '.'