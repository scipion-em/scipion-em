# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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

from matplotlib.ticker import FuncFormatter
from matplotlib.widgets import Button

from pyworkflow import Config
import pyworkflow.viewer as pwviewer

import pwem.objects as emobj
import pwem.protocols as emprot

from .plotter import EmPlotter, plt


class FscViewer(pwviewer.Viewer):
    """ Viewer for plot a FSC object. """
    _targets = [emobj.FSC, emobj.SetOfFSCs]
    _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self.threshold = kwargs.get('threshold', 0.143)
        self.legend = ""
        self.plotter = EmPlotter(x=1, y=1, windowTitle='FSC',
                                 figure=kwargs.get('figure', 'active'))
        self._lastSubPlot = None
        self._addButton = kwargs.get('addButton', False)

    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1 / value
        return "1/%0.2f" % inv

    def show(self, block=True):
        self.plotter.show(block=block)

    def _plotButton(self):
        if self._addButton:
            axcreateFSC = plt.axes([0.75, 0.02, 0.2, 0.050])
            # Button does not allow to define text color so
            # I write it directly
            axcreateFSC.text(0.5, 0.5, 'create FSC',
                             verticalalignment='center',
                             horizontalalignment='center',
                             transform=axcreateFSC.transAxes, color='white')
            bcreateFSC = Button(axcreateFSC, '',  # leave label empty
                                color=Config.SCIPION_MAIN_COLOR,
                                hovercolor='maroon')
            bcreateFSC.on_clicked(self.createFSCObject)
            self._addButton = False
            return bcreateFSC

    def visualize(self, obj, **kwargs):
        # Keep input object in case we need to launch
        # a new protocol and set dependencies
        self.fscList = []
        if isinstance(obj, emobj.SetOfFSCs):
            for i, fsc in enumerate(obj):
                self.fscList.append(fsc.clone())
                fscLabel = fsc.getObjLabel() or 'FSC %d' % (i + 1)
                self.plotFsc(fsc, label=fscLabel, **kwargs)
        else:  # single FSC
            fscLabel = obj.getObjLabel() or 'FSC 1'
            self.plotFsc(obj, label=fscLabel, **kwargs)
            self.fscList.append(obj.clone())

        bcreateFSC = self._plotButton()  # if you do not assign the result
        # to something it will be
        # delete
        self.show()
        # return [self.plotter]

    def createFSCObject(self, event):
        self.project = self.getProject()
        prot = self.getProject().newProtocol(emprot.ProtCreateFSC)
        prot.setObjLabel('FSC-%s' % self.label)  # label)
        prot.setInputObj(self.protocol)
        prot.setInputFscList(self.fscList)
        self.project.launchProtocol(prot)

    def plotFsc(self, obj, **kwargs):

        x, y = obj.getData()

        if self._lastSubPlot is None:
            a = self.plotter.createSubPlot("FSC", "frequency (1/A)", "fsc")
            a.set_ylim([-0.1, 1.1])
            a.plot([0, x[-1]],
                   [self.threshold, self.threshold],
                   'k--')
            a.grid(True)
            a.xaxis.set_major_formatter(FuncFormatter(self._formatFreq))
            self._lastSubPlot = a
        resolution = obj.calculateResolution(self.threshold)
        self.label = kwargs.get('label', self.protocol.getRunName()) + " (" + str(resolution) + "Å)"
        self.plotter.plotData(x, y, '-', label=self.label)
        self.plotter.legend()
        # if I run this
