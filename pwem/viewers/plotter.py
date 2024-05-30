# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
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
import logging

from pwem.convert.trigonometry import FibonacciSphere

logger = logging.getLogger(__name__)
from math import radians, degrees
import numpy as np
import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter

from pwem.convert.transformations import euler_from_matrix
from pyworkflow.gui.plotter import Plotter, plt
import pwem.emlib.metadata as md
import numbers
from math import atan2, sqrt, pi
from pwem import emlib

PLOT_EULER_ANGLES = 1
PLOT_PROJ_ANGLES = 2
PLOT_PROJ_DIR = 3

class EmPlotter(Plotter):
    """ Class to create several plots. """

    def __init__(self, x=1, y=1, mainTitle="", **kwargs):
        Plotter.__init__(self, x, y, mainTitle, **kwargs)

    def plotAngularDistribution(self, title, rot,
                                tilt, weight=[], max_p=40,
                                min_p=5, color='blue', colormap=None, subtitle=None):
        """ Create a special type of subplot, representing the angular
        distribution in 2d of weighted projections. """
        if weight:
            max_w = max(weight)
            min_w = 0
            if subtitle is None:
                subtitle = 'Min weight=%(min_w).2f, Max weight=%(max_w).2f' % locals()
            a = self.createSubPlot(title,
                                   subtitle,
                                    '', projection='polar')

            label_position = a.get_rlabel_position()
            a.text(np.radians(label_position + 10), a.get_rmax() / 2., 'Tilt',
                    rotation=label_position, ha='center', va='center')

            pointSizes = []
            for r, t, w in zip(rot, tilt, weight):
                 pointsize = int((w - min_w) / (max_w - min_w + 0.001) * (max_p - min_p) + min_p)
                 pointSizes.append(pointsize)

            if colormap:
                sc = a.scatter(rot, tilt, s=20, c=pointSizes, cmap=colormap, marker='.')
                plt.colorbar(sc)
            else:
                a.scatter(rot, tilt, s=pointSizes, c=color, marker='.')
        else:
            a = self.createSubPlot(title, 'Non weighted plot', '', projection='polar')
            a.scatter(rot, tilt, s=10, c=color, marker='.')
        return a


    def scatter3DPlot(self, x, y,z, title="3D scatter plot", drawsphere=True, markerSize=1, colormap=cm.jet, subtitle=None):

        ax =self.createSubPlot(title, xlabel=None, ylabel=None, projection='3d', subtitle=subtitle)
        ax.set_box_aspect(aspect=(1, 1, 1))
        sc = ax.scatter(x, y, z, markerSize,c=markerSize, s=60, cmap=colormap, alpha=1)
        plt.colorbar(sc)


        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        # ax.set_xlim([-1,1])
        # ax.set_ylim([-1, 1])
        # ax.set_zlim([-1, 1])

        # Add a color bar which maps values to colors.
        # self.figure.colorbar(cm.ScalarMappable(norm=norm, cmap=cmhot), shrink=0.5, aspect=5)

        if drawsphere:
            # draw sphere
            u, v = np.mgrid[0:2 * np.pi:21j, 0:np.pi:11j]
            x1 = np.cos(u) * np.sin(v)
            y1 = np.sin(u) * np.sin(v)
            z1 = np.cos(v)
            ax.plot_wireframe(x1, y1, z1, color="black", alpha=0.1)
            ax.plot_surface(x1, y1, z1, color="red", alpha=0.05)

            # 3d axis
            ax.plot([0, 1.25], [0, 0], [0, 0], color="red")
            ax.plot([0, 0], [0, 1.25], [0, 0], color="green")
            ax.plot([0, 0], [0, 0], [0, 1.25], color="blue")

        return ax

    def plotAngularDistribution3D(self, title, x,y,z, weights, subtitle, colormap=cm.jet):

        return self.scatter3DPlot(x,y,z, title=title,markerSize=weights, colormap=colormap, subtitle=subtitle)


    def plotAngularDistributionHistogram(self, title, data, eulerAnglesGetterCallback, colormap=cm.jet, subtitle=None):
        """ Create a special type of subplot, representing the angular
        distribution of projections. """

        # Extract the rot and tilt from the data
        thetas = []
        phis = []
        for item in data:
            rot, tilt, psi = eulerAnglesGetterCallback(item)
            thetas.append(rot)
            phis.append(tilt)

        thetas.append(-180)
        thetas.append(180)
        phis.append(0)
        phis.append(180)

        heatmap, xedges, yedges = np.histogram2d(thetas, phis, bins=1000)
        sigma = min(max(xedges) - min(xedges), max(yedges) - min(yedges)) / 20
        heatmap = gaussian_filter(heatmap, sigma=sigma)
        heatmapImage = heatmap.T
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        a = self.createSubPlot(title, 'Angular distribution', '', subtitle=subtitle)
        mappable = a.imshow(heatmapImage, extent=extent, origin='lower', cmap=colormap, aspect='auto')
        a.set_xlabel('Rotational angle')
        a.set_ylabel('Tilt angle')
        plt.colorbar(mappable)
        return mappable

    def plotAngularDistributionFromSet(self, mdSet, title, type=PLOT_PROJ_ANGLES, colormap=cm.jet, subtitle="", **kwargs):
        """ Read the values of the transformation matrix
                 and plot its histogram or its angular distribution.

                 :param type: Type of plot 1=histogram, 2=2D polar, 3=3D
                 :param mdSet: Set with alignment information at with item.getTransform()
                 :param title: Title of the plot
                 :param colormap: matplotlib color map
                """

        def eulerAnglesGetter (item):

            matrix = item.getTransform().getRotationMatrix()
            matrixI = np.linalg.inv(matrix)
            rot, tilt, psi = euler_from_matrix(matrix=matrixI, axes='szyz')

            if tilt<0:
                tilt= -tilt
                rot = -rot

            return degrees(rot), degrees(tilt), degrees(psi)


        self.plotAngularDistributionBase(mdSet, eulerAnglesGetter, title, type, colormap, subtitle=subtitle,**kwargs)


    def plotAngularDistributionBase(self, data, eulerAnglesGetterCallback, title,
                                    type=PLOT_PROJ_ANGLES, colormap=cm.jet, subtitle="", **kwargs):
        """ Read the values of the transformation matrix
         and plot its histogram or its angular distribution.

         :param colormap: matplotlib color map
         :param data: any particles iterator containing particles with alignment information. sqlites or star files , ...
         :param eulerAnglesGetterCallback: a callback to extract rot, tilt, psi IN DEGREES from each row, receiving the row/item.
         :param title: Title of the plot.
         :param type: 1 for histogram, 2 for polar plot, 3 for 3d plot
         :param subtitle: subtitle of the plot
        """

        if type==PLOT_EULER_ANGLES:
            return self.plotAngularDistributionHistogram(title, data , eulerAnglesGetterCallback,
                                                         colormap=colormap, subtitle=subtitle)

        else:

            rots =[]
            tilts= []

            # Get the euler angles
            for item in data:
                rot, tilt, psi = eulerAnglesGetterCallback(item)
                rots.append(rot)
                tilts.append(tilt)

            if type==PLOT_PROJ_ANGLES:
                # Weight (group) rots and tilts
                rots, tilts, weights = self.weightEulerAngles(rots, tilts, rotInRadians=type == PLOT_PROJ_ANGLES)

                return self.plotAngularDistribution(title, rots, tilts, weight=weights, colormap=colormap, subtitle=subtitle, **kwargs)

            else:
                # Create discrete point of a fibonacci sphere (5000 points)
                fiSph = FibonacciSphere()

                Xs, Ys, Zs = self._anglesToSphereCoords(rots, tilts)

                for  x,y,z in zip(Xs,Ys,Zs):
                    fiSph.append(x,y,z)

                fiSph.cleanWeights()

                return self.plotAngularDistribution3D(title, fiSph.sphX, fiSph.sphY, fiSph.sphZ,
                                                      fiSph.weights, subtitle,colormap=colormap)

    def weightEulerAngles(self, rots, tilts, delta=3, rotInRadians=False):
        """ Receives the list of rots and tilts angles (in deg) and returns
         a reduced list of rots, tilts and weights lists

         :param rotCaster: method that receives the rot in degrees and returns it converted to something else, radians?"""

        # Holds pairs of rot, tilt
        projectionList = []

        def getCloseProjectionIndex(angleRot, angleTilt):
            """ Get an existing projection close to angleRot, angleTilt.
            Return None if not found close enough.
            """
            for index, projection in enumerate(projectionList):
                if (abs(projection[0] - angleRot) <= delta and
                        abs(projection[1] - angleTilt) <= delta):
                    return index
            return None

        weight = 1 #1. / len(rots)

        new_rots = []
        new_tilts = []
        weights = []

        if rotInRadians:
            rotCaster = np.radians
        else:
            rotCaster = lambda rot: rot

        # Weight the rots and tilts
        for rot, tilt in zip(rots, tilts):
            projectionIndex = getCloseProjectionIndex(rot, tilt)

            if projectionIndex is None:
                projectionList.append([rot, tilt])
                new_rots.append(rotCaster(rot))
                new_tilts.append(tilt)
                weights.append(weight)
            else:
                weights[projectionIndex] += weight

        return new_rots, new_tilts, weights

    def _anglesToSphereCoords(self, rots, tilts):
        """ Converts euler angles (rot and tilts) to spherical coordinates."""

        X=[]
        Y=[]
        Z=[]

        for rot, tilt in zip(rots, tilts):

            # Converts to euler direction
            x, y, z = emlib.Euler_direction(rot, tilt, 0)

            X.append(x)
            Y.append(y)
            Z.append(z)

        return X, Y, Z

    def plotAngularDistributionFromMd(self, mdFile, title, **kwargs):
        """ Read the values of rot, tilt and weights from
        the metadata and plot the angular distribution.
        ANGLES are in DEGREES
        In the metadata:
            rot: MDL_ANGLE_ROT
            tilt: MDL_ANGLE_TILT
            weight: MDL_WEIGHT
        """
        angMd = md.MetaData(mdFile)

        if 'histogram' in kwargs:
            def eulerAnglesGetter(row):
                return row.getValue(md.MDL_ANGLE_ROT), row.getValue(md.MDL_ANGLE_TILT), None

            class MDIter:
                def __init__(self, mdObj):
                    self.mdObj = mdObj
                def __iter__(self):
                    for row in md.iterRows(self.mdObj):
                        yield row

            return self.plotAngularDistributionHistogram(title, MDIter(angMd), eulerAnglesGetter)
        else:
            rot = []
            tilt = []
            weight = []
            for row in md.iterRows(angMd):
                rot.append(radians(row.getValue(md.MDL_ANGLE_ROT)))
                tilt.append(row.getValue(md.MDL_ANGLE_TILT))
                weight.append(row.getValue(md.MDL_WEIGHT))
            return self.plotAngularDistribution(title, rot, tilt, weight, **kwargs)

    def plotHist(self, yValues, nbins, color='blue', **kwargs):
        """ Create a histogram. """
        # In some cases yValues is a generator, which cannot be indexed
        self.hist(list(yValues), nbins, facecolor=color, **kwargs)

    def plotScatter(self, xValues, yValues, color='blue', **kwargs):
        """ Create a scatter plot. """
        self.scatterP(xValues, yValues, c=color, **kwargs)

    def plotMatrix(self, img
                   , matrix
                   , vminData
                   , vmaxData
                   , cmap='jet'
                   , xticksLablesMajor=None
                   , yticksLablesMajor=None
                   , rotationX=90.
                   , rotationY=0.
                   , **kwargs):
        interpolation = kwargs.pop('interpolation', "none")
        plot = img.imshow(matrix, interpolation=interpolation, cmap=cmap,
                          vmin=vminData, vmax=vmaxData, **kwargs)
        if xticksLablesMajor is not None:
            plt.xticks(range(len(xticksLablesMajor)),
                       xticksLablesMajor[:len(xticksLablesMajor)],
                       rotation=rotationX)
        if yticksLablesMajor is not None:
            plt.yticks(range(len(yticksLablesMajor)),
                       yticksLablesMajor[:len(yticksLablesMajor)],
                       rotation=rotationY)
        return plot

    def plotData(self, xValues, yValues, color='blue', **kwargs):
        """ Shortcut function to plot some values.
        Params:
            xValues: list of values to show in x-axis
            yValues: list of values to show as values in y-axis
            color: color for the plot.
            **kwargs: keyword arguments that accepts:
                marker, linestyle
        """

        self.plot(xValues, yValues, color, **kwargs)

    def plotDataBar(self, xValues, yValues, width, color='blue', **kwargs):
        """ Shortcut function to plot some values.
        Params:
            xValues: list of values to show in x-axis
            yValues: list of values to show as values in y-axis
            color: color for the plot.
            **kwargs: keyword arguments that accepts:
                marker, linestyle
        """

        self.bar(xValues, yValues, width=width, color=color, **kwargs)

    @classmethod
    def createFromFile(cls, dbName, dbPreffix, plotType, columnsStr, colorsStr, linesStr,
                       markersStr, xcolumn, ylabel, xlabel, title, bins, orderColumn,
                       orderDirection):
        columns = columnsStr.split()
        colors = colorsStr.split()
        lines = linesStr.split()
        markers = markersStr.split()
        data = PlotData(dbName, dbPreffix, orderColumn, orderDirection)
        plotter = Plotter(windowTitle=title)
        ax = plotter.createSubPlot(title, xlabel, ylabel)
        xvalues = data.getColumnValues(xcolumn) if xcolumn else range(0, data.getSize())

        for i, col in enumerate(columns):
            yvalues = data.getColumnValues(col)
            color = colors[i]
            line = lines[i]
            colLabel = col if not col.startswith("_") else col[1:]
            if bins:
                yvalues = data._removeInfinites(yvalues)
                ax.hist(yvalues, bins=int(bins), color=color, linestyle=line, label=colLabel)
            else:
                if plotType == 'Plot':
                    marker = (markers[i] if not markers[i] == 'none' else None)
                    ax.plot(xvalues, yvalues, color, marker=marker, linestyle=line, label=colLabel)
                else:
                    ax.scatter(xvalues, yvalues, c=color, label=col, alpha=0.5)
        ax.legend()

        return plotter


class PlotData:
    """ Small wrapper around table data such as: sqlite or metadata
    files. """

    def __init__(self, fileName, tableName, orderColumn, orderDirection):
        self._orderColumn = orderColumn
        self._orderDirection = orderDirection

        if fileName.endswith(".db") or fileName.endswith(".sqlite"):
            self._table = self._loadSet(fileName, tableName)
            self.getColumnValues = self._getValuesFromSet
            self.getSize = self._table.getSize
        else:  # assume a metadata file
            self._table = self._loadMd(fileName, tableName)
            self.getColumnValues = self._getValuesFromMd
            self.getSize = self._table.size

    def _loadSet(self, dbName, dbPreffix):
        from pyworkflow.mapper.sqlite import SqliteFlatDb
        db = SqliteFlatDb(dbName=dbName, tablePrefix=dbPreffix)
        if dbPreffix:
            setClassName = "SetOf%ss" % db.getSelfClassName()
        else:
            setClassName = db.getProperty('self')  # get the set class name

        # FIXME: Check why the import is here
        from pwem import Domain
        setObj = Domain.getObjects()[setClassName](filename=dbName, prefix=dbPreffix)
        return setObj

    def _getValuesFromSet(self, columnName):
        return [self._getValue(obj, columnName)
                for obj in self._table.iterItems(orderBy=self._orderColumn,
                                                 direction=self._orderDirection)]

    @staticmethod
    def _removeInfinites(values):
        newValues = []
        for value in values:
            if isinstance(value, numbers.Number) and value < float("Inf"):
                newValues.append(value)
        return newValues

    def _loadMd(self, fileName, tableName):
        label = md.str2Label(self._orderColumn)
        tableMd = md.MetaData('%s@%s' % (tableName, fileName))
        tableMd.sort(label)  # FIXME: use order direction
        # TODO: sort metadata by self._orderColumn
        return tableMd

    def _getValuesFromMd(self, columnName):
        label = md.str2Label(columnName)
        return [self._table.getValue(label, objId) for objId in self._table]

    def _getValue(self, obj, column):
        if column == 'id':
            return obj.getObjId()

        return obj.getNestedValue(column)


# Functions for the angular distribution. Maybe they could go to other place?
def magnitude(x, y, z):
    """Returns the magnitude of the vector."""
    return sqrt(x * x + y * y + z * z)

def to_spherical(x, y, z):
    """Converts a cartesian coordinate (x, y, z) into a spherical one (radius, theta, phi) in radians.
        theta ranges from 0 t PI
        phi ranges from -PI to PI

    """
    radius = magnitude(x, y, z)
    theta = atan2(sqrt(x * x + y * y), z)
    phi = atan2(y, x)
    return (radius, theta, phi)

def eulerAngles_to_2D(rot, tilt, psi):
    """Converts euler angles to their 2D representation for a polar plot"""

    x, y, z = emlib.Euler_direction(rot,
                                    tilt,
                                    psi)
    # f.write(f".sphere {x} {y} {z} .01\n")
    # radius, theta, phi = to_spherical(x, y, z)
    ## may be radius, theta, phi = to_spherical(x, z, y)
    radius, theta, phi = to_spherical(y, z, x)
    if phi > pi:
        phi -= 2 * pi

    if phi < 0:
        phi = - phi
        theta += pi

    return theta, phi
