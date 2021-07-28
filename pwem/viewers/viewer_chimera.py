# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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

import os
from multiprocessing.connection import Client
import threading

import pyworkflow.utils as pwutils
import pyworkflow.viewer as pwviewer
import pwem.constants as emcts
import pwem.emlib.metadata as md
import pwem.objects as emobj
from pwem import emlib
from pwem import Config as emConfig

chimeraPdbTemplateFileName = "Atom_struct__%06d.cif"
chimeraMapTemplateFileName = "Map__%06d.mrc"
chimeraScriptFileName = "chimeraScript.cxc"
chimeraConfigFileName = "chimera.ini"
sessionFile = "SESSION.cxs"

symMapperScipionchimera = {}
symMapperScipionchimera[emcts.SYM_CYCLIC] = "Cn"
symMapperScipionchimera[emcts.SYM_DIHEDRAL] = "Dn"
symMapperScipionchimera[emcts.SYM_TETRAHEDRAL] = "T"
symMapperScipionchimera[emcts.SYM_OCTAHEDRAL] = "O"
symMapperScipionchimera[emcts.SYM_I222] = "222"
symMapperScipionchimera[emcts.SYM_I222r] = "222r"
symMapperScipionchimera[emcts.SYM_In25] = "n25"
symMapperScipionchimera[emcts.SYM_In25r] = "n25r"
symMapperScipionchimera[emcts.SYM_I2n3] = "2n3"
symMapperScipionchimera[emcts.SYM_I2n3r] = "2n3r"
symMapperScipionchimera[emcts.SYM_I2n5] = "2n5"
symMapperScipionchimera[emcts.SYM_I2n5r] = "2n5r"


class Chimera:
    """ Helper class to execute chimera and handle its environment. """

    # Map symmetries from Scipion convention to Chimera convention
    _symmetryMap = {
        emcts.SYM_CYCLIC: 'Cn',
        emcts.SYM_DIHEDRAL: 'Dn',
        emcts.SYM_TETRAHEDRAL: 'T',
        emcts.SYM_OCTAHEDRAL: 'O',
        emcts.SYM_I222: '222',
        emcts.SYM_I222r: '222r',
        emcts.SYM_In25: 'n25',
        emcts.SYM_In25r: 'n25r',
        emcts.SYM_I2n3: '2n3',
        emcts.SYM_I2n3r: '2n3r',
        emcts.SYM_I2n5: '2n5',
        emcts.SYM_I2n5r: '2n5r'
    }

    @classmethod
    def getSymmetry(cls, scipionSym):
        """ Return the equivalent Chimera symmetry from Scipion one. """
        return cls._symmetryMap[scipionSym]

    @classmethod
    def getHome(cls):
        """ Returns chimera home, trying first to take it from chimera plugin. If it fails it will return, the default value in the config"""
        with pwutils.weakImport("chimera"):
            from chimera import Plugin as chimeraPlugin
            return chimeraPlugin.getHome()

        return os.path.join(emConfig.EM_ROOT, 'chimerax-1.2.5')

    @classmethod
    def getEnviron(cls):
        """ Return the proper environ to launch chimera.
        CHIMERA_HOME variable is read from the ~/.config/scipion.conf file.
        """
        environ = pwutils.Environ(os.environ)
        environ.set('PATH', os.path.join(cls.getHome(), 'bin'),
                    position=pwutils.Environ.BEGIN)

        if "REMOTE_MESA_LIB" in os.environ:
            environ.set('LD_LIBRARY_PATH', os.environ['REMOTE_MESA_LIB'],
                        position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def getProgram(cls, progName="ChimeraX"):
        """ Return the program binary that will be used. """
        home = cls.getHome()
        if home is None:
            return None
        return os.path.join(home, 'bin', os.path.basename(progName))

    @classmethod
    def runProgram(cls, program=None, args="", cwd=None):
        """ Internal shortcut function to launch chimera program. """
        prog = program or cls.getProgram()
        pwutils.runJob(None, prog, args, env=cls.getEnviron(),
                       cwd=cwd)

    @classmethod
    def createCoordinateAxisFile(cls, dim, bildFileName="/tmp/axis.bild",
                                 sampling=1, r1=0.1):
        """ Create a coordinate system, Along each dimension
        we place a small sphere in the negative axis. In this way
        chimera shows the system of coordinates origin in the
        window center"""

        ff = open(bildFileName, "w")
        arrowDict = {}
        arrowDict["x"] = arrowDict["y"] = arrowDict["z"] = \
            sampling * dim * 1. / 2.
        arrowDict["r1"] = r1 * dim / 50.
        arrowDict["r2"] = 2 * arrowDict["r1"]
        arrowDict["rho"] = 0.75  # axis thickness

        ff.write(".color 1 0 0\n"
                 ".arrow 0 0 0 %(x)f 0 0 %(r1)f %(r2)f %(rho)f\n"
                 ".color 0 0 0\n.sphere -%(x)f 0 0 0.00001\n"
                 ".color 1 1 0\n"
                 ".arrow 0 0 0 0 %(y)f 0 %(r1)f %(r2)f %(rho)f\n"
                 ".color 0 0 0\n.sphere 0 -%(y)f 0 0.00001\n"
                 ".color 0 0 1\n"
                 ".arrow 0 0 0 0 0 %(z)f %(r1)f %(r2)f %(rho)f\n"
                 ".color 0 0 0\n.sphere 0 0 -%(z)f 0.00001\n" %
                 arrowDict)
        ff.close()


def printCmd(cmd):
    """ Debug function. """
    pass
    # print cmd


class ChimeraAngDist(pwviewer.CommandView):

    def __init__(self, volFile, tmpFilesPath, **kwargs):
        self.kwargs = kwargs
        self.volfile = os.path.abspath(os.path.join(volFile))
        self.cmdFile = os.path.join(tmpFilesPath,
                                    "chimera_angular_distribution.cxc")
        self.axis = os.path.abspath(os.path.join(tmpFilesPath, "axis.bild"))
        self.spheres = os.path.abspath(os.path.join(tmpFilesPath,
                                                    "spheres.bild"))

        self.angularDistributionList = []
        self.angularDistFile = kwargs.get('angularDistFile', None)
        if self.angularDistFile:
            fileName = self.angularDistFile
            if '@' in self.angularDistFile:
                fileName = self.angularDistFile.split("@")[1]
            if not (os.path.exists(fileName)):  # check blockname:
                raise Exception("Path %s does not exists" %
                                self.angularDistFile)

        self.volOrigin = kwargs.get('volOrigin', None)
        self.spheresColor = kwargs.get('spheresColor', 'orange')
        spheresDistance = kwargs.get('spheresDistance', -1)
        spheresMaxRadius = kwargs.get('spheresMaxRadius', None)
        self.initVolumeData(self.volfile)
        self.spheresDistance = (float(spheresDistance) if spheresDistance != -1
                                else (0.75 * self.voxelSize * max(self.xdim,
                                                                  self.ydim,
                                                                  self.zdim)) / 2)
        self.spheresMaxRadius = (float(spheresMaxRadius) if spheresMaxRadius
                                 else 0.02 * self.spheresDistance)
        self.readAngularDistFile()
        format = kwargs.get('format', None)
        self._createChimeraScript(self.cmdFile, format)
        program = Chimera.getProgram()
        pwviewer.CommandView.__init__(self, '%s "%s" &' % (program,
                                                           self.cmdFile),
                                      env=Chimera.getEnviron(), **kwargs)

    def _createChimeraScript(self, scriptFile, format):

        Chimera.createCoordinateAxisFile(self.xdim,
                                         bildFileName=self.axis,
                                         sampling=self.voxelSize)
        self.createAngularDistributionFile(bildFileName=self.spheres)
        fhCmd = open(scriptFile, 'w')
        fhCmd.write("open %s\n" % self.axis)
        fhCmd.write("open %s\n" % self.spheres)
        fhCmd.write("cofr 0,0,0\n")
        if format is not None:
            fhCmd.write("open %s format %s\n" % (self.volfile.replace(":mrc", ""),
                                                 format))
        else:
            fhCmd.write("open %s\n" % self.volfile.replace(":mrc", ""))
        fhCmd.write("volume #3 voxelSize %s\n" % self.voxelSize)
        x, y, z = self.volOrigin
        fhCmd.write("volume #3 origin %s,%s,%s\n" % (x, y, z))
        fhCmd.write("view\n")

        fhCmd.close()

    def createAngularDistributionFile(self, bildFileName="/tmp/spheres.bild"):
        ff = open(bildFileName, "w")

        for angulardist in self.angularDistributionList:
            ff.write("%s\n" % angulardist)

        ff.close()

    def initVolumeData(self, volfile):
        if volfile.endswith('.mrc'):
            volfile += ':mrc'

        if volfile is None:
            raise ValueError(volfile)
        if '@' in volfile:
            [index, file] = volfile.split('@')
        else:
            file = volfile
        if ':' in file:
            file = file[0: file.rfind(':')]
        if not os.path.exists(file):
            raise Exception("File %s does not exists" % file)

        self.voxelSize = self.kwargs.get('voxelSize', 1.0)
        self.image = emlib.Image(self.volfile)
        self.image.convert2DataType(md.DT_DOUBLE)
        self.xdim, self.ydim, self.zdim, self.n = self.image.getDimensions()
        self.vol = self.image.getData()

    def readAngularDistFile(self):
        angleRotLabel = md.MDL_ANGLE_ROT
        angleTiltLabel = md.MDL_ANGLE_TILT
        anglePsiLabel = md.MDL_ANGLE_PSI
        mdAngDist = md.MetaData(self.angularDistFile)
        if not mdAngDist.containsLabel(md.MDL_ANGLE_PSI):
            anglePsiLabel = None
            if mdAngDist.containsLabel(md.RLN_ORIENT_PSI):
                angleRotLabel = md.RLN_ORIENT_ROT
                angleTiltLabel = md.RLN_ORIENT_TILT
                anglePsiLabel = md.RLN_ORIENT_PSI

        if not mdAngDist.containsLabel(md.MDL_WEIGHT):
            mdAngDist.fillConstant(md.MDL_WEIGHT, 1.)

        maxweight = mdAngDist.aggregateSingle(md.AGGR_MAX, md.MDL_WEIGHT)
        minweight = mdAngDist.aggregateSingle(md.AGGR_MIN, md.MDL_WEIGHT)
        interval = maxweight - minweight
        x2 = self.xdim / 2
        y2 = self.ydim / 2
        z2 = self.zdim / 2
        self.angularDistributionList.append('.color %s\n' % self.spheresColor)
        for id in mdAngDist:
            rot = mdAngDist.getValue(angleRotLabel, id)
            tilt = mdAngDist.getValue(angleTiltLabel, id)
            psi = mdAngDist.getValue(anglePsiLabel, id) if anglePsiLabel else 0
            weight = mdAngDist.getValue(md.MDL_WEIGHT, id)
            # Avoid zero division
            weight = 0 if interval == 0 else (weight - minweight) / interval
            weight = weight + 0.5  # add 0.5 to avoid zero weight
            x, y, z = emlib.Euler_direction(rot, tilt, psi)
            radius = weight * self.spheresMaxRadius

            x = x * self.spheresDistance  # + x2
            y = y * self.spheresDistance  # + y2
            z = z * self.spheresDistance  # + z2
            command = ('.sphere %.1f %.1f %.1f %.1f' %
                       (x, y, z, radius))
            self.angularDistributionList.append(command)


class ChimeraView(pwviewer.CommandView):
    """ View for calling an external command. """

    def __init__(self, inputFile, **kwargs):
        program = Chimera.getProgram()
        inputFile = inputFile.replace(":mrc", "")
        pwviewer.CommandView.__init__(self, '%s "%s" &' % (program, inputFile),
                                      env=Chimera.getEnviron(), **kwargs)


class ChimeraClientView(ChimeraView):
    """ View for calling an external command. """
    pass


class ChimeraDataView(ChimeraView):

    def __init__(self, dataview, vol, viewParams={}, **kwargs):
        self.dataview = dataview
        ChimeraClientView.__init__(self, vol.getFileName())

    def show(self):
        self.dataview.show()
        ChimeraClientView.show(self)


class ChimeraViewer(pwviewer.Viewer):
    """ Wrapper to visualize PDB object with Chimera. """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [emobj.AtomStruct, emobj.PdbFile]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)

    def visualize(self, obj, **kwargs):
        cls = type(obj)
        if issubclass(cls, emobj.AtomStruct):
            # if attribute _chimeraScript exists then protocol
            # has create a script file USE IT
            if hasattr(obj, '_chimeraScript'):
                fn = obj._chimeraScript.get()
                ChimeraView(fn).show()
                return
            # if not create a script file with: coordinates axis, PDB and
            # volume (if available)
            else:
                fn = obj.getFileName()
                extraPath = self.protocol._getExtraPath()
                fnCmd = os.path.join(extraPath, "chimera.cxc")
                f = open(fnCmd, 'w')
                f.write("cofr 0,0,0\n")  # set center of coordinates
                if obj.hasVolume():
                    volID = 1
                    volumeObject = obj.getVolume()
                    dim = volumeObject.getDim()[0]
                    sampling = volumeObject.getSamplingRate()
                    f.write("open %s\n" % os.path.abspath(
                        emlib.image.ImageHandler.removeFileType(volumeObject.getFileName())))
                    f.write("volume #%d style surface voxelSize %f\n"
                            % (volID, sampling))
                    x, y, z = volumeObject.getShiftsFromOrigin()
                    f.write("volume #%d origin %0.2f,%0.2f,%0.2f\n"
                            % (volID, x, y, z))
                else:
                    dim = 150  # eventually we will create a PDB library that
                    # computes PDB dim
                    sampling = 1.
                # Construct the coordinate file
                bildFileName = os.path.abspath(
                    os.path.join(extraPath, "axis.bild"))
                Chimera.createCoordinateAxisFile(dim,
                                                 bildFileName=bildFileName,
                                                 sampling=sampling)
                f.write("open %s\n" % bildFileName)
                f.write("open %s\n" % os.path.abspath(fn))
                f.write("style stick\n")
                f.write("view\n")
                f.close()
                ChimeraView(fnCmd).show()
            # FIXME: there is an asymmetry between ProtocolViewer and Viewer
            # for the first, the visualize method return a list of View's
            # (that are shown)
            # for the second, the visualize directly shows the objects.
            # the first approach is better
        else:
            raise Exception('ChimeraViewer.visualize: '
                            'can not visualize class: %s' % obj.getClassName())


def mapVolsWithColorkey(displayVolFileName,
                        mapVolFileName,
                        stepColors,  # List with resolution values
                        colorList,  # List with colors
                        voldim,  # pixels
                        volOrigin=None,  # in A.
                        step=-1,  # chimera volume step (display resolution)
                        sampling=1.,
                        scriptFileName='/tmp/chimeraColor.py',
                        bgColorImage='white',  # background color
                        showAxis=True,
                        fontSize=12):  # default chimera font size
    """ colors surface of volume 'displayVolFileName' using values from
    'mapVolFikeName'. A colorkey is created using the values in stepColors
    and the color in colorList. """
    scriptFile = scriptFileName
    fhCmd = open(scriptFile, 'w')
    fhCmd.write("from chimerax.core.commands import run\n")
    fhCmd.write("from chimerax.graphics.windowsize import window_size\n")
    fhCmd.write("from PyQt5.QtGui import QFontMetrics\n")
    fhCmd.write("from PyQt5.QtGui import QFont\n")
    fhCmd.write("run(session, 'set bgColor %s')\n" % bgColorImage)

    chimeraVolId = 1

    bildFileName = scriptFile.replace(".py", ".bild")
    Chimera.createCoordinateAxisFile(voldim[0],
                                     bildFileName=bildFileName,
                                     sampling=sampling)
    # axis
    fhCmd.write("run(session, 'open %s')\n" % bildFileName)
    if not showAxis:
        fhCmd.write("run(session, 'hide #1')\n")

    fhCmd.write("run(session, 'cofr 0,0,0')\n")  # set center of coordinates

    chimeraVolId += 1

    # first volume
    if volOrigin is None:

        x = -voldim[0] * sampling // 2
        y = -voldim[1] * sampling // 2
        z = -voldim[2] * sampling // 2
    else:
        # TODO, not sure about sign
        x = volOrigin[0]
        y = volOrigin[1]
        z = volOrigin[2]

    fhCmd.write("run(session, 'open %s')\n" % displayVolFileName)
    if step == -1:
        fhCmd.write("run(session, 'volume #%d voxelSize %s')\n" %
                    (chimeraVolId, str(sampling)))
    else:
        fhCmd.write("run(session, 'volume #%d voxelSize %s step %d')\n" %
                    (chimeraVolId, str(sampling), step))
    fhCmd.write("run(session, 'volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                % (chimeraVolId, x, y, z))
    # second volume
    chimeraVolId += 1
    fhCmd.write("run(session, 'open %s')\n" % mapVolFileName)
    # TODO: Why no step here?
    fhCmd.write("run(session, 'volume #%d voxelSize %f')\n" % (chimeraVolId, sampling))
    fhCmd.write("run(session, 'volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                % (chimeraVolId, x, y, z))
    fhCmd.write("run(session, 'vol #%d hide')\n" % chimeraVolId)

    # replace scolor + colorkey which has been discontinued in chimerax
    scolorStr = ''
    for step, color in zip(stepColors, colorList):
        scolorStr += '%s,%s:' % (step, color)
    scolorStr = scolorStr[:-1]
    fhCmd.write("run(session, 'color sample #%d map #%d "
                "palette " % (chimeraVolId - 1, chimeraVolId) +
                scolorStr + "')\n"
                )
    # get window size so we can place labels properly

    fhCmd.write("v = session.main_view\n")
    fhCmd.write("vx,vy=v.window_size\n")

    # Calculate heights and Y positions: font, scale height and firstY
    ptSize = fontSize  # default size chimera font
    fhCmd.write('font = QFont("Ariel", %d)\n' % ptSize)
    fhCmd.write('f = QFontMetrics(font)\n')
    fhCmd.write('_height =  1 * f.height()/vy\n')  # Font height
    fhCmd.write('_half_scale_height = _height * %f/2\n' % len(stepColors))  # Full height of the scale
    fhCmd.write("_firstY= 0.5 + _half_scale_height\n")  # Y location for first label

    fhCmd.write("step = ")
    # place labels in right place
    # unfortunately chimera has no colorbar
    labelCount = 0
    for step, color in zip(stepColors, colorList):
        if step > 99.9:
            step = "%.2f" % step
        else:
            step = "{:05.2f}".format(step)
        command = 'run(session, "2dlabel text ' + step + \
                  ' bgColor ' + color + \
                  ' xpos 0.01 ypos %f' + \
                  ' size ' + str(ptSize) + \
                  '" % ' + \
                  '(_firstY - %f*_height))\n' % (labelCount)

        fhCmd.write(command)
        labelCount += 1
    fhCmd.write("run(session, 'view')\n")
    fhCmd.close()
