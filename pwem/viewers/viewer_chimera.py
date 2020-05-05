# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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

import os
from multiprocessing.connection import Client
import threading

import pyworkflow.utils as pwutils
import pyworkflow.viewer as pwviewer
import pwem.constants as emcts
import pwem.emlib.metadata as md
import pwem.objects as emobj
from pwem import emlib
from pwem import getCmdPath, Config as emConfig


chimeraPdbTemplateFileName = "Atom_struct__%s.pdb"
chimeraMapTemplateFileName = "Map__%s.mrc"
chimeraScriptFileName = "chimeraScript.py"
sessionFile = "SESSION.py"

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
        return emConfig.CHIMERA_HOME

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
    def getProgram(cls, progName="chimera"):
        """ Return the program binary that will be used. """
        home = cls.getHome()
        if home is None:
            return None
        return os.path.join(home, 'bin', os.path.basename(progName))

    @classmethod
    def runProgram(cls, program=None, args=""):
        """ Internal shortcut function to launch chimera program. """
        prog = program or cls.getProgram()
        pwutils.runJob(None, prog, args, env=cls.getEnviron())

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
            sampling * dim * 3. / 4.
        arrowDict["r1"] = r1 * dim / 50.
        arrowDict["r2"] = 4 * r1
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
    """ Debug funtion. """
    pass
    # print cmd


class ChimeraClient:

    def openVolumeOnServer(self, volume, sendEnd=True):
        self.send('open_volume', volume)
        if self.voxelSize is not None:
            self.send('voxel_size', self.voxelSize)
        if sendEnd:
            self.client.send('end')

    def __init__(self, volfile, sendEnd=True, **kwargs):
        if volfile.endswith('.mrc'):
            volfile += ':mrc'

        self.kwargs = kwargs
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

        self.volfile = volfile
        self.voxelSize = self.kwargs.get('voxelSize', None)
        # ChimeraServer is the basic server. There are other
        # than inherit from it.
        serverName = self.kwargs.get('ChimeraServer', 'ChimeraServer')
        self.address = ''
        self.port = pwutils.getFreePort()

        serverfile = getCmdPath('chimera_server.py')
        command = pwviewer.CommandView("chimera --script '%s %s %s' &" %
                                       (serverfile, self.port, serverName),
                                       env=Chimera.getEnviron(), ).show()
        self.authkey = 'test'
        self.client = Client((self.address, self.port), authkey=self.authkey)
        self.initVolumeData()
        # self.openVolumeOnServer(self.vol,sendEnd)
        self.openVolumeOnServer(self.vol)
        self.initListenThread()

    def send(self, cmd, data):
        self.client.send(cmd)
        self.client.send(data)

    def initListenThread(self):
        self.listen_thread = threading.Thread(name="ChimeraCli.listenTh",
                                              target=self.listen)
        # self.listen_thread.daemon = True
        self.listen_thread.start()

    def listen(self):

        self.listen = True
        try:
            while self.listen:
                msg = self.client.recv()
                self.answer(msg)
        except EOFError:
            print('Lost connection to server')
        finally:
            self.exit()

    def exit(self):
        self.client.close()

    def initVolumeData(self):
        self.image = emlib.Image(self.volfile)
        self.image.convert2DataType(md.DT_DOUBLE)
        self.xdim, self.ydim, self.zdim, self.n = self.image.getDimensions()
        self.vol = self.image.getData()

    def answer(self, msg):
        if msg == 'exit_server':
            self.listen = False


class ChimeraAngDistClient(ChimeraClient):

    def __init__(self, volfile, **kwargs):
        self.angularDistFile = kwargs.get('angularDistFile', None)

        if self.angularDistFile:
            fileName = self.angularDistFile
            if '@' in self.angularDistFile:
                fileName = self.angularDistFile.split("@")[1]
            if not (os.path.exists(fileName)):  # check blockname:
                raise Exception("Path %s does not exists" %
                                self.angularDistFile)
        self.spheresColor = kwargs.get('spheresColor', 'red')
        spheresDistance = kwargs.get('spheresDistance', None)
        spheresMaxRadius = kwargs.get('spheresMaxRadius', None)
        ChimeraClient.__init__(self, volfile, sendEnd=False, **kwargs)
        self.spheresDistance = float(spheresDistance) \
            if spheresDistance else 0.75 * max(self.xdim, self.ydim, self.zdim)
        self.spheresMaxRadius = float(spheresMaxRadius) \
            if spheresMaxRadius else 0.02 * self.spheresDistance
        self.loadAngularDist(True)

    def loadAngularDist(self, sendEnd=True):
        if self.angularDistFile:
            self.readAngularDistFile()
            self.send('command_list', self.angulardist)
        if sendEnd:
            self.client.send('end')

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

        self.angulardist = []
        x2 = self.xdim / 2
        y2 = self.ydim / 2
        z2 = self.zdim / 2
        # cofr does not seem to work!
        # self.angulardist.append('cofr %d,%d,%d'%(x2,y2,z2))
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

            x = x * self.spheresDistance + x2
            y = y * self.spheresDistance + y2
            z = z * self.spheresDistance + z2
            command = 'shape sphere radius %s center %s,%s,%s color %s ' % \
                      (radius, x, y, z, self.spheresColor)
            self.angulardist.append(command)
            # printCmd(command)



class ChimeraView(pwviewer.CommandView):
    """ View for calling an external command. """

    def __init__(self, inputFile, **kwargs):
        pwviewer.CommandView.__init__(self, 'chimera "%s" &' % inputFile,
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
                # check if tmp dir exists, if not use /tmp
                # tmp does not exists if you try to visualize something  (eye)
                # before irunning the protocol
                tmpPath = self.protocol._getTmpPath()
                if not os.path.exists(tmpPath):
                    tmpPath = "/tmp"
                fnCmd = os.path.join(tmpPath, "chimera.cmd")
                f = open(fnCmd, 'w')
                f.write("cofr 0,0,0\n")  # set center of coordinates
                if obj.hasVolume():
                    volID = 0
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
                    os.path.join(tmpPath, "axis.bild"))
                Chimera.createCoordinateAxisFile(dim,
                                                 bildFileName=bildFileName,
                                                 sampling=sampling)
                f.write("open %s\n" % bildFileName)
                f.write("open %s\n" % os.path.abspath(fn))
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
