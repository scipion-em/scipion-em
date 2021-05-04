# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
Configuration definition for launching the ShowJ visualization utility.
The name of ShowJ is due historical reasons. This data visualization tools
was first called xmipp_show. And later was re-written using Java and 
became xmipp_showj.
"""

import os
import sys
from os.path import join
from collections import OrderedDict
import subprocess

from pyworkflow import Config
from pyworkflow.utils import getFreePort
from pwem import Config as emConfig
import threading
import shlex
import socketserver
from pyworkflow.constants import SCIPION_DOMAIN, SCIPION_HOME_VAR, SCIPION_DEBUG


# ----------------------- Showj constants ---------------------------
COL_RENDER_NONE = 0
COL_RENDER_ID = 1
COL_RENDER_TEXT = 2
COL_RENDER_IMAGE = 3
COL_RENDER_CHECKBOX = 4
COL_RENDER_VOLUME = 5

PATH = 'path'
DATASET = 'dataset'

TABLE_NAME = 'blockComboBox'
COLS_CONFIG = 'tableLayoutConfiguration'
COLS_CONFIG_DEFAULT = 'defaultColumnsLayoutProperties'
LABEL_SELECTED = 'labelsToRenderComboBox'

MODE = 'mode'
MODE_GALLERY = 'gallery'
MODE_TABLE = 'table'
MODE_MD = 'metadata'
MODE_VOL_ASTEX = 'volume_astex'
MODE_VOL_CHIMERA = 'volume_chimera'
MODE_VOL_JSMOL = 'volume_jsmol'
RENDER = 'render'
ORDER = 'order'
VISIBLE = 'visible'
ZOOM = 'zoom'
SORT_BY = 'sortby'
LABELS = 'labels'
CLASSIFIER = 'classifier'

SAMPLINGRATE = 'sampling_rate'
CHIMERA_PORT = 'chimera_port'
INVERTY = 'inverty'

OBJCMDS = 'object_commands'

GOTO = 'goto'
ROWS = 'rows'
COLS = 'cols'
ALLOW_RENDER = 'allowRender'
MANUAL_ADJUST = 'colRowMode'

SELECTEDITEMS = 'listSelectedItems'
ENABLEDITEMS = 'listEnabledItems'
CHANGES = 'listChangesItems'
OLDMODE = 'oldModeShowj'
RENDERITEMS = 'renderItems'

VOL_SELECTED = 'volumesToRenderComboBox'
VOL_TYPE = 'typeVolume'
VOL_VIEW = 'resliceComboBox'

IMG_DIMS = 'imageDimensions'
IMG_ZOOM = 'zoom'
IMG_ZOOM_DEFAULT = 'defaultZoom'
IMG_MIRRORY = 'mirrorY'
IMG_APPLY_TRANSFORM = 'applyTransformMatrix'
IMG_ONLY_SHIFTS = 'onlyShifts'
IMG_WRAP = 'wrap'
IMG_MAX_WIDTH = 'imageMaxWidth'
IMG_MIN_WIDTH = 'imageMinWidth'
IMG_MAX_HEIGHT = 'imageMaxHeight'
IMG_MIN_HEIGHT = 'imageMinHeight'

PROJECT_NAME = 'projectName'
PROJECT_PATH = 'projectPath'
OBJECT_ID = 'objectId'


def getJvmMaxMemory():
    # Memory in GB
    return int(emConfig.JAVA_MAX_MEMORY)


class ColumnsConfig:
    """ Store the configuration of the columns for a given table in a dataset.
    The order of the columns will be stored and configuration for each columns.
    For each column, we store properties:
    - visible
    - allowSetVisible
    - renderable
    - allowSetRenderable
    - editable
    - allowSetEditable
    - renderFunc
    - renderFuncExtra
    """
    def __init__(self, table, allowRender=True, defaultColumnsLayout={}):

        self._columnsDict = OrderedDict()

        for col in table.iterColumns():
            colDefaultLayout = defaultColumnsLayout.get(col.getLabel(), {})
            col_properties = ColumnProperties(col, allowRender, colDefaultLayout)
            self._columnsDict[col.getName()] = col_properties

    def getRenderableColumns(self, extra=None):
        """ Return a list with the name of renderable columns.
            extra: parameter used to keep some rendering columns."""
        columnsName = []
        columnsLabel = []

        for col in self._columnsDict.values():
            if col.isRenderable():
                if extra is not None:
                    if col.getName() in extra:
                        columnsName.append(col.getName())
                        columnsLabel.append(col.getLabel())

                else:
                    columnsName.append(col.getName())
                    columnsLabel.append(col.getLabel())
#       columns = [col.getName() for col in self._columnsDict.values() if col.isRenderable()]
        return columnsName, columnsLabel

    def hasEnableColumn(self):
        for columnLayout in self._columnsDict.values():
            if "enable" == columnLayout.getLabel():
                return True
        return False

    def getColumnProperty(self, colName, propName):
        """ Get some property value of a given column. """
        col = self._columnsDict[colName]
        return getattr(col, propName)

    def configColumn(self, colName, **kwargs):
        """ Configure properties of a given column. """
        col = self._columnsDict[colName]
        for k, v in kwargs.items():
            setattr(col, k, v)

    def printColumns(self):
        for col in self._columnsDict.values():
            print("column: ", col.getLabel())
            print("  values: ", col.getValues())


class ColumnProperties:
    """ Store some properties to customize how each column
    will be display in the table.
    """
    def __init__(self, col, allowRender, defaultColumnLayoutProperties):
        self._column = col
        self.columnType = col.getRenderType()

        if 'visible' in defaultColumnLayoutProperties:
            self.visible = defaultColumnLayoutProperties['visible']
        else:
            self.visible = not (self.columnType == COL_RENDER_ID)

        self.allowSetVisible = True

        self.editable = (self.columnType == COL_RENDER_TEXT)
        self.allowSetEditable = self.editable

        self.renderable = ('renderable' in defaultColumnLayoutProperties
                           and defaultColumnLayoutProperties['renderable'].lower() == 'true')

        self.allowSetRenderable = ((self.columnType == COL_RENDER_IMAGE or
                                   self.columnType == COL_RENDER_VOLUME) and allowRender)

        self.renderFunc = "get_image"
        self.extraRenderFunc = ""

    def getLabel(self):
        return self._column.getLabel()

    def getName(self):
        return self._column.getName()

    def getColumnType(self):
        return self.columnType

    def allowsRenderable(self):
        return self.renderable or self.allowSetRenderable

    def isRenderable(self):
        return self.renderable or self.allowSetRenderable

    def setValues(self, defaultColumnLayoutProperties):
        for key in defaultColumnLayoutProperties:
            setattr(self, key, defaultColumnLayoutProperties[key])

    def getValues(self):
        return {"visible": self.visible,
                "allowSetVisible": self.allowSetVisible,
                "editable": self.editable,
                "allowSetEditable": self.allowSetEditable,
                "renderable": self.renderable,
                "allowSetRenderable": self.allowSetRenderable,
                "renderFunc": self.renderFunc,
                "extraRenderFunc": self.extraRenderFunc,
                'columnType': self.columnType,
                'columnName': self.getName(),
                'columnLabel': self.getLabel()
                }


def getArchitecture():
    import platform
    arch = platform.architecture()[0]
    for a in ['32', '64']:
        if a in arch:
            return a
    return 'NO_ARCH'


def getJavaIJappArguments(memory, appName, appArgs):
    """ Build the command line arguments to launch
    a Java application based on ImageJ.
    memory: the amount of memory passed to the JVM.
    appName: the qualified name of the java application.
    appArgs: the arguments specific to the application.
    """
    if memory is None:
        memory = getJvmMaxMemory()

    jdkLib = join(emConfig.JAVA_HOME, 'lib')
    javaBind = join(emConfig.XMIPP_HOME, "bindings", "java")
    imagej_home = join(javaBind, "imagej")
    lib = join(emConfig.XMIPP_HOME, "lib")
    javaLib = join(javaBind, 'lib')
    plugins_dir = os.path.join(imagej_home, "plugins")
    arch = getArchitecture()

    import subprocess
    version = subprocess.check_output(['java', '-version'], stderr=subprocess.STDOUT)
    majorVersion = int(version.decode("utf-8").split('"')[1].split('.')[0])
    if majorVersion > 9:
        args = "-Xmx%(memory)sg -Djava.library.path=%(lib)s -Dplugins.dir=%(plugins_dir)s -cp %(jdkLib)s/*:%(imagej_home)s/*:%(javaLib)s/* %(appName)s %(appArgs)s" % locals()
    else:
        args = "-Xmx%(memory)sg -d%(arch)s -Djava.library.path=%(lib)s -Dplugins.dir=%(plugins_dir)s -cp %(jdkLib)s/*:%(imagej_home)s/*:%(javaLib)s/* %(appName)s %(appArgs)s" % locals()

    return args


def runJavaIJapp(memory, appName, args, env=None):
    env = env or {}
    getEnviron = Config.getDomain().importFromPlugin('xmipp3', 'Plugin').getEnviron
    env.update(getEnviron(xmippFirst=False))

    # Add scipion variables
    env[SCIPION_DOMAIN] = Config.SCIPION_DOMAIN
    env[SCIPION_HOME_VAR] = Config.SCIPION_HOME

    if Config.debugOn():
        env[SCIPION_DEBUG] = "1"
    
    # Add out python path to the PATH
    env["PATH"] = ":".join([sys.executable, env.get("PATH", "")])
    # Try chimera
    try:
        chimeraPlugin = Config.getDomain().importFromPlugin('chimera', 'Plugin', doRaise=True)
        env["CHIMERA_LAUNCHER"] = chimeraPlugin.getProgram() + " %s"
    except Exception as e:
        pass

    args = getJavaIJappArguments(memory, appName, args)
    print('java %s' % args)
    # return subprocess.Popen('java ' + args, shell=True, env=env)
    cmd = ['java'] + shlex.split(args)
    return subprocess.Popen(cmd, env=env)


def launchSupervisedPickerGUI(micsFn, outputDir, protocol,
                              mode=None, memory=None,
                              pickerProps=None, inTmpFolder=False):
    app = "xmipp.viewer.particlepicker.training.SupervisedPickerRunner"
    args = "--input %s --output %s" % (micsFn, outputDir)

    if mode:
        args += " --mode %s" % mode

    if pickerProps:
        args += " --classifier " + pickerProps
    else:
        port = initProtocolTCPServer(protocol)
        args += " --scipion %s" % port

    if inTmpFolder:
        args += " --tmp true"

    return runJavaIJapp(memory, app, args)


def launchTiltPairPickerGUI(micsFn, outputDir, protocol, mode=None):
    port = initProtocolTCPServer(protocol)
    app = "xmipp.viewer.particlepicker.tiltpair.TiltPairPickerRunner"
    args = "--input %s --output %s  --scipion %s" % (micsFn, outputDir, port)
    if mode:
        args += " --mode %s" % mode
    return runJavaIJapp(None, app, args)


class ProtocolTCPRequestHandler(socketserver.BaseRequestHandler):

    def handle(self):  # FIXME: RUNNING FOREVER
        protocol = self.server.protocol
        msg = str(self.request.recv(1024), "utf-8")
        tokens = shlex.split(msg)
        if msg.startswith('run function'):
            functionName = tokens[2]
            # try:
            functionPointer = getattr(protocol, functionName)
            functionPointer(*tokens[3:])
            self.request.sendall(b'done\n')
            self.server.end = True
            # except:
            # print 'protocol %s must implement %s'%(protocol.getName(), functionName)

        else:
            answer = b'no answer available\n'
            self.request.sendall(answer)


class MySocketServer(socketserver.TCPServer):

    def serve_forever(self):
        self.end = False
        while not self.end:
            self.handle_request()


def initProtocolTCPServer(protocol):
    address = ''
    port = getFreePort()
    server = MySocketServer((address, port), ProtocolTCPRequestHandler)
    server.protocol = protocol
    server_thread = threading.Thread(name="ProtocolTCPServer",
                                     target=server.serve_forever)
    server_thread.daemon = True
    server_thread.start()
    return port
