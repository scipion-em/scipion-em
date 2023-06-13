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
"""
This modules contains classes related with EM
"""

import os
from pkg_resources import parse_version

import pyworkflow as pw
from pyworkflow.protocol import Protocol
from pyworkflow.utils import  getSubclasses
from pyworkflow.viewer import Viewer
from pyworkflow.wizard import Wizard
import pyworkflow.plugin

from .constants import *
from .objects import EMObject
from .tests import defineDatasets
from .utils import *

__version__ = '3.1.0'
NO_VERSION_FOUND_STR = "0.0"
CUDA_LIB_VAR = 'CUDA_LIB'

_logo = "scipion_icon.png"
_references = ["delaRosaTrevin201693"]


class Config(pw.Config):
    _get = pw.Config._get
    _join = pw.Config._join

    EM_ROOT = _join(_get(EM_ROOT_VAR, _join(pw.Config.SCIPION_SOFTWARE, 'em')))

    # Default XMIPP_HOME: needed here for ShowJ viewers
    XMIPP_HOME = _join(_get('XMIPP_HOME', os.path.join(EM_ROOT, 'xmipp')))

    # Get java home, we might need to provide correct default value. Use SCIPION_JAVA_HOME to force it when there is other JAVA_HOME you don't want/cant to change: e.g. pycharm debugging.
    JAVA_HOME = _get('SCIPION_JAVA_HOME', _get('JAVA_HOME', ''))
    JAVA_MAX_MEMORY = _get('JAVA_MAX_MEMORY', '4')

    # MPI
    MPI_LIBDIR = _get('MPI_LIBDIR', '/usr/lib64/mpi/gcc/openmpi/lib')
    MPI_BINDIR = _get('MPI_BINDIR', '/usr/lib64/mpi/gcc/openmpi/bin')

    # CUDA
    CUDA_LIB = _get(CUDA_LIB_VAR, '/usr/local/cuda/lib64')
    CUDA_BIN = _get('CUDA_BIN', '/usr/local/cuda/bin')
    MAX_PREVIEW_FILE_SIZE = float(_get("MAX_PREVIEW_FILE_SIZE", DEFAULT_MAX_PREVIEW_FILE_SIZE))

    # OLD CHIMERA variable
    CHIMERA_OLD_BINARY_PATH = _get("CHIMERA_OLD_BINARY_PATH",'')


class Domain(pyworkflow.plugin.Domain):
    _name = __name__
    _objectClass = EMObject
    _protocolClass = Protocol
    _viewerClass = Viewer
    _wizardClass = Wizard
    _baseClasses = getSubclasses(EMObject, globals())


class Plugin(pyworkflow.plugin.Plugin):
    _url = URL

    @classmethod
    def _defineEmVar(cls, varName, defaultValue):
        """ Shortcut method to define variables prepending EM_ROOT if variable is not absolute"""

        # Get the value, either whatever is in the environment or a join of EM_ROOT + defaultValue
        value = os.environ.get(varName, os.path.join(Config.EM_ROOT, defaultValue))

        # CASE-1 : Users might have used ~ and that has to be expanded
        value = os.path.expanduser(value)

        # CASE-2 :Old configs 2.0 will likely have:
        # EM_ROOT = software/em
        # CHIMERA_HOME = %(EM_ROOT)s/chimera-13.0.1
        #  ...
        #
        # this end up in the environment resolved like: software/em/chimera-13.0.1
        # whereas as default values will come as absolute:
        # /<scipion_home>/software/em/chimera-13.0.1
        # In any case we join it (absolute paths will not join)
        value = os.path.join(pw.Config.SCIPION_HOME, value)

        cls._addVar(varName, value)

    @classmethod
    def getMaxitHome(cls):
        return cls.getVar(MAXIT_HOME)

    @classmethod
    def getMaxitBin(cls):
        return os.path.join(cls.getMaxitHome(), 'bin', MAXIT)

    @classmethod
    def _defineVariables(cls):
        # Avoid defining variables from children that does not define variables.
        if cls == Plugin:
            cls._defineVar(EM_ROOT_VAR, pwem.Config.EM_ROOT)
            cls._defineEmVar(MAXIT_HOME, 'maxit-10.1')

            # Take this initialization event to define own datasets
            defineDatasets()
            # Register filehandlers too
            cls._registerFileHandlers()

    @classmethod
    def defineBinaries(cls, env):
        cls.defineBinariesMaxit(False, env)

    @classmethod
    def defineBinariesMaxit(cls, default, env):
        # If not defined already (several plugins needs this and call this but has to be added once
        if not env.hasPackage(MAXIT):
            MAXIT_URL = 'https://sw-tools.rcsb.org/apps/MAXIT/maxit-v10.100-prod-src.tar.gz'
            MAXIT_TAR = 'maxit-v10.100-prod-src.tar.gz'
            maxit_commands = [('make -j 1 binary ', ['bin/maxit'])]
            env.addPackage(MAXIT, version='10.1',
                           tar=MAXIT_TAR,
                           url=MAXIT_URL,
                           neededProgs=['gcc', 'flex', 'make', 'bison', 'tcsh'],
                           commands=maxit_commands,
                           default=default)  # scipion installb maxit
            # requirements bison, flex, gcc

        maxit10 = env.getTarget(env._getExtName(MAXIT, "10.1"))
        maxit10.setDefault(default)


    @classmethod
    def guessCudaVersion(cls, variable, default="10.1"):
        """ Guesses the cuda version from a variable

        :param variable: Name of a variable defined by a plugin
        :param default: Default value in case cuda is not found
        :return Version after parsing the variable or default"""

        return cls.getVersionFromVariable(variable, default=default, pattern="cuda")

    @classmethod
    def getVersionFromVariable(cls, variable, default=NO_VERSION_FOUND_STR, pattern=None):
        """ Returns a Version instance from parsing the content of
         a variable

         :param variable: Variable name to use and parse
         :param default: Optional. Default value to return in case version is not found
         :param pattern: Optional. If passed, the pattern will be used to locate the folder
         :return Version after parsing the variable"""

        value = cls.getVar(variable)
        return cls.getVersionFromPath(value, default= default, pattern=pattern)


    @classmethod
    def getVersionFromPath(cls, path, separator="-", default=NO_VERSION_FOUND_STR, pattern=None):
        """ Resolves path to the realpath (links) and returns the last part
         of the path after the separator as a Version object.
         If separator is not present returns None

         :param path: string containing a path with a version as part of its name
         :param separator: Optional. defaults to "-".
         :param default: Optional. Default version ("X.Y" string) to return if no version is found
         :param pattern: Optional. If passed, the pattern will be used to locate the folder
         :return Version object filled with version found or 0.0

         """

        path = os.path.realpath(path)
        if pattern is not None:
            path = findFolderWithPattern(path, pattern)
            if path is None:
                return parse_version(default)

        parts = path.split(separator)
        if len(parts)>=2:

            # Version should be the last bit
            versionStr = parts[-1]
            return parse_version(versionStr)
        else:
            return parse_version(default)
    @classmethod
    def _registerFileHandlers(cls):
        # register file handlers to preview info in the Filebrowser....
        from pyworkflow.gui.browser import FileTreeProvider, STANDARD_IMAGE_EXTENSIONS
        from .viewers.filehandlers import MdFileHandler, ParticleFileHandler, VolFileHandler, StackHandler, ChimeraHandler

        register = FileTreeProvider.registerFileHandler
        register(MdFileHandler(), '.xmd', '.star', '.pos', '.ctfparam', '.doc')
        register(ParticleFileHandler(),
                 '.xmp', '.tif', '.tiff', '.spi', '.mrc', '.map', '.raw',
                 '.inf', '.dm3', '.em', '.pif', '.psd', '.spe', '.ser', '.img',
                 '.hed', *STANDARD_IMAGE_EXTENSIONS)
        register(VolFileHandler(), '.vol', '.hdf', '.rec')
        register(StackHandler(), '.stk', '.mrcs', '.st', '.pif', '.dm4', '.ali')
        register(ChimeraHandler(), '.bild', '.mrc', '.pdb', '.vol', '.hdf', '.cif', '.mmcif')


def findFolderWithPattern(path, pattern):
    """
    Returns the first folder found from the right that matches the pattern

    :param path: path to do the search on
    :param pattern: pattern to look for. No regex for now, just a contains

    :return the folder matching the pattern or None

    """

    previous, last = os.path.split(path)

    # If end reached: windows Paths and linux paths seems to end in a different way
    if previous == "" or last=="":
        return  None
    elif pattern in last:
        return last
    else:
        return findFolderWithPattern(previous, pattern)

