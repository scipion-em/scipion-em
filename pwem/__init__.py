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

import pyworkflow as pw
from pyworkflow.protocol import Protocol
from pyworkflow.viewer import Viewer
from pyworkflow.wizard import Wizard
import pyworkflow.plugin

from .constants import *
from .objects import EMObject
from .tests import defineDatasets
from .utils import *

__version__ = '3.0.14'
_logo = "scipion_icon.gif"
_references = ["delaRosaTrevin201693"]


class Config(pw.Config):
    _get = pw.Config._get
    _join = pw.Config._join

    EM_ROOT = _join(_get(EM_ROOT_VAR, _join(pw.Config.SCIPION_SOFTWARE, 'em')))

    # Default XMIPP_HOME: needed here for ShowJ viewers
    XMIPP_HOME = _join(_get('XMIPP_HOME', os.path.join(EM_ROOT, 'xmipp')))

    # Get java home, we might need to provide correct default value
    JAVA_HOME = _get('JAVA_HOME', '')
    JAVA_MAX_MEMORY = _get('JAVA_MAX_MEMORY', '4')

    # MPI
    MPI_LIBDIR = _get('MPI_LIBDIR', '/usr/lib64/mpi/gcc/openmpi/lib')
    MPI_BINDIR = _get('MPI_BINDIR', '/usr/lib64/mpi/gcc/openmpi/bin')

    # CUDA
    CUDA_LIB = _get('CUDA_LIB', '/usr/local/cuda/lib64')
    CUDA_BIN = _get('CUDA_BIN', '/usr/local/cuda/bin')
    MAX_PREVIEW_FILE_SIZE = float(_get("MAX_PREVIEW_FILE_SIZE", DEFAULT_MAX_PREVIEW_FILE_SIZE))


class Domain(pyworkflow.plugin.Domain):
    _name = __name__
    _objectClass = EMObject
    _protocolClass = Protocol
    _viewerClass = Viewer
    _wizardClass = Wizard
    _baseClasses = globals()


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
        cls._defineVar(EM_ROOT_VAR, pwem.Config.EM_ROOT)
        cls._defineEmVar(MAXIT_HOME, 'maxit-10.1')

        # Take this initialization event to define own datasets
        defineDatasets()

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
                           commands=maxit_commands,
                           default=default)  # scipion installb maxit
            # requirements bison, flex, gcc

        maxit10 = env.getTarget(env._getExtName(MAXIT, "10.1"))
        maxit10.setDefault(default)
