# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
"""
This modules contains classes related with EM
"""

import os
import sys

import pyworkflow as pw
from pyworkflow.protocol import Protocol
from pyworkflow.viewer import Viewer
from pyworkflow.wizard import Wizard
import pyworkflow.plugin

from .constants import *
from .objects import EMObject
from .utils import *

_logo = "scipion_icon.gif"
_references = ["delaRosaTrevin201693"]

class Config:
    __get = os.environ.get  # shortcut

    # This will prepend SCIPION_HOME in case the variable is not absolute
    __prefixHome = lambda path: path if os.path.isabs(path) else os.path.join(pw.Config.SCIPION_HOME, path)

    EM_ROOT = __prefixHome(__get('EM_ROOT', os.path.join(pw.Config.SCIPION_SOFTWARE, 'em')))
    # Default XMIPP_HOME: needed here for ShowJ viewers
    XMIPP_HOME = __prefixHome(__get('XMIPP_HOME', os.path.join(EM_ROOT, 'xmipp')))

    # Needed by Chimera viewer.
    CHIMERA_HOME = __get('CHIMERA_HOME', os.path.join(EM_ROOT,'chimera-1.13.1'))

    # Get java home, we might need to provide correct default value
    JAVA_HOME = __get('JAVA_HOME', '')
    JAVA_MAX_MEMORY = __get('JAVA_MAX_MEMORY', '2')

    # CUDA
    CUDA_LIB = __get('CUDA_LIB', '/usr/local/cuda/lib64')
    CUDA_BIN = __get('CUDA_BIN', '/usr/local/cuda/bin')

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
        """ Shortcut method to define variables by prepending EM_ROOT
        to the default value.

        If the varName is defined, we will make the value absolute adding SCIPION_HOME
        """
        prefixedDefault = os.path.join(Config.EM_ROOT, defaultValue)
        varValue = os.environ.get(varName)

        # If en value is not absolute, then we prepend SCIPION_HOME
        # This case addresses users reusing old config files where variable comes like:
        # CHIMERA_HEADLESS_HOME="software/em/chimera_headless"
        # CHIMERA_HOME="software/em/chimera-1.10.1"
        varValue = prefixedDefault if varValue is None else os.path.join(pw.Config.SCIPION_HOME, varValue)

        cls._addVar( varName, varValue)


    @classmethod
    def getMaxitHome(cls):
        return cls.getVar(MAXIT_HOME)

    @classmethod
    def getMaxitBin(cls):
        return os.path.join(cls.getMaxitHome(), 'bin', MAXIT)

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(MAXIT_HOME, 'maxit-10.1')

    @classmethod
    def defineBinaries(cls, env):
        cls.defineBinariesMaxit(False, env)

    @classmethod
    def defineBinariesMaxit(cls, default, env):
        MAXIT_URL = 'https://sw-tools.rcsb.org/apps/MAXIT/maxit-v10.100-prod-src.tar.gz'
        MAXIT_TAR = 'maxit-v10.100-prod-src.tar.gz'
        maxit_commands = [('make -j 1 binary ', ['bin/maxit'])]
        env.addPackage(MAXIT, version='10.1',
                       tar=MAXIT_TAR,
                       url=MAXIT_URL,
                       commands=maxit_commands,
                       default=default)  # scipion installb maxit
        # requirements bison, flex, gcc

