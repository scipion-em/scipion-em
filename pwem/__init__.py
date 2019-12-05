# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
This modules contains classes related with EM
"""

from pyworkflow.protocol import Protocol
from pyworkflow.viewer import Viewer
from pyworkflow.wizard import Wizard
import pyworkflow.plugin
import os

from .constants import *
from pwem.objects import EMObject
from .utils import *


class Domain(pyworkflow.plugin.Domain):
    _name = __name__
    _objectClass = EMObject
    _protocolClass = Protocol
    _viewerClass = Viewer
    _wizardClass = Wizard
    _baseClasses = globals()


class Plugin(pyworkflow.plugin.Plugin):

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
        MAXIT_URL = 'https://sw-tools.rcsb.org/apps/MAXIT/maxit-v10.100-prod-src.tar.gz'
        MAXIT_TAR = 'maxit-v10.100-prod-src.tar.gz'
        maxit_commands = [('make -j 1 binary ' , ['bin/maxit'])]
        env.addPackage(MAXIT, version='10.1',
                       tar=MAXIT_TAR,
                       url=MAXIT_URL,
                       commands=maxit_commands,
                       default=False)  # scipion installb maxit
                                       # requirements bison, flex, gcc


Domain.registerPlugin(__name__)

