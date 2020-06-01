# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
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
from os.path import join, dirname, basename
import pyworkflow.utils as pwutils
import pwem


def loadSetFromDb(dbName, dbPrefix=''):
    from pwem import Domain
    from pyworkflow.mapper.sqlite import SqliteFlatDb
    db = SqliteFlatDb(dbName=dbName, tablePrefix=dbPrefix)
    setClassName = db.getProperty('self')  # get the set class name
    setObj = Domain.getObjects()[setClassName](filename=dbName, prefix=dbPrefix)
    return setObj


def runProgram(program, params):
    """ Runs a em program setting its environment matching a prefix"""
    env = None

    # Allow passing absolute paths
    programName = basename(program)

    from pwem import Domain
    if programName.startswith('xmipp'):
        xmipp3 = Domain.getPlugin('xmipp3').Plugin
        env = xmipp3.getEnviron()
    if programName.startswith('relion'):
        relion = Domain.getPlugin("relion").Plugin
        env = relion.getEnviron()
        print(env)
    elif (programName.startswith('e2') or
          programName.startswith('sx')):
        eman2 = Domain.importFromPlugin('eman2', 'Plugin')
        env = eman2.getEnviron()
    elif programName.startswith('b'):
        bsoft = Domain.importFromPlugin('bsoft', 'Plugin')
        env = bsoft.getEnviron()

    pwutils.runJob(None, program, params, env=env)


def getPWEMPath(*paths):
    return join(dirname(pwem.__file__), *paths)


def getTemplatePath(*paths):
    return join(getPWEMPath('templates'), *paths)


def getCmdPath(*paths):
    return join(getPWEMPath('cmd'), *paths)


def convertPixToLength(samplingRate, length):
    return samplingRate * length
