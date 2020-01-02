# ***************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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
# ***************************************************************************/

import os


import pyworkflow.tests as pwtests
import pyworkflow.protocol as pwprot

import pwem.protocols as emprot


# Load the number of movies for the simulation, by default equal 5, but
# can be modified in the environement
MICS = os.environ.get('SCIPION_TEST_MICS', 3)


class TestCtfStream(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)

    def _updateProtocol(self, prot):
        prot2 = pwprot.getProtocolFromDb(prot.getProject().path,
                                         prot.getDbPath(),
                                         prot.getObjId())
        # Close DB connections
        prot2.getProject().closeMapper()
        prot2.closeMappers()
        return prot2

    def test_pattern(self):
        """ Import several Particles from a given pattern.
        """
        kwargs = {'xDim': 4096,
                  'yDim': 4096,
                  'nDim': MICS,
                  'samplingRate': 1.25,
                  'creationInterval': 5,
                  'delay': 0,
                  'setof': emprot.SET_OF_RANDOM_MICROGRAPHS  # SetOfMicrographs
                  }

        # put some stress on the system
        protStream = self.newProtocol(emprot.ProtCreateStreamData, **kwargs)
        protStream.setObjLabel('create Stream Mic')
        self.proj.launchProtocol(protStream, wait=False)

        self._waitOutput(protStream, 'outputMicrographs')

        # then introduce monitor, checking all the time ctf and saving to
        # database
        kwargs = {
            'useCtffind4': True,
            'ctfDownFactor': 2,
            'highRes': 0.4,
            'lowRes': 0.05,
            'numberOfThreads': 4
        }
        from pwem import Domain

        ProtCTFFind = Domain.importFromPlugin('cistem.protocols',
                                              'CistemProtCTFFind', doRaise=True)
        protCTF = self.newProtocol(ProtCTFFind, **kwargs)
        protCTF.inputMicrographs.set(protStream.outputMicrographs)
        self.proj.launchProtocol(protCTF, wait=False)

        self._waitOutput(protCTF, 'outputCTF')

        kwargs = {'samplingInterval': 10,
                  'interval': 300,
                  'maxDefocus': 40000,
                  'minDefocus': 1000,
                  'astigmatism': 0.2,
                  'monitorTime': 5
                  }

        protMonitor = self.newProtocol(emprot.ProtMonitorCTF, **kwargs)
        protMonitor.inputProtocol.set(protCTF)
        self.launchProtocol(protMonitor)

        baseFn = protMonitor._getPath(emprot.CTF_LOG_SQLITE)
        self.assertTrue(os.path.isfile(baseFn))
