# ***************************************************************************
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# ***************************************************************************/
import os

import pyworkflow.tests as pwtests

import pwem.protocols as emprot
import pwem.constants as emcts
import pwem.convert as emconv
from pwem.convert import Ccp4Header
from pyworkflow.object import Pointer


class TestImportBase(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.dsRelion = pwtests.DataSet.getDataSet('relion_tutorial')


class TestImportVolumes2(TestImportBase):

    def test_import_volume(self):
        # Import three volumes (mrc stack) and set origin in userDefined position
        args = {'filesPath': self.dsRelion.getFile('import/case2/relion_volumes.mrc'),
                'filesPattern': '',
                'setHalfMaps': False,
                'setOrigCoord': True,
                'samplingRate': 2.1,
                'x': 16.8,
                'y': 33.6,
                'z': 50.4
                }
        prot2 = self.newProtocol(emprot.ProtImportVolumes, **args)
        prot2.setObjLabel('import 3 vols from mrc stack,\n user origin,'
                          '\n chimera no axis')
        self.launchProtocol(prot2)
        # Check the number of output volumes and dimensions
        self.assertEqual(3, prot2.outputVolumes.getSize())
        self.assertEqual(60, prot2.outputVolumes.getDim()[0])

        # TODO: Daniel - add a call to the corresponding methods of the test centralized
        # layer you are writing.

