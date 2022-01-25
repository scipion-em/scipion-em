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
        cls.dsModBuild = pwtests.DataSet.getDataSet('model_building_tutorial')


class TestImportAtomStruct(TestImportBase):

    def test_import_atomStruct(self):
        """ 1) Import single atomstruct form file
        """
        args = {'inputPdbData': emprot.ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBx_mmCIF/1ake_mut1.pdb'),
                }

        prot1 = self.newProtocol(emprot.ProtImportPdb, **args)
        prot1.setObjLabel('import atomstruct from file')
        self.launchProtocol(prot1)
        self.assertIsNotNone(prot1.outputPdb)


        # """ 2) Import single atomstruct from pdb id
        # """
        args = {'inputPdbData': emprot.ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': '1ake',
                }

        # Id's should be set increasing from 1 if ### is not in the
        # pattern
        prot2 = self.newProtocol(emprot.ProtImportPdb, **args)
        prot2.setObjLabel('import atomstruct from pdb id')
        self.launchProtocol(prot2)
        self.assertIsNotNone(prot2.outputPdb)


    def test_import_setOfAtomStruct(self):
        # """ 1) Import set of atomstruct from files
        # """
        args = {'inputPdbData': emprot.ProtImportSetOfAtomStructs.IMPORT_FROM_FILES,
                'filesPath': self.dsModBuild.getFile('PDBx_mmCIF'),
                'filesPattern': self.dsModBuild.getFile('PDBx_mmCIF/1ake_*.pdb')
                }
        prot1 = self.newProtocol(emprot.ProtImportSetOfAtomStructs, **args)
        prot1.setObjLabel('import 1akes structures from files\n')
        self.launchProtocol(prot1)
        self.assertIsNotNone(getattr(prot1, prot1._OUTNAME))

        # """ 2) Import set of atomstruct form pdb id
        # """
        args = {'inputPdbData': emprot.ProtImportSetOfAtomStructs.IMPORT_FROM_ID,
                'pdbIds': '1ake,5ni1',
                }
        prot2 = self.newProtocol(emprot.ProtImportSetOfAtomStructs, **args)
        prot2.setObjLabel('import 1ake and 5ni1 structures from pdb ids\n')
        self.launchProtocol(prot2)
        self.assertIsNotNone(getattr(prot2, prot2._OUTNAME))