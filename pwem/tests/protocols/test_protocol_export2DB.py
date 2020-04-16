# **************************************************************************
# *
# * Authors:    Amaya Jimenez (ajimenez@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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

import pyworkflow.tests as pwtest

import pwem.protocols as emprot
import pwem.constants as emcts
from .test_protocols_import_volumes import createFeatVolume

class TestExport2DataBases(pwtest.BaseTest):
    @classmethod
    def runImportVolumes(cls, pattern, samplingRate, label):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(emprot.ProtImportVolumes,
                                         objLabel=label,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                         )
        cls.launchProtocol(cls.protImport)
        return cls.protImport

    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = pwtest.DataSet.getDataSet(dataProject)
        cls.half1 = cls.dataset.getFile('betagal_half1')
        cls.half2 = cls.dataset.getFile('betagal_half2')
        cls.dsModBuild = pwtest.DataSet.getDataSet('model_building_tutorial')

    @classmethod
    def setUpClass(cls):
        pwtest.setupTestProject(cls)
        cls.setData()
        cls.protImportHalf1 = cls.runImportVolumes(cls.half1, 3.54,
                                                   'import half1')
        cls.protImportHalf2 = cls.runImportVolumes(cls.half2, 3.54,
                                                   'import half2')

    @classmethod
    def _importAtomStructCIF(self):
        args = {'inputPdbData': emprot.ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/5ni1.cif'),
                }
        protImportPDB = self.newProtocol(emprot.ProtImportPdb, **args)
        protImportPDB.setObjLabel('import atom struct\nmmCIF\n5ni1.cif')
        self.launchProtocol(protImportPDB)
        structure_mmCIF = protImportPDB.outputPdb
        return structure_mmCIF

    def test_volume1_AtomStructCIF(self):
        from pwem import Domain
        XmippProtResolution3D = Domain.importFromPlugin('xmipp3.protocols',
                                                        'XmippProtResolution3D',
                                                        doRaise=True)

        prot = self.newProtocol(XmippProtResolution3D)
        prot.inputVolume.set(self.protImportHalf1.outputVolume)
        prot.referenceVolume.set(self.protImportHalf2.outputVolume)
        self.launchProtocol(prot)

        protExp = self.newProtocol(emprot.ProtExportDataBases)
        protExp.setObjLabel("export to DB\ndir1")
        protExp.exportVolume.set(self.protImportHalf1.outputVolume)
        protExp.exportFSC.set(prot.outputFSC)
        protExp.exportAtomStruct.set(self._importAtomStructCIF())

        protExp.filesPath.set(os.getcwd() + "/dir1")
        self.launchProtocol(protExp)

        # Check the files were generated properly.
        # protExp._createFileNamesTemplates()
        dirName = protExp.filesPath.get()
        nameVolume = os.path.join(dirName, protExp.VOLUMENAME)
        nameFsc = os.path.join(dirName, "fsc_%02d.xml" % 1)
        nameAtomStruct = os.path.join(dirName, protExp.COORDINATEFILENAME)
        self.assertTrue(os.path.exists(nameVolume))
        self.assertTrue(os.path.exists(nameFsc))
        self.assertTrue(os.path.exists(nameAtomStruct))

        # Check if the files have the correct data
        orig_list_x, orig_list_y = protExp.exportFSC.get().getData()
        fo = open(nameFsc, "rU")
        saved_x = []
        orig_x = []
        count = 0
        for line in fo:
            if line[0:3] == '<x>':
                saved_x.append(int(float(line[3:-5])*1000))
                orig_x.append(int(orig_list_x[count]*1000))
                count = count+1

        self.assertListEqual(orig_x, saved_x)

    def test_volume1_AtomStructPDB(self):
        from pwem import Domain
        XmippProtResolution3D = Domain.importFromPlugin('xmipp3.protocols',
                                                        'XmippProtResolution3D',
                                                        doRaise=True)
        prot = self.newProtocol(XmippProtResolution3D)
        prot.inputVolume.set(self.protImportHalf1.outputVolume)
        prot.referenceVolume.set(self.protImportHalf2.outputVolume)
        self.launchProtocol(prot)

        protExp = self.newProtocol(emprot.ProtExportDataBases)
        protExp.setObjLabel("export to DB\ndir2")
        protExp.exportVolume.set(self.protImportHalf1.outputVolume)
        protExp.exportFSC.set(prot.outputFSC)

        # run Chimera rigid fit ti get a PDB file
        extraCommands = ""
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 prefix DONOTSAVESESSION_')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': self.protImportHalf1.outputVolume,
                'pdbFileToBeRefined': self._importAtomStructCIF()
                }
        from pwem import Domain

        ChimeraProtRigidFit = Domain.importFromPlugin('chimera.protocols',
                                                      'ChimeraProtRigidFit',
                                                      doRaise=True)
        protChimera = self.newProtocol(ChimeraProtRigidFit, **args)
        protChimera.setObjLabel('chimera fit\n volume and PDB\n save model')
        self.launchProtocol(protChimera)

        protExp.exportAtomStruct.set(protChimera.DONOTSAVESESSION_Atom_struct__2)

        protExp.filesPath.set(os.getcwd() + "/dir2")
        self.launchProtocol(protExp)

        # Check the files were generated properly.
        # protExp._createFileNamesTemplates()
        dirName = protExp.filesPath.get()
        nameVolume = os.path.join(dirName, protExp.VOLUMENAME)
        nameFsc = os.path.join(dirName, "fsc_%02d.xml" % 1)
        nameAtomStruct = os.path.join(dirName, protExp.COORDINATEFILENAME)
        self.assertTrue(os.path.exists(nameVolume))
        self.assertTrue(os.path.exists(nameFsc))
        self.assertTrue(os.path.exists(nameAtomStruct))

        # Check if the files have the correct data
        orig_list_x, orig_list_y = protExp.exportFSC.get().getData()
        fo = open(nameFsc, "rU")
        saved_x = []
        orig_x = []
        count = 0
        for line in fo:
            if line[0:3] == '<x>':
                saved_x.append(int(float(line[3:-5])*1000))
                orig_x.append(int(orig_list_x[count]*1000))
                count = count+1

        self.assertListEqual(orig_x, saved_x)

    def import_volume_halfmaps(self):
        """
        Test to import a full map (Icosahedron) and two maps (half1 and half2)
        to compute the FSC
        """
        volFeatName = '/tmp/Icosahedron_map.txt'
        volMapNamefull = '/tmp/Icosahedron_map_full.mrc'
        volMapNamehalf1 = '/tmp/Icosahedron_map_half1.mrc'
        volMapNamehalf2 = '/tmp/Icosahedron_map_half2.mrc'
        volMaskName = '/tmp/mask.mrc'
        createFeatVolume(volFeatName, volMapNamefull, sym=emcts.SYM_I222r)
        createFeatVolume(volFeatName, volMapNamehalf1, sym=emcts.SYM_I222r)
        createFeatVolume(volFeatName, volMapNamehalf2, sym=emcts.SYM_I222r, factor=1.005)
        createFeatVolume(volFeatName, volMaskName, sym=emcts.SYM_I222r, factor=1.2)
        _samplingRate = 1.0

        # import 3d map
        args = {'filesPath': volMapNamefull,
                'filesPattern': '',
                'setHalfMaps': True,
                'half1map': volMapNamehalf1,
                'half2map': volMapNamehalf2,
                'samplingRate': _samplingRate
                }
        prot5 = self.newProtocol(emprot.ProtImportVolumes, **args)
        prot5.setObjLabel('import phantom icosahedron,\n half1 and half2')
        self.launchProtocol(prot5)
        self.volume = prot5.outputVolume
        self.volume.setOrigin(None)
        # The volume has no origin
        t = self.volume.getOrigin(force=True)
        x, y, z = t.getShifts()
        self.assertEqual(-90.0, x)
        self.assertEqual(-90.0, y)
        self.assertEqual(-90.0, z)
        # import half maps
        args = {'filesPath': volMapNamehalf1,
                'filesPattern': '',
                'samplingRate': _samplingRate
                }
        prot5_1 = self.newProtocol(emprot.ProtImportVolumes, **args)
        prot5_1.setObjLabel('import phantom\nicosahedron,\nhalf1')
        self.launchProtocol(prot5_1)
        self.half1 = prot5_1.outputVolume
        self.half1.setOrigin(None)
            # The volume has no origin
        t = self.volume.getOrigin(force=True)
        x, y, z = t.getShifts()
        self.assertEqual(-90.0, x)
        self.assertEqual(-90.0, y)
        self.assertEqual(-90.0, z)

        args = {'filesPath': volMapNamehalf2,
                'filesPattern': '',
                'samplingRate': _samplingRate
                }
        prot5_2 = self.newProtocol(emprot.ProtImportVolumes, **args)
        prot5_2.setObjLabel('import phantom\nicosahedron,\nhalf2')
        self.launchProtocol(prot5_2)
        self.half2 = prot5_2.outputVolume
        self.half2.setOrigin(None)
        # The volume has no origin
        t = self.volume.getOrigin(force=True)
        x, y, z = t.getShifts()
        self.assertEqual(-90.0, x)
        self.assertEqual(-90.0, y)
        self.assertEqual(-90.0, z)
        # import mask
        args = {'maskPath': volMaskName,
                'samplingRate': _samplingRate
                }
        prot6 = self.newProtocol(emprot.ProtImportMask, **args)
        prot6.setObjLabel('import mask') #  add origin to the mask
        self.launchProtocol(prot6)
        self.mask = prot6.outputMask
        self.mask.setOrigin(None)
        # The volume has no origin
        t = self.mask.getOrigin(force=True)
        x, y, z = t.getShifts()
        self.assertEqual(-90.0, x)
        self.assertEqual(-90.0, y)
        self.assertEqual(-90.0, z)


    def test_halfmaps_mask(self):
        """ If the input is a 3D map with half volumes associated
            Save them to the output directory
            Also create a mask and export it
        """
        # first create the 3 3dmaps and the mask
        self.import_volume_halfmaps()
        # run the export protocol
        from pwem import Domain
        XmippProtResolution3D = Domain.importFromPlugin('xmipp3.protocols',
                                                        'XmippProtResolution3D',
                                                        doRaise=True)

        prot = self.newProtocol(XmippProtResolution3D)
        prot.inputVolume.set(self.half1)
        prot.referenceVolume.set(self.half2)
        self.launchProtocol(prot)

        protExp = self.newProtocol(emprot.ProtExportDataBases)
        protExp.setObjLabel("export to DB\ndir3")
        protExp.exportVolume.set(self.volume)
        protExp.additionalVolumesToExport.set(True)
        protExp.exportAdditionalVolumes.set([self.protImportHalf1.outputVolume,
                                            self.protImportHalf2.outputVolume])
        protExp.exportFSC.set(prot.outputFSC)
        protExp.masksToExport.set(True)
        protExp.exportMasks.set([self.mask])
        protExp.exportAtomStruct.set(self._importAtomStructCIF())

        protExp.filesPath.set(os.getcwd() + "/dir3")
        self.launchProtocol(protExp)

        # Check the files were generated properly.
        # protExp._createFileNamesTemplates()
        dirName = protExp.filesPath.get()
        nameVolume = os.path.join(dirName, protExp.VOLUMENAME)
        nameFsc = os.path.join(dirName, "fsc_%02d.xml" % 1)
        nameAtomStruct = os.path.join(dirName, protExp.COORDINATEFILENAME)
        self.assertTrue(os.path.exists(nameVolume))
        self.assertTrue(os.path.exists(nameFsc))
        self.assertTrue(os.path.exists(nameAtomStruct))
        # assert results