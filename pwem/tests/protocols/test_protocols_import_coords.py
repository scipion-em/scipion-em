# ***************************************************************************
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es)
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

import pyworkflow.tests as pwtests

import pwem.protocols as emprot


class TestImportBase(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.dsXmipp = pwtests.DataSet.getDataSet('xmipp_tutorial')
        # cls.dsRelion = DataSet.getDataSet('relion_tutorial')
        cls.dsGroel = pwtests.DataSet.getDataSet('groel')
        
    def checkOutput(self, prot, outputName, conditions=[]):
        """ Check that an ouput was generated and
        the condition is valid. 
        """
        o = getattr(prot, outputName, None)
        locals()[outputName] = o 
        self.assertIsNotNone(o, "Output: %s is None" % outputName)
        for cond in conditions:
            self.assertTrue(eval(cond), 'Condition failed: ' + cond)
        
    
class TestImportCoordinates(TestImportBase):

    def testImportCoordinates(self):
        # First, import a set of micrographs
        protImport = self.newProtocol(emprot.ProtImportMicrographs,
                                      filesPath=self.dsXmipp.getFile('allMics'),
                                      samplingRate=1.237, voltage=300)
        protImport.setObjLabel('import micrographs from xmipp tutorial ')
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs.getFileName(),
                             "There was a problem with the import")

        prot1 = self.newProtocol(emprot.ProtImportCoordinates,
                                 importFrom=emprot.ProtImportCoordinates.IMPORT_FROM_XMIPP,
                                 filesPath=self.dsXmipp.getFile('pickingXmipp'),
                                 filesPattern='*.pos', boxSize=550,
                                 scale=3.,
                                 invertX=False,
                                 invertY=False
                                 )
        prot1.inputMicrographs.set(protImport.outputMicrographs)
        prot1.setObjLabel('import coords from xmipp ')
        self.launchProtocol(prot1)
        
        # Make sure that all 264 coordinates where correctly imported
        self.assertTrue(prot1.outputCoordinates.getSize() == 264)

        filepath = prot1.outputCoordinates.getFileName()
        prot2 = self.newProtocol(emprot.ProtImportCoordinates,
                                 importFrom=emprot.ProtImportCoordinates.IMPORT_FROM_SCIPION,
                                 filesPath=filepath, 
                                 filesPattern='',
                                 boxSize=110,
                                 )
        prot2.inputMicrographs.set(protImport.outputMicrographs)
        prot2.setObjLabel('import coords from scipion ')
        self.launchProtocol(prot2)
        # Make sure that all 264 coordinates where correctly imported
        self.assertTrue(prot2.outputCoordinates.getSize() == 264)

        # prot2 = self.newProtocol(ProtImportCoordinates,
        #                          importFrom=ProtImportCoordinates.IMPORT_FROM_RELION,
        #                          filesPath=self.dsXmipp.getFile('boxingDir'),#no dataset with picking
        #                          filesPattern='info/*_info.json',
        #                          boxSize=110,
        #                          scale=2,
        #                          invertX=False,
        #                          invertY=False
        #                          )
        # prot2.inputMicrographs.set(protImport.outputMicrographs)
        # prot2.setObjLabel('import coords from relion ')

        prot3 = self.newProtocol(emprot.ProtImportCoordinates,
                                 importFrom=emprot.ProtImportCoordinates.IMPORT_FROM_EMAN,
                                 filesPath=self.dsXmipp.getFile('boxingDir'),
                                 filesPattern='*_info.json', boxSize=550,
                                 scale=5.,
                                 invertX=False,
                                 invertY=False)
        prot3.inputMicrographs.set(protImport.outputMicrographs)
        prot3.setObjLabel('import coords from eman ')

        self.launchProtocol(prot3)

        # First, import a set of micrographs
        protImportGroel = self.newProtocol(emprot.ProtImportMicrographs,
                                           filesPath=self.dsGroel.getFile('mic1'),
                                           samplingRate=1)
        protImportGroel.setObjLabel('import micrographs from groel')
        self.launchProtocol(protImportGroel)
        self.assertIsNotNone(protImportGroel.outputMicrographs.getFileName(),
                             "There was a problem with the import")

        protPickGroel = self.newProtocol(emprot.ProtImportCoordinates,
                                         importFrom=emprot.ProtImportCoordinates.IMPORT_FROM_DOGPICKER,
                                         filesPath=self.dsGroel.getFile('pickingDogpicker'),
                                         filesPattern='*.txt', boxSize=10,
                                         threshold=0.7)
        protPickGroel.inputMicrographs.set(protImportGroel.outputMicrographs)
        protPickGroel.setObjLabel('import coords from dogpicker ')

        self.launchProtocol(protPickGroel)


class TestImportCoordinatesPairs(TestImportBase):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.dsRct = pwtests.DataSet.getDataSet('rct')
        cls.micsFn = cls.dsRct.getFile('positions/input_micrographs.xmd')
        cls.patternU1 = cls.dsRct.getFile('positions/F_rct_u_*.pos')
        cls.patternT1 = cls.dsRct.getFile('positions/F_rct_t_*.pos')
        cls.micsUFn1 = cls.dsRct.getFile('untilted')
        cls.micsTFn1 = cls.dsRct.getFile('tilted')

        cls.dsEman = pwtests.DataSet.getDataSet('eman')
        cls.micsUFn2 = cls.dsEman.getFile('micU')
        cls.micsTFn2 = cls.dsEman.getFile('micT')
        cls.patternU2 = cls.dsEman.getFile("coords/ip3r10252011-0005_0-2_info.json")
        cls.patternT2 = cls.dsEman.getFile("coords/ip3r10252011-0005_10_info.json")

    def testImportCoordinatesPairs(self):
        # First, import a set of micrograph pairs
        protImport1 = self.newProtocol(emprot.ProtImportMicrographsTiltPairs,
                                       patternUntilted=self.micsUFn1,
                                       patternTilted=self.micsTFn1,
                                       samplingRate=2.28, voltage=100,
                                       sphericalAberration=2.9)
        protImport1.setObjLabel('import tilt pair micrographs from rct tutorial ')
        self.launchProtocol(protImport1)
        self.assertIsNotNone(protImport1.outputMicrographsTiltPair.getFileName(),
                             "There was a problem with the import")

        protImport2 = self.newProtocol(emprot.ProtImportMicrographsTiltPairs,
                                       patternUntilted=self.micsUFn2,
                                       patternTilted=self.micsTFn2,
                                       samplingRate=2.8, voltage=200,
                                       sphericalAberration=2.0)
        protImport2.setObjLabel('import tilt pair micrographs from eman2 tutorial ')
        self.launchProtocol(protImport2)
        self.assertIsNotNone(protImport2.outputMicrographsTiltPair.getFileName(),
                             "There was a problem with the import")

        prot1 = self.newProtocol(emprot.ProtImportCoordinatesPairs,
                                 importFrom=emprot.ProtImportCoordinatesPairs.IMPORT_FROM_XMIPP,
                                 xmippMdFn=self.micsFn,
                                 patternUntilted=self.patternU1,
                                 patternTilted=self.patternT1,
                                 boxSize=100)
        prot1.inputMicrographsTiltedPair.set(protImport1.outputMicrographsTiltPair)
        prot1.setObjLabel('import coord pairs from xmipp ')
        self.launchProtocol(prot1)

        # Make sure that all 1901 coordinates where correctly imported
        self.assertTrue(prot1.outputCoordinatesTiltPair.getSize() == 1901)

        prot2 = self.newProtocol(emprot.ProtImportCoordinatesPairs,
                                 importFrom=emprot.ProtImportCoordinatesPairs.IMPORT_FROM_EMAN,
                                 patternUntilted=self.patternU2,
                                 patternTilted=self.patternT2,
                                 boxSize=256)
        prot2.inputMicrographsTiltedPair.set(protImport2.outputMicrographsTiltPair)
        prot2.setObjLabel('import coord pairs from eman2 ')
        self.launchProtocol(prot2)

        # Make sure that all 104 coordinates where correctly imported
        self.assertTrue(prot2.outputCoordinatesTiltPair.getSize() == 104)
