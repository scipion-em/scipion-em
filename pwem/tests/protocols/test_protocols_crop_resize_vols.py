# **************************************************************************
# *
# * Authors:    Scipion (scipion@cnb.csic.es)
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
from pwem.protocols import ProtImportVolumes, exists
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwem.protocols.protocol_crop_resize_vols import ProtCropResizeVols, OUTPUT_VOLUME



class TestCropResizeBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.map3D = cls.dataset.getFile('betagal')
        cls.half1 = cls.dataset.getFile('betagal_half1')
        cls.half2 = cls.dataset.getFile('betagal_half2')


    @classmethod
    def runImportVolumes(cls, pattern, samplingRate, setHalfMaps=False,
                         half1map='', half2map=''):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate,
                                         setHalfMaps=setHalfMaps,
                                         half1map=half1map,
                                         half2map=half2map
                                         )
        cls.launchProtocol(cls.protImport)
        return cls.protImport
    
    @classmethod
    def runCropResize(cls, inVol, **kwargs):
        protCropResize = cls.newProtocol(ProtCropResizeVols,
                                          inVolume=inVol,
                                           **kwargs)
        cls.launchProtocol(protCropResize)
        return getattr(protCropResize, OUTPUT_VOLUME, None)
    

class TestCropResize(TestCropResizeBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.sampling = 3.54
        TestCropResizeBase.setData()
        cls.protImportVol = cls.runImportVolumes(cls.map3D, cls.sampling)
        cls.protImportWithHalves = cls.runImportVolumes(cls.map3D, cls.sampling, 
                        setHalfMaps=True, half1map=cls.half1, half2map=cls.half2)
    

    def test_resize_sampling(self):
        targetFactors = [2.0, 0.5]
        for factor in targetFactors:
            targetSR = factor * self.sampling
            objLabel=f'resize to sampling {targetSR}',
            vol = self.runCropResize(self.protImportVol.outputVolume,
                                     objLabel=objLabel,
                                     doResize=True,
                                     resizeOption=ProtCropResizeVols.RESIZE_SAMPLINGRATE,
                                     resizeSamplingRate=targetSR,
                                     doCropPad=False)

            self.assertTrue(vol.getDim()[0]==self.protImportVol.outputVolume.getDim()[0]*(1/factor), 
                            f'resizing has failed with factor {factor} or sampling {targetSR}')
 

    def test_resize_dims(self):
        targetSize = [50, 200]
        for boxsize in targetSize:
            objLabel=f'resize to {boxsize}',
            vol = self.runCropResize(self.protImportVol.outputVolume,
                                     objLabel=objLabel,
                                     doResize=True,
                                     resizeOption=ProtCropResizeVols.RESIZE_DIMENSIONS,
                                     resizeDim=boxsize,
                                     doCropPad=False)
            self.assertTrue(vol.getDim()[0]==boxsize, 
                            f'resizing has failed setting a boxsize of {boxsize}')
    
    def test_cropPad(self):
        targetSize = [50, 200]
        for boxsize in targetSize:
            objLabel=f'crop to {boxsize}',
            vol = self.runCropResize(self.protImportVol.outputVolume,
                                     objLabel=objLabel,
                                     doResize=False,
                                     doCropPad=True,
                                     cropPadDim=boxsize
                                     )
            self.assertTrue(vol.getDim()[0]==boxsize, 
                            f'crop/padding has failed setting a boxsize of {boxsize}')
    
    def test_resize_cropPad(self):
        targetSize = [50, 200]
        targetFactors = [2.0, 0.5]
        for factor in targetFactors:
            targetSR = factor * self.sampling
            for boxsize in targetSize:
                objLabel=f'resize sampling {targetSR} crop {boxsize}'
                vol = self.runCropResize(self.protImportWithHalves.outputVolume,
                                        objLabel=objLabel,
                                        doResize=True,
                                        resizeOption=ProtCropResizeVols.RESIZE_SAMPLINGRATE,
                                        resizeSamplingRate=targetSR,
                                        doCropPad=True,
                                        cropPadDim=boxsize
                                        )
                self.assertTrue(vol.getDim()[0]==boxsize, 
                                f'protocol failed with resize {targetSR} and cropPadding of {boxsize}')
                
                fnHalf1, fnHalf2 = vol.getHalfMaps().split(',')
                self.assertTrue(exists(fnHalf1), 'The half map 1 were not cropped/resized')
                self.assertTrue(exists(fnHalf2), 'The half map 2 were not cropped/resized')
    