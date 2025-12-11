from typing import Optional, List, Union, Tuple
from os.path import exists

import mrcfile
import numpy as np

from pyworkflow.tests import BaseTest
from pwem.objects import Volume, Image, SetOfVolumes, SetOfImages, Acquisition, Movie, Transform

class TestBaseCentralizedLayer(BaseTest):

    #def checkTomoAcquisition(self,
    #                         testAcqObj: Acquisition = None,
    #                         tiltAnglesTolDeg: float = 0.01,None
    #                         rotAngleTolDeg: float = 0.01) -> None:
    #    self.assertAlmostEqual(testAcq.getVoltage(), currentAcq.getVoltage(), delta=1)
        
    def _checkVolumeOrigin(self,
                           vol: Volume,
                           expectedOriginShifts: Union[List[float], None] = None) -> None:
        if expectedOriginShifts is not None:
            #x, y, z = vol.getTransform().getShifts()
            x, y, z = vol.getOrigin().getShifts()
            for i, j in zip([x, y, z], expectedOriginShifts):
                self.assertAlmostEqual(i, j, delta=0.5, msg="Expected and resulting volume shifts are different")

    def checkTransformMatrix(self,
                             outMatrix: np.ndarray,
                             alignment: bool = False,
                             is2d: bool = False):
        """Checks the shape and coarsely the contents of the transformation matrix provided.

        :param outMatrix: transformation matrix of a subtomogram or coordinate.
        :param alignment: False by default. Used to specify if the expected transformation matrix should be an
        eye matrix (False) or not (True).
        :param is2d: False by default. Used to indicate if the expected transformation matrix should be of size 3 x 3
        (True) or 4 x 4 (False).
        """
        size = 3 if is2d else 4
        transfMatrixShape = (size, size)
        identityMatrix = np.eye(size)
        self.assertIsNotNone(outMatrix)
        if type(outMatrix) is not np.ndarray:
            outMatrix = np.array(outMatrix)
        self.assertIsNotNone(outMatrix)
        self.assertEqual(outMatrix.shape, transfMatrixShape)
        
        if alignment:
            self.assertFalse(np.array_equal(outMatrix, identityMatrix))
        else:
            self.assertTrue(np.array_equal(outMatrix, identityMatrix))
      
    def checkSetGeneralProps(self, inSet, expectedSetSize: int, expectedSRate: float,
                             streamState: int = 2,
                             sRateAngsPixTol: float = 0.01) -> None:
        """
        :param inSet: A set of Scipion Tomo objects.
        :param expectedSetSize: expected set site to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param streamState: expected stream state, being 2 (default) a stream that is closed.
        :param sRateAngsPixTol: tolerance, in angtroms/pixel, of the sampling rate.
        """
        if expectedSetSize > 0:
            self.assertSetSize(inSet, expectedSetSize)
        self.assertAlmostEqual(inSet.getSamplingRate(), expectedSRate, delta=sRateAngsPixTol)
        self.assertEqual(inSet.getStreamState(), streamState)

        print("Checking properties are persisted")
        self.assertTrue(inSet.hasProperty("self"),f"Set {inSet.getFileName()} does not have 'self' in properties table. Probably properties are not persisted.")

    def checkHeaderSRate(self,
                         inObj: Union[SetOfVolumes, Volume, SetOfImages, Image],
                         expectedSRate: float,
                         sRateAngsPixTol: float = 0.01) -> None:
        """Checks the headers of a voume or image.

        :param TODO.
        """
        volumeFile = inObj.getFileName()
        with mrcfile.open(volumeFile, permissive=True, header_only=True) as mrc:
            vs = mrc.voxel_size
            vs = [float(vs.x), float(vs.y), float(vs.z)] if type(inObj) in [Volume, SetOfVolumes] else [float(vs.x), float(vs.y)]
            for voxelSize in vs:
                self.assertAlmostEqual(voxelSize,
                                       expectedSRate,
                                       delta=sRateAngsPixTol,
                                       msg=f"tolerance {sRateAngsPixTol}\n"
                                           f"File = {volumeFile}\n"
                                           f"Expected sampling rate = {expectedSRate} angst / pix\n"
                                           f"Sampling rate/s in header  = {vs}")

    def checkVolume(self,
                    vol: Volume,
                    expectedSRate: float = -1.,
                    expectedBoxSize: int = -1,
                    hasCtf: bool = False,
                    hasHalves: bool = True,
                    sRateAngsPixTol: float = 0.01,
                    expectedOriginShifts: Union[List[float], None] = None) -> None:
        """Checks the main properties of a voume, which can be TODO.

        :param volume: Volume object to check.
        :param expectedBoxSize: expected box size, in pixels, to check.
        :param expectedSRate: expected sampling rate, in Å/pix, to check.
        :param hasHalves: True by default. Used to indicate if the average is expected to have halves associated.
        :param sRateAngsPixTol: tolerance, in angstroms/pixel, of the sampling rate.
        """
        print(type(vol.getTransform()))
        print(type(vol.getOrigin()))
        testBoxSize = (expectedBoxSize, expectedBoxSize, expectedBoxSize)
        self.assertTrue(exists(vol.getFileName()), f"Volume {vol.getFileName()} does not exists")
        self.assertTrue(vol.getFileName().endswith(".mrc"), f"Volume {vol.getFileName()} is not an MRC file.")
        self.assertTrue(abs(vol.getSamplingRate() - expectedSRate) <= 1e-2, "Volume doesn't have the expected rate.")
        self.assertEqual(vol.getDim(), testBoxSize, "Volume doesn't have the expected dimensions.")
        # Check the sampling rate value in the header
        self.checkHeaderSRate(vol, expectedSRate=expectedSRate, sRateAngsPixTol=sRateAngsPixTol)
        # Check the halves
        if hasHalves:
            self.assertTrue(vol.hasHalfMaps(), "Halves not registered.")
            half1, half2 = vol.getHalfMaps().split(',')
            self.assertTrue(exists(half1), msg="Average 1st half %s does not exists" % half1)
            self.assertTrue(exists(half2), msg="Average 2nd half %s does not exists" % half2)

        # hasctf
        self.assertEqual(hasCtf, vol.hasCTF(), "The volume doesn't have a CTF set.")
        # transform -> matrix
        self._checkVolumeOrigin(vol, expectedOriginShifts = expectedOriginShifts)
        # Ver si comprobar matriz

    def checkVolumeSet(self,
                       inVolumeSet: SetOfVolumes,
                       expectedSetSize: int,
                       expectedSRate: float,
                       expectedBoxSize: int,
                       expectedOriginShifts: Union[List[float], None] = None,
                       hasCtf: bool = False,
                       hasHalves: bool = False,
                       testAcqObj: Acquisition = None,
                       checkHeaderApix: bool = True,
                       sRateAngsPixTol: float = 0.01) -> None:

        # check header

        self.checkSetGeneralProps(inVolumeSet,
                                  expectedSetSize=expectedSetSize,
                                  expectedSRate=expectedSRate,
                                  sRateAngsPixTol=sRateAngsPixTol)

        # check hasctf
        self.assertEqual(hasCtf, inVolumeSet.hasCTF(), "The set of volumes doesn't have a CTF set.")
        # Comprobar acq del set
        #if testAcqObj is not None:
        #    self.checkTomoAcquisition()
        
        # Comprobar volumenes
        for volume in inVolumeSet:
            volId = volume.getClassId()
            print(f'---> checking the Volume Id = {volId}')
            self.checkVolume(
                vol = volume,
                expectedSRate = expectedSRate,
                expectedBoxSize = expectedBoxSize,
                hasCtf = hasCtf,
                hasHalves = hasHalves,
                sRateAngsPixTol = sRateAngsPixTol,
                expectedOriginShifts = expectedOriginShifts)
            
    def checkMovie(self,
                   mov: Movie,
                   samplingRate: float,
                   voltage: float,
                   moviesId: List[int] = [],
                   moviesNames: List[str] = [],
                   size: (float, float, float) = None,
                   dim: (int, int, int) = None) -> None:
        self.assertIsNotNone(mov)
        self.checkImage(mov)
        self.assertEqual(mov.getSize(), size)

        for i, m in enumerate(mov):
            if moviesId:
                self.assertEqual(m.getObjId(), moviesId[i])

            if moviesNames:
                self.assertEqual(os.path.basename(m.getFileName()),
                                 moviesNames[i])

            self.assertAlmostEqual(m.getSamplingRate(),
                                   samplingRate)
            a = m.getAcquisition()
            self.assertAlmostEqual(a.getVoltage(), voltage)

            if dim is not None:  # Check if dimensions are the expected ones
                x, y, n = m.getDim()
                self.assertEqual(dim, (x, y, n))


    def checkAcquisition(
        self,
        acquisition,
        voltage: Optional[float] = None,
        sphericalAberration: Optional[float] = None,
        amplitudeContrast: Optional[float] = None,
        magnification: Optional[float] = None,
        doseInitial: Optional[float] = None,
        dosePerFrame: Optional[float] = None):
        """
        Check that an acquisition has the expected parameters.
        
        :param acquisition: Acquisition object to check
        :param voltage: Expected voltage in kV
        :param sphericalAberration: Expected spherical aberration in mm
        :param amplitudeContrast: Expected amplitude contrast (0-1)
        :param magnification: Expected magnification
        """
        self.assertIsNotNone(acquisition, "Acquisition is None")
        
        if voltage is not None:
            self.assertAlmostEqual(
                acquisition.getVoltage(),
                voltage,
                places=1,
                msg=f"Voltage mismatch: expected {voltage}, got {acquisition.getVoltage()}"
            )
        
        if sphericalAberration is not None:
            self.assertAlmostEqual(
                acquisition.getSphericalAberration(),
                sphericalAberration,
                places=2,
                msg=f"Spherical aberration mismatch: expected {sphericalAberration}, got {acquisition.getSphericalAberration()}"
            )
        
        if amplitudeContrast is not None:
            self.assertAlmostEqual(
                acquisition.getAmplitudeContrast(),
                amplitudeContrast,
                places=2,
                msg=f"Amplitude contrast mismatch: expected {amplitudeContrast}, got {acquisition.getAmplitudeContrast()}"
            )
        
        if magnification is not None:
            self.assertAlmostEqual(
                acquisition.getMagnification(),
                magnification,
                places=1,
                msg=f"Magnification mismatch: expected {magnification}, got {acquisition.getMagnification()}"
            )

        if doseInitial is not None:
            self.assertAlmostEqual(
                acquisition.getDoseInitial(),
                doseInitial,
                places=2,
                msg=f"Initial dose mismatch: expected {doseInitial}, got {acquisition.getDoseInitial()}"
            )

        if dosePerFrame is not None:
            self.assertAlmostEqual(
                acquisition.getDosePerFrame(),
                dosePerFrame,
                places=2,
                msg=f"Dose per frame mismatch: expected {dosePerFrame}, got {acquisition.getDosePerFrame()}"
            )

    def checkTransform(
        self,
        shiftObject: Transform,
        shifts: Optional[Tuple[float, float, float]] = None):
        """
        Check that a transform has the expected shifts.
        
        :param shiftObject: Object to check
        :param shifts: Expected shifts as tuple (x, y, z) in Å
        """
        self.assertIsNotNone(shiftObject, "Transform is None")
        
        if shifts is not None:
            objectShift = shiftObject.getShifts()
            self.assertAlmostEqual(
                transformShifts[0],
                shifts[0],
                places=2,
                msg=f"Transform X shift mismatch: expected {shifts[0]}, got {transformShifts[0]}"
            )

            self.assertAlmostEqual(
                transformShifts[1],
                shifts[1],
                places=2,
                msg=f"Transform Y shift mismatch: expected {shifts[1]}, got {transformShifts[1]}"
            )

            self.assertAlmostEqual(
                transformShifts[2],
                shifts[2],
                places=2,
                msg=f"Transform Z shift mismatch: expected {shifts[2]}, got {transformShifts[2]}"
            )
    
    def checkImage(
        self, 
        img: Image,
        imageId: Optional[int] = None,
        samplingRate: Optional[float] = None,
        voltage: Optional[float] = None,
        sphericalAberration: Optional[float] = None,
        amplitudeContrast: Optional[float] = None,
        magnification: Optional[float] = None,
        doseInitial: Optional[float] = None,
        dosePerFrame: Optional[float] = None,
        dim: Optional[Tuple[int, int, int]] = None,
        transform: Optional[Tuple[float, float, float]] = None,
        origin: Optional[Tuple[float, float, float]] = None,
        hasCTF: bool = False):
        """
        Check that an image has the expected parameters.
        
        :param img: Image object to check
        :param imageId: Expected object ID
        :param samplingRate: Expected sampling rate in Å/px
        :param voltage: Expected voltage in kV
        :param sphericalAberration: Expected spherical aberration in mm
        :param amplitudeContrast: Expected amplitude contrast (0-1)
        :param magnification: Expected magnification
        :param dim: Expected dimensions as tuple (x, y, z)
        :param transform: Expected transform shifts as tuple (x, y, z) in Å
        :param origin: Expected origin shifts as tuple (x, y, z) in Å
        :param hasCTF: Whether image should have a CTF model
        """
        self.assertIsNotNone(img, "Image is None")
        
        # Check image ID
        if imageId is not None:
            self.assertEqual(
                img.getObjId(),
                imageId,
                f"Image ID mismatch: expected {imageId}, got {img.getObjId()}"
            )
        
        # Check filename
        if imageName is not None: # Se hace siempre
            # Comprobar que existe 
            self.assertEqual(
                os.path.basename(img.getFileName()),
                imageName,
                f"Image name mismatch: expected {imageName}, got {os.path.basename(img.getFileName())}"
            )
        
        # Check sampling rate
        if samplingRate is not None:
            self.assertAlmostEqual(
                img.getSamplingRate(),
                samplingRate,
                places=2,
                msg=f"Sampling rate mismatch: expected {samplingRate}, got {img.getSamplingRate()}"
            )
        
        # Check acquisition parameters
        if img.hasAcquisition():
            acquisition = img.getAcquisition()
            self.checkAcquisition(
                acquisition,
                voltage=voltage,
                sphericalAberration=sphericalAberration,
                amplitudeContrast=amplitudeContrast,
                magnification=magnification,
                doseInitial=doseInitial,
                dosePerFrame=dosePerFrame
            )

        # Check dimensions
        if dim is not None:
            imgDim = img.getDim()
            self.assertIsNotNone(imgDim, "Image dimensions are None")
            self.assertEqual( imgDim, dim, f"Dimension mismatch: expected {dim}, got {imgDim}")
        
        # Check transform shifts
        if transform is not None:
            self.assertTrue(img.hasTransform(), "Image should have a transform but doesn't")
            self.checkShift(transform=img.getTransform(), shifts=shifts)
        
        # Check origin shifts
        if origin is not None:
            self.assertTrue(img.hasOrigin(), "Image should have an origin but doesn't")
            self.checkShift(transform=img.getOrigin(), shifts=shifts)
            
        # Check CTF
        if hasCTF:
            self.assertTrue(img.hasCTF(), "Image should have a CTF model but doesn't")
        else:
            self.assertFalse(img.hasCTF(), "Image should not have a CTF model but does")

