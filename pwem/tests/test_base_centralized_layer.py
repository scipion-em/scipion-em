from typing import Optional, List, Union, Tuple
from os.path import exists

import mrcfile
import numpy as np
import os

from pyworkflow.tests import BaseTest
from pwem.objects import (Volume, Micrograph, Particle, Image, SetOfVolumes, SetOfImages, Acquisition, Movie,
                          Transform, SetOfParticles, Coordinate, SetOfCoordinates, CTFModel, AtomStruct,
                          SetOfClasses, SetOfAtomStructs, MovieAlignment, FramesRange, SetOfMovies)


class TestBaseCentralizedLayer(BaseTest):
    """
    
    """

    # ==========================================
    # 1. MÉTODOS DE CLASES PRINCIPALES
    # ==========================================

    def checkImage(self,
                   img: Image,
                   imageId: Optional[int] = None,
                   imageName: Optional[str] = None,
                   samplingRate: Optional[float] = None,
                   voltage: Optional[float] = None,
                   sphericalAberration: Optional[float] = None,
                   amplitudeContrast: Optional[float] = None,
                   magnification: Optional[float] = None,
                   doseInitial: Optional[float] = None,
                   dosePerFrame: Optional[float] = None,
                   dim: Optional[Tuple[int, int, int]] = None,
                   transformShifts: Optional[Tuple[float, float, float]] = None,
                   origin: Optional[Tuple[float, float, float]] = None,
                   hasCTF: bool = False,
                   sRateAngsPixTol: float = 0.01) -> None:
        """
        Validate that an image has all expected parameters.

        :param img: Image object to validate
        :param imageId: Expected image ID
        :param imageName: Expected image filename
        :param samplingRate: Expected sampling rate (Å/pix)
        :param voltage: Expected electron microscope voltage (kV)
        :param sphericalAberration: Expected spherical aberration (mm)
        :param amplitudeContrast: Expected amplitude contrast value
        :param magnification: Expected microscope magnification
        :param doseInitial: Expected initial dose (e⁻/Å²)
        :param dosePerFrame: Expected dose per frame (e⁻/Å²)
        :param dim: Expected image dimensions (x, y, z)
        :param transformShifts: Expected transformation shifts (x, y, z)
        :param origin: Expected origin shifts (x, y, z)
        :param hasCTF: Whether image should have CTF model
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)
        """
        # Validate object type
        if not isinstance(img, Image):
            self.fail(f"Expected Image object, got {type(img)}")

        # Check image ID
        if imageId is not None:
            print(f'---> Checking image ID = {imageId}')
            self.assertEqual(
                img.getObjId(),
                imageId,
                f"Image ID mismatch: expected {imageId}, got {img.getObjId()}"
            )

        # Check filename
        if imageName is not None:
            self.assertEqual(
                os.path.basename(img.getFileName()),
                imageName,
                f"Image name mismatch: expected {imageName}, got {os.path.basename(img.getFileName())}"
            )

        # Validate input parameters
        if samplingRate is not None:
            self.assertGreater(
                samplingRate, 0,
                f"Sampling rate must be positive, got {samplingRate}"
            )

        if dim is not None:
            self.assertEqual(
                len(dim), 3,
                f"Dimensions must have 3 elements (x, y, z), got {len(dim)}"
            )
            self.assertTrue(
                all(d > 0 for d in dim),
                f"All dimensions must be positive, got {dim}"
            )

        # Check sampling rate
        if samplingRate is not None:
            self.assertAlmostEqual(
                img.getSamplingRate(),
                samplingRate,
                places=2,
                msg=f"Sampling rate mismatch: expected {samplingRate}, got {img.getSamplingRate()}"
            )
            self.checkHeaderSRate(img, expectedSRate=samplingRate, sRateAngsPixTol=sRateAngsPixTol)

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
            self.assertIsNotNone(
                imgDim,
                f"Image dimensions are None for {img.getFileName()}"
            )
            self.assertEqual(
                imgDim, dim,
                f"Dimension mismatch: expected {dim}, got {imgDim}"
            )
        else:
            # Verify image has at least some dimensions
            imgDim = img.getDim()
            self.assertIsNotNone(
                imgDim,
                f"Image should have dimensions: {img.getFileName()}"
            )

        # Check transform shifts
        if transformShifts is not None:
            self.assertTrue(
                img.hasTransform(),
                "Image should have a transform but doesn't"
            )
            self.checkTransform(shiftObject=img.getTransform(), shifts=transformShifts)

        # Check origin shifts
        if origin is not None:
            self.assertTrue(
                img.hasOrigin(),
                "Image should have an origin but doesn't"
            )
            self.checkTransform(shiftObject=img.getOrigin(), shifts=origin)

        # Check CTF
        if hasCTF:
            self.assertTrue(
                img.hasCTF(),
                "Image should have a CTF model but doesn't"
            )
        else:
            self.assertFalse(
                img.hasCTF(),
                "Image should not have a CTF model but does"
            )

    def checkVolume(self,
                    vol: Volume,
                    expectedSRate: float = -1.,
                    expectedBoxSize: int = -1,
                    hasCtf: bool = False,
                    hasHalves: bool = True,
                    sRateAngsPixTol: float = 0.01,
                    expectedOriginShifts: Union[List[float], None] = None) -> None:
        """
        Validate the main properties of a volume.

        :param vol: Volume object to validate
        :param expectedSRate: Expected sampling rate (Å/pix), -1 to skip check
        :param expectedBoxSize: Expected box size in pixels, -1 to skip check
        :param hasCtf: Whether volume should have CTF model
        :param hasHalves: Whether volume should have half-maps
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)
        :param expectedOriginShifts: Expected origin shifts (x, y, z)
        """
        # Validate object type
        if not isinstance(vol, Volume):
            self.fail(f"Expected Volume object, got {type(vol)}")

        # Prepare dimensions if provided
        dim = None
        if expectedBoxSize != -1:
            dim = (expectedBoxSize, expectedBoxSize, expectedBoxSize)

        # Prepare origin if provided
        origin = None
        if expectedOriginShifts:
            origin = tuple(expectedOriginShifts)

        # Check shared properties using checkImage
        self.checkImage(
            img=vol,
            samplingRate=expectedSRate if expectedSRate != -1 else None,
            dim=dim,
            hasCTF=hasCtf,
            origin=origin
        )

        # Verify volume file exists
        volFileName = vol.getFileName()
        self.assertTrue(
            exists(volFileName),
            f"Volume {volFileName} does not exist"
        )

        # Check half-maps if required
        if hasHalves:
            self.assertTrue(
                vol.hasHalfMaps(),
                "Volume should have half-maps registered"
            )
            half1, half2 = vol.getHalfMaps().split(',')
            self.assertTrue(
                exists(half1),
                f"Volume first half {half1} does not exist"
            )
            self.assertTrue(
                exists(half2),
                f"Volume second half {half2} does not exist"
            )

    def checkMovie(self,
                   mov: Movie,
                   samplingRate: float,
                   voltage: float,
                   moviesId: Optional[List[int]] = None,
                   moviesNames: Optional[List[str]] = None,
                   size: Tuple[float, float, float] = None,
                   dim: Tuple[int, int, int] = None,
                   framesRange: Optional[Tuple[int, int, int]] = None,
                   numFrames: Optional[int] = None) -> None:
        """
        Validate a Movie object with all its properties.

        :param mov: Movie object to validate
        :param samplingRate: Expected sampling rate (Å/pix)
        :param voltage: Expected electron microscope voltage (kV)
        :param moviesId: List of expected movie IDs (optional)
        :param moviesNames: List of expected movie filenames (optional)
        :param size: Expected size dimensions (optional)
        :param dim: Expected dimensions (x, y, z) (optional)
        :param framesRange: Expected frames range (firstFrame, lastFrame, firstFrameIndex) (optional)
        :param numFrames: Expected number of frames (optional)
        """
        # Validate object type
        if not isinstance(mov, Movie):
            self.fail(f"Expected Movie object, got {type(mov)}")

        # Check shared image properties
        self.checkImage(mov, samplingRate=samplingRate, voltage=voltage)

        # Check size if provided
        if size is not None:
            self.assertEqual(
                mov.getSize(), size,
                f"Movie size mismatch: expected {size}, got {mov.getSize()}"
            )

        # Check sampling rate
        self.assertAlmostEqual(
            mov.getSamplingRate(), samplingRate,
            msg=f"Sampling rate mismatch: expected {samplingRate}, got {mov.getSamplingRate()}"
        )

        # Check acquisition voltage
        if mov.hasAcquisition():
            self.assertAlmostEqual(
                mov.getAcquisition().getVoltage(), voltage,
                msg=f"Voltage mismatch: expected {voltage}, got {mov.getAcquisition().getVoltage()}"
            )

        # Check dimensions if provided
        if dim is not None:
            self.assertEqual(
                mov.getDim(), dim,
                f"Dimension mismatch: expected {dim}, got {mov.getDim()}"
            )

        # Verify movie has FramesRange
        self.assertIsNotNone(
            mov.getFramesRange(),
            "Movie should have FramesRange"
        )

        # Check frames range details if provided
        if framesRange is not None:
            firstFrame, lastFrame, firstFrameIndex = framesRange
            actualRange = mov.getFramesRange()
            self.assertEqual(
                actualRange.getFirstFrame(), firstFrame,
                f"First frame mismatch: expected {firstFrame}, got {actualRange.getFirstFrame()}"
            )
            self.assertEqual(
                actualRange.getLastFrame(), lastFrame,
                f"Last frame mismatch: expected {lastFrame}, got {actualRange.getLastFrame()}"
            )
            self.assertEqual(
                actualRange.getFirstFrameIndex(), firstFrameIndex,
                f"First frame index mismatch: expected {firstFrameIndex}, got {actualRange.getFirstFrameIndex()}"
            )

        # Check number of frames if provided
        if numFrames is not None:
            actualNumFrames = mov.getNumberOfFrames()
            self.assertEqual(
                actualNumFrames, numFrames,
                f"Frame count mismatch: expected {numFrames}, got {actualNumFrames}"
            )

    def checkMicrograph(self,
                        mic: Micrograph,
                        micId: Optional[int] = None,
                        micName: Optional[str] = None,
                        imageName: Optional[str] = None,
                        samplingRate: Optional[float] = None,
                        voltage: Optional[float] = None,
                        sphericalAberration: Optional[float] = None,
                        amplitudeContrast: Optional[float] = None,
                        magnification: Optional[float] = None,
                        doseInitial: Optional[float] = None,
                        dosePerFrame: Optional[float] = None,
                        dim: Optional[Tuple[int, int, int]] = None,
                        transformShifts: Optional[Tuple[float, float, float]] = None,
                        origin: Optional[Tuple[float, float, float]] = None,
                        hasCTF: bool = False,
                        sRateAngsPixTol: float = 0.01) -> None:
        """
        Validate that a micrograph has all expected parameters.

        :param mic: Micrograph object to validate
        :param micId: Expected micrograph ID
        :param micName: Expected micrograph name
        :param imageName: Expected image filename
        :param samplingRate: Expected sampling rate (Å/pix)
        :param voltage: Expected electron microscope voltage (kV)
        :param sphericalAberration: Expected spherical aberration (mm)
        :param amplitudeContrast: Expected amplitude contrast value
        :param magnification: Expected microscope magnification
        :param doseInitial: Expected initial dose (e⁻/Å²)
        :param dosePerFrame: Expected dose per frame (e⁻/Å²)
        :param dim: Expected image dimensions (x, y, z)
        :param transformShifts: Expected transformation shifts (x, y, z)
        :param origin: Expected origin shifts (x, y, z)
        :param hasCTF: Whether micrograph should have CTF model
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)
        """
        # Validate object type
        if not isinstance(mic, Micrograph):
            self.fail(f"Expected Micrograph object, got {type(mic)}")

        # Check inherited image properties
        self.checkImage(
            img=mic,
            imageId=micId,
            imageName=imageName,
            samplingRate=samplingRate,
            voltage=voltage,
            sphericalAberration=sphericalAberration,
            amplitudeContrast=amplitudeContrast,
            magnification=magnification,
            doseInitial=doseInitial,
            dosePerFrame=dosePerFrame,
            dim=dim,
            transformShifts=transformShifts,
            origin=origin,
            hasCTF=hasCTF,
            sRateAngsPixTol=sRateAngsPixTol
        )

        # Check micrograph-specific name
        if micName is not None:
            self.assertEqual(
                mic.getMicName(), micName,
                f"Micrograph name mismatch: expected {micName}, got {mic.getMicName()}"
            )

    def checkAtomStruct(self,
                        atomStruct: AtomStruct,
                        hasVolume: bool = False,
                        pseudoatoms: bool = False) -> None:
        """
        Validate an AtomStruct (PDB) object.

        :param atomStruct: AtomStruct object to validate
        :param hasVolume: Whether atom structure should have associated volume
        :param pseudoatoms: Whether atom structure uses pseudoatoms representation
        """
        # Validate object type
        if not isinstance(atomStruct, AtomStruct):
            self.fail(f"Expected AtomStruct object, got {type(atomStruct)}")

        # Verify file exists
        self.assertTrue(
            exists(atomStruct.getFileName()),
            f"AtomStruct file {atomStruct.getFileName()} does not exist"
        )

        # Check volume association
        self.assertEqual(
            atomStruct.hasVolume(), hasVolume,
            f"HasVolume mismatch: expected {hasVolume}, got {atomStruct.hasVolume()}"
        )

        # Check pseudoatom representation
        self.assertEqual(
            atomStruct.getPseudoAtoms(), pseudoatoms,
            f"Pseudoatoms mismatch: expected {pseudoatoms}, got {atomStruct.getPseudoAtoms()}"
        )

    def checkCoordinate(self,
                        coord: Coordinate,
                        expectedX: int,
                        expectedY: int,
                        expectedMicId: Optional[int] = None) -> None:
        """
        Validate a Coordinate object.

        :param coord: Coordinate object to validate
        :param expectedX: Expected X coordinate (pixels)
        :param expectedY: Expected Y coordinate (pixels)
        :param expectedMicId: Expected micrograph ID (optional)
        """
        # Validate object type
        if not isinstance(coord, Coordinate):
            self.fail(f"Expected Coordinate object, got {type(coord)}")

        # Check X coordinate
        self.assertEqual(
            coord.getX(), expectedX,
            f"X coordinate mismatch: expected {expectedX}, got {coord.getX()}"
        )

        # Check Y coordinate
        self.assertEqual(
            coord.getY(), expectedY,
            f"Y coordinate mismatch: expected {expectedY}, got {coord.getY()}"
        )

        # Check micrograph ID if provided
        if expectedMicId is not None:
            self.assertEqual(
                coord.getMicId(), expectedMicId,
                f"Micrograph ID mismatch: expected {expectedMicId}, got {coord.getMicId()}"
            )

    def checkParticle(self,
                      particle: Particle,
                      classId: Optional[int] = None,
                      micId: Optional[int] = None,
                      particleId: Optional[int] = None,
                      imageName: Optional[str] = None,
                      samplingRate: Optional[float] = None,
                      voltage: Optional[float] = None,
                      sphericalAberration: Optional[float] = None,
                      amplitudeContrast: Optional[float] = None,
                      magnification: Optional[float] = None,
                      doseInitial: Optional[float] = None,
                      dosePerFrame: Optional[float] = None,
                      dim: Optional[Tuple[int, int, int]] = None,
                      transformShifts: Optional[Tuple[float, float, float]] = None,
                      origin: Optional[Tuple[float, float, float]] = None,
                      hasCTF: bool = False,
                      sRateAngsPixTol: float = 0.01,
                      corExpectedX: Optional[int] = None,
                      corExpectedY: Optional[int] = None,
                      corExpectedMicId: Optional[int] = None) -> None:
        """
        Validate a Particle object and its associated properties.

        Image Properties:
        :param particle: Particle object to validate
        :param particleId: Expected particle ID
        :param imageName: Expected image filename
        :param samplingRate: Expected sampling rate (Å/pix)
        :param voltage: Expected electron microscope voltage (kV)
        :param sphericalAberration: Expected spherical aberration (mm)
        :param amplitudeContrast: Expected amplitude contrast value
        :param magnification: Expected microscope magnification
        :param doseInitial: Expected initial dose (e⁻/Å²)
        :param dosePerFrame: Expected dose per frame (e⁻/Å²)
        :param dim: Expected image dimensions (x, y, z)
        :param transformShifts: Expected transformation shifts (x, y, z)
        :param origin: Expected origin shifts (x, y, z)
        :param hasCTF: Whether particle should have CTF model
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)

        Particle-specific Properties:
        :param classId: Expected class ID
        :param micId: Expected micrograph ID

        Coordinate Properties:
        :param corExpectedX: Expected X coordinate (pixels)
        :param corExpectedY: Expected Y coordinate (pixels)
        :param corExpectedMicId: Expected micrograph ID for coordinate
        """
        # Validate object type
        if not isinstance(particle, Particle):
            self.fail(f"Expected Particle object, got {type(particle)}")

        # Check particle class ID
        if classId is not None:
            self.assertEqual(
                particle.getClassId(), classId,
                f"Class ID mismatch: expected {classId}, got {particle.getClassId()}"
            )

        # Check associated micrograph ID
        if micId is not None:
            self.assertEqual(
                particle.getMicId(), micId,
                f"Micrograph ID mismatch: expected {micId}, got {particle.getMicId()}"
            )

        # Check inherited image properties
        self.checkImage(
            img=particle,
            imageId=particleId,
            imageName=imageName,
            samplingRate=samplingRate,
            voltage=voltage,
            sphericalAberration=sphericalAberration,
            amplitudeContrast=amplitudeContrast,
            magnification=magnification,
            doseInitial=doseInitial,
            dosePerFrame=dosePerFrame,
            dim=dim,
            transformShifts=transformShifts,
            origin=origin,
            hasCTF=hasCTF,
            sRateAngsPixTol=sRateAngsPixTol
        )

        # Check associated coordinate
        self.checkCoordinate(
            coord=particle.getCoordinate(),
            expectedX=corExpectedX,
            expectedY=corExpectedY,
            expectedMicId=corExpectedMicId
        )

    # ==========================================
    # 2. MÉTODOS DE SETS DE CLASES
    # ==========================================

    def checkImageSet(self,
                      inImageSet: SetOfImages,
                      expectedSetSize: int,
                      expectedSRate: float,
                      hasCtf: bool = False,
                      testAcqObj: Acquisition = None,
                      sRateAngsPixTol: float = 0.01,
                      checkItemsInSet: bool = false) -> None:
        """
        Validate a SetOfImages with its general properties.

        More generic than checkVolumeSet().

        :param inImageSet: SetOfImages object to validate
        :param expectedSetSize: Expected number of images in set
        :param expectedSRate: Expected sampling rate (Å/pix)
        :param hasCtf: Whether set should have CTF information
        :param testAcqObj: Acquisition object with expected parameters (optional)
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)
        """
        # Validate object type
        if not isinstance(inImageSet, SetOfImages):
            self.fail(f"Expected SetOfImages object, got {type(inImageSet)}")

        # Check general set properties
        self.checkSetGeneralProps(
            inSet=inImageSet,
            expectedSetSize=expectedSetSize,
            expectedSRate=expectedSRate,
            sRateAngsPixTol=sRateAngsPixTol
        )

        # Check CTF presence
        self.assertEqual(
            hasCtf, inImageSet.hasCTF(),
            "SetOfImages CTF presence mismatch"
        )

        # Check acquisition properties if provided
        if testAcqObj is not None and inImageSet.hasAcquisition():
            self.checkAcquisition(
                inImageSet.getAcquisition(),
                voltage=testAcqObj.getVoltage(),
                sphericalAberration=testAcqObj.getSphericalAberration(),
                amplitudeContrast=testAcqObj.getAmplitudeContrast(),
                magnification=testAcqObj.getMagnification()
            )

        # Validate each item in the set
        if checkItemsInSet:
            for image in inImageeSet:
                imId = image.getClassId()
                self.checkImage(image)

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
        """
        Validate a SetOfVolumes with all its properties.

        :param inVolumeSet: SetOfVolumes object to validate
        :param expectedSetSize: Expected number of volumes in set
        :param expectedSRate: Expected sampling rate (Å/pix)
        :param expectedBoxSize: Expected volume box size (pixels)
        :param expectedOriginShifts: Expected origin shifts (x, y, z) (optional)
        :param hasCtf: Whether set should have CTF information
        :param hasHalves: Whether volumes should have half-maps
        :param testAcqObj: Acquisition object with expected parameters (optional)
        :param checkHeaderApix: Whether to check voxel size in file headers
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)
        """
        # Validate object type
        if not isinstance(inVolumeSet, SetOfVolumes):
            self.fail(f"Expected SetOfVolumes object, got {type(inVolumeSet)}")

        # Check general set properties
        self.checkSetGeneralProps(
            inSet=inVolumeSet,
            expectedSetSize=expectedSetSize,
            expectedSRate=expectedSRate,
            sRateAngsPixTol=sRateAngsPixTol
        )

        # Check CTF presence
        self.assertEqual(
            hasCtf, inVolumeSet.hasCTF(),
            "SetOfVolumes CTF presence mismatch"
        )

        # Check acquisition properties if provided
        if testAcqObj is not None and inVolumeSet.hasAcquisition():
            self.checkAcquisition(
                inVolumeSet.getAcquisition(),
                voltage=testAcqObj.getVoltage(),
                sphericalAberration=testAcqObj.getSphericalAberration(),
                amplitudeContrast=testAcqObj.getAmplitudeContrast(),
                magnification=testAcqObj.getMagnification()
            )

        # Validate each volume in the set
        for volume in inVolumeSet:
            volId = volume.getClassId()
            print(f'---> Checking volume ID = {volId}')
            self.checkVolume(
                vol=volume,
                expectedSRate=expectedSRate,
                expectedBoxSize=expectedBoxSize,
                hasCtf=hasCtf,
                hasHalves=hasHalves,
                sRateAngsPixTol=sRateAngsPixTol,
                expectedOriginShifts=expectedOriginShifts
            )

    def checkSetOfMovies(self,
                         movieSet: SetOfMovies,
                         expectedSetSize: int,
                         expectedSRate: float,
                         expectedGain: Optional[str] = None,
                         expectedDark: Optional[str] = None,
                         sRateAngsPixTol: float = 0.01) -> None:
        """
        Validate a SetOfMovies with all its properties.

        :param movieSet: SetOfMovies object to validate
        :param expectedSetSize: Expected number of movies in set
        :param expectedSRate: Expected sampling rate (Å/pix)
        :param expectedGain: Expected gain reference file (optional)
        :param expectedDark: Expected dark reference file (optional)
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)
        """
        # Validate object type
        if not isinstance(movieSet, SetOfMovies):
            self.fail(f"Expected SetOfMovies object, got {type(movieSet)}")

        # Check general image set properties
        self.checkImageSet(
            inImageSet=movieSet,
            expectedSetSize=expectedSetSize,
            expectedSRate=expectedSRate,
            sRateAngsPixTol=sRateAngsPixTol
        )

        # Check gain reference if provided
        if expectedGain is not None:
            self.assertEqual(
                movieSet.getGain(), expectedGain,
                "Gain file mismatch"
            )

        # Check dark reference if provided
        if expectedDark is not None:
            self.assertEqual(
                movieSet.getDark(), expectedDark,
                "Dark file mismatch"
            )

        # Verify FramesRange of first movie
        if movieSet.getSize() > 0:
            firstMovie = movieSet.getFirstItem()
            self.assertIsNotNone(
                firstMovie.getFramesRange(),
                "First movie should have FramesRange"
            )

    def checkSetOfCoordinates(self,
                              coordSet: SetOfCoordinates,
                              expectedSize: int,
                              expectedBoxSize: Optional[int] = None) -> None:
        """
        Validate a SetOfCoordinates with all its properties.

        :param coordSet: SetOfCoordinates object to validate
        :param expectedSize: Expected number of coordinates in set
        :param expectedBoxSize: Expected box size for picking (optional)
        """
        # Validate object type
        if not isinstance(coordSet, SetOfCoordinates):
            self.fail(f"Expected SetOfCoordinates object, got {type(coordSet)}")

        # Check set size
        self.assertSetSize(coordSet, expectedSize)

        # Check box size if provided
        if expectedBoxSize is not None:
            self.assertEqual(
                coordSet.getBoxSize(), expectedBoxSize,
                f"Box size mismatch: expected {expectedBoxSize}, got {coordSet.getBoxSize()}"
            )

        # Validate each coordinate in the set
        for coord in coordSet.iterCoordinates():
            self.checkCoordinate(coord=coord)

    def checkSetOfClasses(self,
                          classesSet: SetOfClasses,
                          expectedSize: int,
                          hasRepresentatives: bool = True) -> None:
        """
        Validate a SetOfClasses with all its properties.

        :param classesSet: SetOfClasses object to validate
        :param expectedSize: Expected number of classes in set
        :param hasRepresentatives: Whether classes should have representative images
        """
        # Validate object type
        if not isinstance(classesSet, SetOfClasses):
            self.fail(f"Expected SetOfClasses object, got {type(classesSet)}")

        # Check set size
        self.assertSetSize(classesSet, expectedSize)

        # Check representative images presence
        self.assertEqual(
            classesSet.hasRepresentatives(), hasRepresentatives,
            f"HasRepresentatives mismatch: expected {hasRepresentatives}, got {classesSet.hasRepresentatives()}"
        )

    def checkSetOfParticles(self,
                            inParticleSet: SetOfParticles,
                            isSubparticle: Optional[bool] = None,
                            expectedSetSize: Optional[int] = None,
                            expectedSRate: Optional[float] = None,
                            hasCtf: Optional[bool] = False,
                            testAcqObj: Optional[Acquisition] = None,
                            sRateAngsPixTol: Optional[float] = 0.01,
                            expectedX: Optional[int] = None,
                            expectedY: Optional[int] = None,
                            expectedMicId: Optional[int] = None) -> None:
        """
        Validate a SetOfParticles with all its properties.

        :param inParticleSet: SetOfParticles object to validate
        :param isSubparticle: Whether particles are subparticles
        :param expectedSetSize: Expected number of particles in set
        :param expectedSRate: Expected sampling rate (Å/pix)
        :param hasCtf: Whether particles should have CTF information
        :param testAcqObj: Acquisition object with expected parameters (optional)
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)
        :param expectedX: Expected X coordinate for particles (pixels)
        :param expectedY: Expected Y coordinate for particles (pixels)
        :param expectedMicId: Expected micrograph ID for particles
        """
        # Validate object type
        if not isinstance(inParticleSet, SetOfParticles):
            self.fail(f"Expected SetOfParticles object, got {type(inParticleSet)}")

        # Check subparticle flag
        self.assertEqual(
            inParticleSet.getIsSubparticles(), isSubparticle,
            f"IsSubparticle mismatch: expected {isSubparticle}, got {inParticleSet.getIsSubparticles()}"
        )

        # Check general image set properties
        self.checkImageSet(
            inImageSet=inParticleSet,
            expectedSetSize=expectedSetSize,
            expectedSRate=expectedSRate,
            hasCtf=hasCtf,
            testAcqObj=testAcqObj,
            sRateAngsPixTol=sRateAngsPixTol
        )

        # Check coordinate properties
        self.checkCoordinate(
            coord=inParticleSet.getCoordinates(),
            expectedX=expectedX,
            expectedY=expectedY,
            expectedMicId=expectedMicId
        )

        # Validate each particle in the set
        for particle in inParticleSet:
            parId = particle.getClassId()
            # TODO: Uncomment and implement individual particle validation
            # print(f'---> Checking particle ID = {parId}')
            # self.checkParticle(particle, ...)

    # ==========================================
    # 3. MÉTODOS AUXILIARES O DE CLASES ACOMPAÑANTES
    # ==========================================

    def checkSetGeneralProps(self,
                             inSet,
                             expectedSetSize: int,
                             expectedSRate: float,
                             streamState: Optional[int] = None,
                             sRateAngsPixTol: float = 0.01) -> None:
        """
        Validate general properties of a Scipion set object.

        :param inSet: Scipion set object to validate
        :param expectedSetSize: Expected number of items in set
        :param expectedSRate: Expected sampling rate (Å/pix)
        :param streamState: Expected stream state (None = don't check, 2 = closed stream)
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)
        """
        # Check set size if greater than 0
        if expectedSetSize > 0:
            self.assertSetSize(inSet, expectedSetSize)

        # Check sampling rate
        self.assertAlmostEqual(
            inSet.getSamplingRate(), expectedSRate,
            delta=sRateAngsPixTol,
            msg=f"Sampling rate mismatch: expected {expectedSRate} ± {sRateAngsPixTol}, got {inSet.getSamplingRate()}"
        )

        # Check stream state if provided
        if streamState is not None:
            self.assertEqual(
                inSet.getStreamState(), streamState,
                f"Stream state mismatch: expected {streamState}, got {inSet.getStreamState()}"
            )

        # Verify properties are persisted in database
        self.assertTrue(
            inSet.hasProperty("self"),
            f"Set {inSet.getFileName()} does not have 'self' in properties table. "
            f"Properties may not be persisted correctly."
        )

    def checkTransform(self,
                       shiftObject: Transform,
                       shifts: Optional[Tuple[float, float, float]] = None,
                       places: int = 2) -> None:
        """
        Validate shift values in a Transform object.

        NOTE: This validates SHIFT VALUES (x, y, z), not the transformation matrix.
        For matrix validation, use checkTransformMatrix().

        :param shiftObject: Transform object to validate
        :param shifts: Expected shift values (x, y, z)
        :param places: Decimal places for floating-point comparison (default: 2)
        """
        # Validate object exists
        self.assertIsNotNone(shiftObject, "Transform is None")

        # Check shift values if provided
        if shifts is not None:
            transformShifts = shiftObject.getShifts()
            for expected, actual in zip(shifts, transformShifts):
                self.assertAlmostEqual(
                    actual, expected,
                    places=places,
                    msg=f"Shift mismatch: expected {expected}, got {actual}"
                )

    def checkTransformMatrix(self,
                             outMatrix: np.ndarray,
                             alignment: bool = False,
                             is2d: bool = False) -> None:
        """
        Validate the shape and contents of a transformation matrix.

        NOTE: This validates the MATRIX itself, not shift values.
        For shift validation, use checkTransform().

        :param outMatrix: Transformation matrix to validate (numpy array)
        :param alignment: Whether matrix represents alignment (non-identity) or default (identity)
        :param is2d: Whether to expect 3x3 (True) or 4x4 (False) matrix
        """
        # Determine expected matrix size
        size = 3 if is2d else 4
        transfMatrixShape = (size, size)
        identityMatrix = np.eye(size)

        # Validate matrix object
        self.assertIsNotNone(outMatrix, "Transform matrix is None")

        # Convert to numpy array if needed
        if type(outMatrix) is not np.ndarray:
            outMatrix = np.array(outMatrix)

        self.assertIsNotNone(outMatrix)

        # Check matrix dimensions
        self.assertEqual(
            outMatrix.shape, transfMatrixShape,
            f"Matrix shape mismatch: expected {transfMatrixShape}, got {outMatrix.shape}"
        )

        # Check matrix content (identity vs. alignment)
        if alignment:
            self.assertFalse(
                np.array_equal(outMatrix, identityMatrix),
                "Matrix should represent alignment but is identity"
            )
        else:
            self.assertTrue(
                np.array_equal(outMatrix, identityMatrix),
                "Matrix should be identity but represents alignment"
            )

    def checkMovieAlignment(self,
                            alignment: MovieAlignment,
                            expectedFirst: int,
                            expectedLast: int,
                            expectedXShifts: Optional[List[float]] = None,
                            expectedYShifts: Optional[List[float]] = None) -> None:
        """
        Validate a MovieAlignment object.

        :param alignment: MovieAlignment object to validate
        :param expectedFirst: Expected first frame number
        :param expectedLast: Expected last frame number
        :param expectedXShifts: Expected X-axis shift values for each frame (optional)
        :param expectedYShifts: Expected Y-axis shift values for each frame (optional)
        """
        # Validate object exists
        self.assertIsNotNone(alignment, "MovieAlignment is None")

        # Check frame range
        first, last = alignment.getRange()
        self.assertEqual(
            first, expectedFirst,
            f"First frame mismatch: expected {expectedFirst}, got {first}"
        )
        self.assertEqual(
            last, expectedLast,
            f"Last frame mismatch: expected {expectedLast}, got {last}"
        )

        # Check shift values if provided
        if expectedXShifts is not None or expectedYShifts is not None:
            xShifts, yShifts = alignment.getShifts()

            # Validate X shifts
            if expectedXShifts is not None:
                self.assertEqual(
                    len(xShifts), len(expectedXShifts),
                    f"X shifts count mismatch: expected {len(expectedXShifts)}, got {len(xShifts)}"
                )
                for i, (expected, actual) in enumerate(zip(expectedXShifts, xShifts)):
                    self.assertAlmostEqual(
                        actual, expected, places=2,
                        msg=f"X shift at frame {i} mismatch: expected {expected}, got {actual}"
                    )

            # Validate Y shifts
            if expectedYShifts is not None:
                self.assertEqual(
                    len(yShifts), len(expectedYShifts),
                    f"Y shifts count mismatch: expected {len(expectedYShifts)}, got {len(yShifts)}"
                )
                for i, (expected, actual) in enumerate(zip(expectedYShifts, yShifts)):
                    self.assertAlmostEqual(
                        actual, expected, places=2,
                        msg=f"Y shift at frame {i} mismatch: expected {expected}, got {actual}"
                    )

    def checkHeaderSRate(self,
                         inObj: Union[SetOfVolumes, Volume, SetOfImages, Image],
                         expectedSRate: float,
                         sRateAngsPixTol: float = 0.01) -> None:
        """
        Validate the sampling rate in file headers of volumes or images.

        Reads MRC file headers and verifies voxel size values.

        :param inObj: Object (Volume, Image, SetOfVolumes, SetOfImages) to validate
        :param expectedSRate: Expected sampling rate (Å/pix)
        :param sRateAngsPixTol: Tolerance for sampling rate comparison (Å/pix)
        """
        # Get file path
        volumeFile = inObj.getFileName()

        # Skip validation if file doesn't exist or is virtual
        if not volumeFile or not os.path.exists(volumeFile):
            return

        # Only validate MRC format files
        if not volumeFile.endswith(".mrc") and not volumeFile.endswith(".mrcs"):
            return

        # Read MRC header and validate voxel size
        with mrcfile.open(volumeFile, permissive=True, header_only=True) as mrc:
            vs = mrc.voxel_size

            # Determine which voxel size values to check
            if isinstance(inObj, (Volume, SetOfVolumes)):
                vs_values = [float(vs.x), float(vs.y), float(vs.z)]
            else:
                vs_values = [float(vs.x), float(vs.y)]

            # Validate each voxel size value
            for voxelSize in vs_values:
                self.assertAlmostEqual(
                    voxelSize, expectedSRate,
                    delta=sRateAngsPixTol,
                    msg=f"Voxel size mismatch in {volumeFile}\n"
                        f"Expected: {expectedSRate} Å/pix (±{sRateAngsPixTol})\n"
                        f"Got: {voxelSize} Å/pix\n"
                        f"File header values: {vs}"
                )

    def checkAcquisition(self,
                         acquisition: Acquisition,
                         voltage: Optional[float] = None,
                         sphericalAberration: Optional[float] = None,
                         amplitudeContrast: Optional[float] = None,
                         magnification: Optional[float] = None,
                         doseInitial: Optional[float] = None,
                         dosePerFrame: Optional[float] = None) -> None:
        """
        Validate microscope acquisition parameters.

        :param acquisition: Acquisition object to validate
        :param voltage: Expected electron microscope voltage (kV)
        :param sphericalAberration: Expected spherical aberration (mm)
        :param amplitudeContrast: Expected amplitude contrast value [0-1]
        :param magnification: Expected objective magnification
        :param doseInitial: Expected initial dose (e⁻/Å²)
        :param dosePerFrame: Expected dose per frame (e⁻/Å²)
        """
        # Validate object exists
        self.assertIsNotNone(acquisition, "Acquisition is None")

        # Check voltage
        if voltage is not None:
            self.assertAlmostEqual(
                acquisition.getVoltage(), voltage,
                places=1,
                msg=f"Voltage mismatch: expected {voltage} kV, got {acquisition.getVoltage()} kV"
            )

        # Check spherical aberration
        if sphericalAberration is not None:
            self.assertAlmostEqual(
                acquisition.getSphericalAberration(), sphericalAberration,
                places=2,
                msg=f"Spherical aberration mismatch: expected {sphericalAberration} mm, got {acquisition.getSphericalAberration()} mm"
            )

        # Check amplitude contrast
        if amplitudeContrast is not None:
            self.assertAlmostEqual(
                acquisition.getAmplitudeContrast(), amplitudeContrast,
                places=2,
                msg=f"Amplitude contrast mismatch: expected {amplitudeContrast}, got {acquisition.getAmplitudeContrast()}"
            )

        # Check magnification
        if magnification is not None:
            self.assertAlmostEqual(
                acquisition.getMagnification(), magnification,
                places=1,
                msg=f"Magnification mismatch: expected {magnification}x, got {acquisition.getMagnification()}x"
            )

        # Check initial dose
        if doseInitial is not None:
            self.assertAlmostEqual(
                acquisition.getDoseInitial(), doseInitial,
                places=2,
                msg=f"Initial dose mismatch: expected {doseInitial} e⁻/Å², got {acquisition.getDoseInitial()} e⁻/Å²"
            )

        # Check dose per frame
        if dosePerFrame is not None:
            self.assertAlmostEqual(
                acquisition.getDosePerFrame(), dosePerFrame,
                places=2,
                msg=f"Dose per frame mismatch: expected {dosePerFrame} e⁻/Å², got {acquisition.getDosePerFrame()} e⁻/Å²"
            )

    def checkCTF(self,
                 ctf: CTFModel,
                 defocusU: float,
                 defocusV: float,
                 defocusAngle: float,
                 resolution: Optional[float] = None) -> None:
        """
        Validate CTF (Contrast Transfer Function) model parameters.

        :param ctf: CTFModel object to validate
        :param defocusU: Expected defocus in U direction (Å)
        :param defocusV: Expected defocus in V direction (Å)
        :param defocusAngle: Expected defocus angle (degrees)
        :param resolution: Expected CTF resolution (Å) (optional)
        """
        # Validate object exists
        self.assertIsNotNone(ctf, "CTFModel is None")

        # Check defocus U
        self.assertAlmostEqual(
            ctf.getDefocusU(), defocusU,
            places=1,
            msg=f"DefocusU mismatch: expected {defocusU} Å, got {ctf.getDefocusU()} Å"
        )

        # Check defocus V
        self.assertAlmostEqual(
            ctf.getDefocusV(), defocusV,
            places=1,
            msg=f"DefocusV mismatch: expected {defocusV} Å, got {ctf.getDefocusV()} Å"
        )

        # Check defocus angle
        self.assertAlmostEqual(
            ctf.getDefocusAngle(), defocusAngle,
            places=1,
            msg=f"DefocusAngle mismatch: expected {defocusAngle}°, got {ctf.getDefocusAngle()}°"
        )

        # Check resolution if provided
        if resolution is not None:
            self.assertAlmostEqual(
                ctf.getResolution(), resolution,
                places=2,
                msg=f"Resolution mismatch: expected {resolution} Å, got {ctf.getResolution()} Å"
            )

    # REFERENCIA PARA HACER
    #def checkExtracted3dCoordinates(self,
    #                                inSet: SetOfCoordinates3D,
    #                                outCoords: SetOfCoordinates3D,
    #                                expectedSetSize: int = -1,
    #                                expectedBoxSize: int = -1,
    #                                expectedSRate: float = -1.0,
    #                                convention: Union[str, None] = TR_SCIPION,
    #                                orientedParticles: bool = False) -> None:
    #    """Checks the results of a coordinate extraction protocol.

    #    :param inSet: input set from which the coordinates were extracted. It can be a SetOf3DCoordinates or a
    #    SetOfSubTomograms.
    #    :param outCoords: the resulting SetOf3DCoordinates after the coordinate extraction.
    #    :param expectedSetSize: expected set site to check.
    #    :param expectedBoxSize: expected box size, in pixels, to check.
    #    :param expectedSRate: expected sampling rate, in Å/pix, to check.
    #    :param convention: TR_SCIPION by default. Convention of the coordinates. See scipion-em-tomo/tomo/constants.py.
    #    :param orientedParticles: False by default. Used to specify if the expected transformation matrix should be
    #    and eye matrix (False) or not (True)."""
    #    if type(inSet) == SetOfSubTomograms:
    #        inSet = inSet.getCoordinates3D()
    #    # First, check the set size, sampling rate, and box size
    #    self.checkCoordsOrPartsSetGeneralProps(outCoords,
    #                                           expectedSetSize=expectedSetSize,
    #                                           expectedSRate=expectedSRate,
    #                                           expectedBoxSize=expectedBoxSize)
    #    # Check the coordinate extremes
    #    inCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(inSet)
    #    outCoordsExtremes = self.getMinAndMaxCoordValuesFromSet(outCoords)
    #    coordScaleFactor = inSet.getSamplingRate() / outCoords.getSamplingRate()
    #    shiftsScaleFactor = 1 / coordScaleFactor
    #    self.assertTrue(np.array_equal(outCoordsExtremes, coordScaleFactor * inCoordsExtremes))
    #    # Other checks per coordinate
    #    for inElement, outCoord in zip(inSet, outCoords):
    #        # Check the transformation matrices and shifts
    #        inSetTrMatrix = inElement.getMatrix(convention=convention)
    #        outCoordTrMatrix = outCoord.getMatrix(convention=convention)
    #        self.checkTransformMatrix(outCoordTrMatrix, alignment=orientedParticles)
    #        self.checkShiftsScaling(inSetTrMatrix, outCoordTrMatrix, shiftsScaleFactor)
    #        # Check the tomoId
    #        self.assertEqual(outCoord.getTomoId(), inElement.getTomoId())

    # def _checkVolumeOrigin(self,
    #                        vol: Volume,
    #                        expectedOriginShifts: Union[List[float], None] = None,
    #                        places: int = 1) -> None:
    #     """
    #     Check volume origin shifts.
    #     :param vol: Volume object
    #     :param expectedOriginShifts: Expected origin shifts
    #     :param places: Decimal places for comparison (default: 1)
    #     """
    #     if expectedOriginShifts is not None:
    #         x, y, z = vol.getOrigin().getShifts()
    #         for i, j in zip([x, y, z], expectedOriginShifts):
    #             self.assertAlmostEqual(
    #                 i, j, places=places,
    #                 msg="Expected and resulting volume shifts are different"
    #             )