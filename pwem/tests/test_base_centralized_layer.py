# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *
# **************************************************************************
"""
Centralized test layer for pwem base data objects.

Provides reusable test methods for validating Scipion EM data objects
and their corresponding sets, meant to be inherited by protocol-specific
test classes.

Methods are organized into three sections:

    1. AUXILIARY METHODS -- shared validation utilities.
    2. INDIVIDUAL OBJECT CHECKS -- per-type validation of a single object.
    3. SET CHECKS -- set-level validation + iteration over items.
"""
from os.path import exists
from typing import List, Optional, Tuple, Union

import mrcfile
import numpy as np

from pyworkflow.tests import BaseTest
from pwem.objects import (Acquisition, AtomStruct, CTFModel, Coordinate,
                           Image, Mask, Micrograph, Movie, MovieAlignment,
                           Particle, SetOfClasses, SetOfCoordinates, SetOfCTF,
                           SetOfImages, SetOfMicrographs, SetOfMovies,
                           SetOfParticles, SetOfVolumes, Transform, Volume,
                           VolumeMask)


class TestBaseCentralizedLayer(BaseTest):

    # =========================================================================
    # 1. AUXILIARY METHODS
    # =========================================================================

    def checkSetGeneralProps(self, inSet,
                             expectedSetSize: int,
                             expectedSRate: Optional[float] = None,
                             streamState: Optional[int] = None,
                             sRateAngsPixTol: float = 0.01) -> None:
        """Validate general properties of a Scipion set.

        Parameters
        ----------
        inSet : Set
            Scipion set object to validate.
        expectedSetSize : int
            Expected number of items (checked only when > 0).
        expectedSRate : float, optional
            Expected sampling rate in A/pix. ``None`` to skip.
        streamState : int, optional
            Expected stream state (2 = closed). ``None`` to skip.
        sRateAngsPixTol : float
            Tolerance for the sampling-rate comparison.
        """
        if expectedSetSize > 0:
            self.assertSetSize(inSet, expectedSetSize)
        if expectedSRate is not None:
            self.assertAlmostEqual(inSet.getSamplingRate(), expectedSRate,
                                   delta=sRateAngsPixTol,
                                   msg=f"Sampling rate mismatch: "
                                       f"expected {expectedSRate} +/- {sRateAngsPixTol}, "
                                       f"got {inSet.getSamplingRate()}")
        if streamState is not None:
            self.assertEqual(inSet.getStreamState(), streamState,
                             msg=f"Stream state mismatch: "
                                 f"expected {streamState}, "
                                 f"got {inSet.getStreamState()}")
        self.assertTrue(inSet.hasProperty("self"),
                        msg=f"Set {inSet.getFileName()} does not have 'self' "
                            f"in the properties table -- properties may not be "
                            f"persisted correctly.")

    def checkTransform(self,
                       shiftObject: Transform,
                       shifts: Optional[Tuple[float, float, float]] = None,
                       places: int = 2) -> None:
        """Validate shift values stored in a Transform object.

        Parameters
        ----------
        shiftObject : Transform
            Transform object to validate (its ``getShifts()`` is checked).
        shifts : tuple of float, optional
            Expected (x, y, z) shifts. ``None`` to skip.
        places : int
            Decimal places for the floating-point comparison.
        """
        self.assertIsNotNone(shiftObject, "Transform object is None.")
        if shifts is not None:
            actual = shiftObject.getShifts()
            for exp, act in zip(shifts, actual):
                self.assertAlmostEqual(act, exp, places=places,
                                       msg=f"Shift mismatch: "
                                           f"expected {exp}, got {act}")

    def checkTransformMatrix(self,
                             outMatrix: np.ndarray,
                             alignment: bool = False,
                             is2d: bool = False) -> None:
        """Validate the shape and coarse content of a transformation matrix.

        Parameters
        ----------
        outMatrix : numpy.ndarray
            Transformation matrix to validate.
        alignment : bool
            If ``True`` the matrix must **not** be the identity.
            If ``False`` it must **be** the identity.
        is2d : bool
            If ``True`` a 3x3 matrix is expected; otherwise 4x4.
        """
        size = 3 if is2d else 4
        expShape = (size, size)
        identity = np.eye(size)
        self.assertIsNotNone(outMatrix, "Transform matrix is None.")
        if not isinstance(outMatrix, np.ndarray):
            outMatrix = np.array(outMatrix)
        self.assertIsNotNone(outMatrix)
        self.assertEqual(outMatrix.shape, expShape,
                         msg=f"Matrix shape mismatch: "
                             f"expected {expShape}, got {outMatrix.shape}")
        if alignment:
            self.assertFalse(np.array_equal(outMatrix, identity),
                             msg="Matrix should represent an alignment "
                                 "(non-identity) but is identity.")
        else:
            self.assertTrue(np.array_equal(outMatrix, identity),
                            msg="Matrix should be identity but represents "
                                "an alignment (non-identity).")

    def checkHeaderSRate(self,
                         inObj: Union['SetOfVolumes', Volume, 'SetOfImages', Image],
                         expectedSRate: float,
                         sRateAngsPixTol: float = 0.01) -> None:
        """Validate the sampling rate stored in the MRC header of a file.

        Parameters
        ----------
        inObj : Volume, Image, SetOfVolumes or SetOfImages
            Object whose backing file's header will be inspected.
        expectedSRate : float
            Expected sampling rate in A/pix.
        sRateAngsPixTol : float
            Tolerance for the comparison.
        """
        fn = inObj.getFileName()
        if not fn or not exists(fn):
            return
        if not fn.endswith(('.mrc', '.mrcs')):
            return
        with mrcfile.open(fn, permissive=True, header_only=True) as mrc:
            vs = mrc.voxel_size
            vals = [float(vs.x), float(vs.y), float(vs.z)] \
                if isinstance(inObj, (Volume, SetOfVolumes)) \
                else [float(vs.x), float(vs.y)]
            for v in vals:
                self.assertAlmostEqual(v, expectedSRate,
                                       delta=sRateAngsPixTol,
                                       msg=f"Header voxel size mismatch in "
                                           f"{fn}: expected {expectedSRate} "
                                           f"+/- {sRateAngsPixTol}, got {v}. "
                                           f"Full header values: {vs}")

    def checkAcquisition(self,
                         acquisition: Acquisition,
                         voltage: Optional[float] = None,
                         sphericalAberration: Optional[float] = None,
                         amplitudeContrast: Optional[float] = None,
                         magnification: Optional[float] = None,
                         doseInitial: Optional[float] = None,
                         dosePerFrame: Optional[float] = None) -> None:
        """Validate microscope acquisition parameters.

        All parameters are optional; only those provided are tested.

        Parameters
        ----------
        acquisition : Acquisition
            Acquisition object to validate.
        voltage : float, optional
            Expected microscope voltage (kV).
        sphericalAberration : float, optional
            Expected spherical aberration (mm).
        amplitudeContrast : float, optional
            Expected amplitude contrast.
        magnification : float, optional
            Expected magnification.
        doseInitial : float, optional
            Expected initial dose (e-/A2).
        dosePerFrame : float, optional
            Expected dose per frame (e-/A2).
        """
        self.assertIsNotNone(acquisition, "Acquisition is None.")
        if voltage is not None:
            self.assertAlmostEqual(acquisition.getVoltage(), voltage,
                                   delta=1,
                                   msg=f"Voltage mismatch: "
                                       f"expected {voltage}, "
                                       f"got {acquisition.getVoltage()}")
        if sphericalAberration is not None:
            self.assertAlmostEqual(acquisition.getSphericalAberration(),
                                   sphericalAberration, delta=0.01,
                                   msg=f"Spherical aberration mismatch: "
                                       f"expected {sphericalAberration}, "
                                       f"got {acquisition.getSphericalAberration()}")
        if amplitudeContrast is not None:
            self.assertAlmostEqual(acquisition.getAmplitudeContrast(),
                                   amplitudeContrast, delta=0.01,
                                   msg=f"Amplitude contrast mismatch: "
                                       f"expected {amplitudeContrast}, "
                                       f"got {acquisition.getAmplitudeContrast()}")
        if magnification is not None:
            self.assertAlmostEqual(acquisition.getMagnification(),
                                   magnification, delta=1,
                                   msg=f"Magnification mismatch: "
                                       f"expected {magnification}, "
                                       f"got {acquisition.getMagnification()}")
        if doseInitial is not None:
            self.assertAlmostEqual(acquisition.getDoseInitial(),
                                   doseInitial, delta=0.01,
                                   msg=f"Initial dose mismatch: "
                                       f"expected {doseInitial}, "
                                       f"got {acquisition.getDoseInitial()}")
        if dosePerFrame is not None:
            self.assertAlmostEqual(acquisition.getDosePerFrame(),
                                   dosePerFrame, delta=0.01,
                                   msg=f"Dose per frame mismatch: "
                                       f"expected {dosePerFrame}, "
                                       f"got {acquisition.getDosePerFrame()}")

    def checkCTF(self,
                 ctf: CTFModel,
                 defocusU: float,
                 defocusV: float,
                 defocusAngle: float,
                 resolution: Optional[float] = None,
                 phaseShift: Optional[float] = None) -> None:
        """Validate CTF model parameters.

        Parameters
        ----------
        ctf : CTFModel
            CTF model to validate.
        defocusU : float
            Expected defocus in U (A).
        defocusV : float
            Expected defocus in V (A).
        defocusAngle : float
            Expected defocus angle (degrees).
        resolution : float, optional
            Expected resolution (A).
        phaseShift : float, optional
            Expected phase shift (degrees).
        """
        self.assertIsNotNone(ctf, "CTFModel is None.")
        if not isinstance(ctf, CTFModel):
            self.fail(f"Expected CTFModel, got {type(ctf)}.")
        self.assertAlmostEqual(ctf.getDefocusU(), defocusU, delta=1,
                               msg=f"DefocusU mismatch: "
                                   f"expected {defocusU}, "
                                   f"got {ctf.getDefocusU()}")
        self.assertAlmostEqual(ctf.getDefocusV(), defocusV, delta=1,
                               msg=f"DefocusV mismatch: "
                                   f"expected {defocusV}, "
                                   f"got {ctf.getDefocusV()}")
        self.assertAlmostEqual(ctf.getDefocusAngle(), defocusAngle,
                               delta=1,
                               msg=f"DefocusAngle mismatch: "
                                   f"expected {defocusAngle}, "
                                   f"got {ctf.getDefocusAngle()}")
        if resolution is not None:
            self.assertAlmostEqual(ctf.getResolution(), resolution,
                                   delta=0.01,
                                   msg=f"Resolution mismatch: "
                                       f"expected {resolution}, "
                                       f"got {ctf.getResolution()}")
        if phaseShift is not None:
            self.assertTrue(ctf.hasPhaseShift(),
                            msg="CTF has no phase shift but one was expected.")
            self.assertAlmostEqual(ctf.getPhaseShift(), phaseShift,
                                   delta=1,
                                   msg=f"PhaseShift mismatch: "
                                       f"expected {phaseShift}, "
                                       f"got {ctf.getPhaseShift()}")

    def checkMovieAlignment(self,
                            alignment: MovieAlignment,
                            expectedFirst: int,
                            expectedLast: int,
                            expectedXShifts: Optional[List[float]] = None,
                            expectedYShifts: Optional[List[float]] = None):
        """Validate a MovieAlignment object.

        Parameters
        ----------
        alignment : MovieAlignment
            The alignment object to validate.
        expectedFirst : int
            Expected first frame used for alignment.
        expectedLast : int
            Expected last frame used for alignment.
        expectedXShifts : list of float, optional
            Expected X shifts per frame.
        expectedYShifts : list of float, optional
            Expected Y shifts per frame.
        """
        self.assertIsNotNone(alignment, "MovieAlignment is None.")
        first, last = alignment.getRange()
        self.assertEqual(first, expectedFirst,
                         msg=f"First frame mismatch: "
                             f"expected {expectedFirst}, got {first}")
        self.assertEqual(last, expectedLast,
                         msg=f"Last frame mismatch: "
                             f"expected {expectedLast}, got {last}")
        if expectedXShifts is not None or expectedYShifts is not None:
            xShifts, yShifts = alignment.getShifts()
            if expectedXShifts is not None:
                self.assertEqual(len(xShifts), len(expectedXShifts),
                                 msg=f"X-shifts length mismatch: "
                                     f"expected {len(expectedXShifts)}, "
                                     f"got {len(xShifts)}")
                for i, (exp, act) in enumerate(zip(expectedXShifts, xShifts)):
                    self.assertAlmostEqual(act, exp, places=2,
                                           msg=f"X shift [{i}] mismatch: "
                                               f"expected {exp}, got {act}")
            if expectedYShifts is not None:
                self.assertEqual(len(yShifts), len(expectedYShifts),
                                 msg=f"Y-shifts length mismatch: "
                                     f"expected {len(expectedYShifts)}, "
                                     f"got {len(yShifts)}")
                for i, (exp, act) in enumerate(zip(expectedYShifts, yShifts)):
                    self.assertAlmostEqual(act, exp, places=2,
                                           msg=f"Y shift [{i}] mismatch: "
                                               f"expected {exp}, got {act}")

    def checkObjectEnabled(self, obj, isExcluded: bool, objLabel: str = ""):
        """Assert that an object's enabled/disabled state is as expected.

        Parameters
        ----------
        obj : EMObject
            Object whose ``isEnabled()`` status is checked.
        isExcluded : bool
            ``True`` if the object is expected to be disabled (excluded).
        objLabel : str
            Optional label for error messages.
        """
        enabled = obj.isEnabled()
        if isExcluded:
            self.assertFalse(enabled,
                             msg=f"{objLabel}: object expected to be "
                                 f"disabled, but it is enabled.")
        else:
            self.assertTrue(enabled,
                            msg=f"{objLabel}: object expected to be "
                                f"enabled, but it is disabled.")

    # =========================================================================
    # 2. INDIVIDUAL OBJECT CHECKS
    # =========================================================================

    def checkImage(self,
                   img: Image,
                   imageId: Optional[int] = None,
                   imageName: Optional[str] = None,
                   samplingRate: Optional[float] = None,
                   dim: Optional[Tuple[int, int, int]] = None,
                   hasCTF: Optional[bool] = None,
                   transformShifts: Optional[Tuple[float, float, float]] = None,
                   origin: Optional[Tuple[float, float, float]] = None,
                   voltage: Optional[float] = None,
                   sphericalAberration: Optional[float] = None,
                   amplitudeContrast: Optional[float] = None,
                   magnification: Optional[float] = None,
                   doseInitial: Optional[float] = None,
                   dosePerFrame: Optional[float] = None,
                   sRateAngsPixTol: float = 0.01,
                   checkHeaderApix: bool = True) -> None:
        """Validate an Image object and its associated attributes.

        Parameters
        ----------
        img : Image
            Image object to validate.
        imageId : int, optional
            Expected object ID.
        imageName : str, optional
            Expected base filename.
        samplingRate : float, optional
            Expected sampling rate (A/pix).
        dim : tuple of int, optional
            Expected (x, y, z) dimensions.
        hasCTF : bool, optional
            ``True``/``False`` to enforce CTF presence/absence;
            ``None`` to skip.
        transformShifts : tuple of float, optional
            Expected (x, y, z) transform shifts.
        origin : tuple of float, optional
            Expected (x, y, z) origin shifts.
        voltage : float, optional
            Expected microscope voltage (kV).
        sphericalAberration : float, optional
            Expected spherical aberration (mm).
        amplitudeContrast : float, optional
            Expected amplitude contrast.
        magnification : float, optional
            Expected magnification.
        doseInitial : float, optional
            Expected initial dose (e-/A2).
        dosePerFrame : float, optional
            Expected dose per frame (e-/A2).
        sRateAngsPixTol : float
            Tolerance for sampling-rate comparisons.
        checkHeaderApix : bool
            If ``True`` (default) also validate the sampling rate in the
            file header (MRC voxel size). Set to ``False`` for protocols
            that only produce metadata without modifying the binary file.
        """
        if not isinstance(img, Image):
            self.fail(f"Expected Image, got {type(img)}.")

        # Object ID
        if imageId is not None:
            self.assertEqual(img.getObjId(), imageId,
                             msg=f"Image ID mismatch: "
                                 f"expected {imageId}, got {img.getObjId()}")

        # File name
        if imageName is not None:
            self.assertEqual(img.getFileName(), imageName,
                             msg=f"Image filename mismatch: "
                                 f"expected {imageName}, "
                                 f"got {img.getFileName()}")

        # Sampling rate
        if samplingRate is not None:
            self.assertAlmostEqual(img.getSamplingRate(), samplingRate,
                                   delta=sRateAngsPixTol,
                                   msg=f"Sampling rate mismatch: "
                                       f"expected {samplingRate}, "
                                       f"got {img.getSamplingRate()}")
            if checkHeaderApix:
                self.checkHeaderSRate(img, expectedSRate=samplingRate,
                                      sRateAngsPixTol=sRateAngsPixTol)

        # Dimensions
        if dim is not None:
            actualDim = img.getDim()
            self.assertIsNotNone(actualDim,
                                 msg=f"Image dimensions are None for "
                                     f"{img.getFileName()}.")
            self.assertEqual(actualDim, dim,
                             msg=f"Dimension mismatch: "
                                 f"expected {dim}, got {actualDim}")

        # CTF presence
        if hasCTF is not None:
            if hasCTF:
                self.assertTrue(img.hasCTF(),
                                msg="Image should have a CTF model but does not.")
            else:
                self.assertFalse(img.hasCTF(),
                                 msg="Image should not have a CTF model but does.")

        # Transform
        if transformShifts is not None:
            self.assertTrue(img.hasTransform(),
                            msg="Image should have a transform but does not.")
            self.checkTransform(img.getTransform(), shifts=transformShifts)

        # Origin
        if origin is not None:
            self.assertTrue(img.hasOrigin(),
                            msg="Image should have an origin but does not.")
            self.checkTransform(img.getOrigin(), shifts=origin)

        # Acquisition
        if img.hasAcquisition():
            self.checkAcquisition(img.getAcquisition(),
                                  voltage=voltage,
                                  sphericalAberration=sphericalAberration,
                                  amplitudeContrast=amplitudeContrast,
                                  magnification=magnification,
                                  doseInitial=doseInitial,
                                  dosePerFrame=dosePerFrame)

    def checkMicrograph(self,
                        mic: Micrograph,
                        micName: Optional[str] = None,
                        micId: Optional[int] = None,
                        samplingRate: Optional[float] = None,
                        dim: Optional[Tuple[int, int, int]] = None,
                        hasCTF: Optional[bool] = None,
                        transformShifts: Optional[Tuple[float, float, float]] = None,
                        origin: Optional[Tuple[float, float, float]] = None,
                        voltage: Optional[float] = None,
                        sphericalAberration: Optional[float] = None,
                        amplitudeContrast: Optional[float] = None,
                        magnification: Optional[float] = None,
                        doseInitial: Optional[float] = None,
                        dosePerFrame: Optional[float] = None,
                        sRateAngsPixTol: float = 0.01) -> None:
        """Validate a Micrograph object.

        See ``checkImage`` for the parameter documentation of inherited
        Image attributes.
        """
        if not isinstance(mic, Micrograph):
            self.fail(f"Expected Micrograph, got {type(mic)}.")

        # Micrograph-specific name
        if micName is not None:
            self.assertEqual(mic.getMicName(), micName,
                             msg=f"Micrograph name mismatch: "
                                 f"expected {micName}, got {mic.getMicName()}")

        self.checkImage(mic, imageId=micId,
                        samplingRate=samplingRate, dim=dim, hasCTF=hasCTF,
                        transformShifts=transformShifts, origin=origin,
                        voltage=voltage,
                        sphericalAberration=sphericalAberration,
                        amplitudeContrast=amplitudeContrast,
                        magnification=magnification,
                        doseInitial=doseInitial,
                        dosePerFrame=dosePerFrame,
                        sRateAngsPixTol=sRateAngsPixTol)

    def checkMovie(self,
                   mov: Movie,
                   movieId: Optional[int] = None,
                   movieName: Optional[str] = None,
                   samplingRate: Optional[float] = None,
                   voltage: Optional[float] = None,
                   dim: Optional[Tuple[int, int, int]] = None,
                   hasCTF: Optional[bool] = None,
                   transformShifts: Optional[Tuple[float, float, float]] = None,
                   origin: Optional[Tuple[float, float, float]] = None,
                   sphericalAberration: Optional[float] = None,
                   amplitudeContrast: Optional[float] = None,
                   magnification: Optional[float] = None,
                   doseInitial: Optional[float] = None,
                   dosePerFrame: Optional[float] = None,
                   framesRange: Optional[Tuple[int, int, int]] = None,
                   numFrames: Optional[int] = None,
                   alignment: Optional[MovieAlignment] = None,
                   alignmentFirst: Optional[int] = None,
                   alignmentLast: Optional[int] = None,
                   alignmentXShifts: Optional[List[float]] = None,
                   alignmentYShifts: Optional[List[float]] = None,
                   sRateAngsPixTol: Optional[float] = None) -> None:# 0.01) -> None:
        """Validate a Movie object.

        Parameters
        ----------
        mov : Movie
            Movie object to validate.
        movieId : int, optional
            Expected movie ID.
        movieName : str, optional
            Expected movie name.
        samplingRate : float, optional
            Expected sampling rate (A/pix).
        voltage : float, optional
            Expected microscope voltage (kV).
        dim : tuple of int, optional
            Expected (x, y, z) dimensions.
        framesRange : tuple of int, optional
            Expected (firstFrame, lastFrame, firstFrameIndex).
        numFrames : int, optional
            Expected number of frames.
        alignment : MovieAlignment, optional
            Optional alignment object for ``checkMovieAlignment``.
        alignmentFirst, alignmentLast : int
            Expected alignment frame range (used when *alignment* is given).
        alignmentXShifts, alignmentYShifts : list of float, optional
            Expected alignment shifts (used when *alignment* is given).
        sRateAngsPixTol : float
            Tolerance for sampling-rate comparisons.
        """
        if not isinstance(mov, Movie):
            self.fail(f"Expected Movie, got {type(mov)}.")

        # Check inherited Micrograph properties
        self.checkMicrograph(mov, micId=movieId, micName=movieName,
                             samplingRate=samplingRate, voltage=voltage,
                             hasCTF=hasCTF,
                             transformShifts=transformShifts, origin=origin,
                             sphericalAberration=sphericalAberration,
                             amplitudeContrast=amplitudeContrast,
                             magnification=magnification,
                             doseInitial=doseInitial,
                             dosePerFrame=dosePerFrame,
                             sRateAngsPixTol=sRateAngsPixTol)

        if dim is not None:
            self.assertEqual(mov.getDimensions(), dim,
                             msg=f"Dimensions mismatch: "
                                 f"expected {dim}, got {mov.getDimensions()}")

        # FramesRange
        fr = mov.getFramesRange()
        self.assertIsNotNone(fr, "Movie should have a FramesRange.")
        if framesRange is not None:
            first, last, firstIdx = framesRange
            self.assertEqual(fr.getFirstFrame(), first,
                             msg=f"First frame mismatch: "
                                 f"expected {first}, got {fr.getFirstFrame()}")
            self.assertEqual(fr.getLastFrame(), last,
                             msg=f"Last frame mismatch: "
                                 f"expected {last}, got {fr.getLastFrame()}")
            self.assertEqual(fr.getFirstFrameIndex(), firstIdx,
                             msg=f"First frame index mismatch: "
                                 f"expected {firstIdx}, "
                                 f"got {fr.getFirstFrameIndex()}")

        # Number of frames
        if numFrames is not None:
            self.assertEqual(mov.getNumberOfFrames(), numFrames,
                             msg=f"Number of frames mismatch: "
                                 f"expected {numFrames}, "
                                 f"got {mov.getNumberOfFrames()}")

        # Movie alignment
        if alignment is not None:
            self.checkMovieAlignment(alignment,
                                     expectedFirst=alignmentFirst,
                                     expectedLast=alignmentLast,
                                     expectedXShifts=alignmentXShifts,
                                     expectedYShifts=alignmentYShifts)

    def checkParticle(self,
                      particle: Particle,
                      classId: Optional[int] = None,
                      micId: Optional[int] = None,
                      particleId: Optional[int] = None,
                      imageName: Optional[str] = None,
                      samplingRate: Optional[float] = None,
                      dim: Optional[Tuple[int, int, int]] = None,
                      hasCTF: Optional[bool] = None,
                      transformShifts: Optional[Tuple[float, float, float]] = None,
                      origin: Optional[Tuple[float, float, float]] = None,
                      voltage: Optional[float] = None,
                      sphericalAberration: Optional[float] = None,
                      amplitudeContrast: Optional[float] = None,
                      magnification: Optional[float] = None,
                      doseInitial: Optional[float] = None,
                      dosePerFrame: Optional[float] = None,
                      sRateAngsPixTol: float = 0.01,
                      corExpectedX: Optional[int] = None,
                      corExpectedY: Optional[int] = None,
                      corExpectedMicId: Optional[int] = None,
                      checkHeaderApix: Optional[bool] = None) -> None:
        """Validate a Particle object.

        In addition to the ``checkImage`` parameters documented above,
        the following Particle-specific parameters are available:
        """
        if not isinstance(particle, Particle):
            self.fail(f"Expected Particle, got {type(particle)}.")

        # Particle-specific attributes
        if classId is not None:
            self.assertEqual(particle.getClassId(), classId,
                             msg=f"Class ID mismatch: "
                                 f"expected {classId}, "
                                 f"got {particle.getClassId()}")
        if micId is not None:
            self.assertEqual(particle.getMicId(), micId,
                             msg=f"Micrograph ID mismatch: "
                                 f"expected {micId}, got {particle.getMicId()}")

        # Inherited Image attributes
        self.checkImage(particle, imageName=imageName,
                        samplingRate=samplingRate, dim=dim, hasCTF=hasCTF,
                        transformShifts=transformShifts, origin=origin,
                        voltage=voltage,
                        sphericalAberration=sphericalAberration,
                        amplitudeContrast=amplitudeContrast,
                        magnification=magnification,
                        doseInitial=doseInitial,
                        dosePerFrame=dosePerFrame,
                        sRateAngsPixTol=sRateAngsPixTol,
                        checkHeaderApix=checkHeaderApix)

        # Associated coordinate
        coord = particle.getCoordinate()
        if any(x is not None for x in [corExpectedX, corExpectedY,
                                       corExpectedMicId]):
            self.checkCoordinate(coord, expectedX=corExpectedX,
                                 expectedY=corExpectedY,
                                 expectedMicId=corExpectedMicId)

    def checkVolume(self,
                    vol: Volume,
                    expectedSRate: Optional[float] = None,
                    expectedBoxSize: Optional[int] = None,
                    hasCTF: Optional[bool] = None,
                    hasHalves: Optional[bool] = None,
                    expectedOriginShifts: Optional[Union[List[float],
                                                         Tuple[float, ...]]] = None,
                    volumeId: Optional[int] = None,
                    volumeName: Optional[str] = None,
                    transformShifts: Optional[Tuple[float, float, float]] = None,
                    origin: Optional[Tuple[float, float, float]] = None,
                    voltage: Optional[float] = None,
                    sphericalAberration: Optional[float] = None,
                    amplitudeContrast: Optional[float] = None,
                    magnification: Optional[float] = None,
                    doseInitial: Optional[float] = None,
                    dosePerFrame: Optional[float] = None,
                    sRateAngsPixTol: Optional[float] = None, # 0.01,
                    checkHeaderApix: Optional[bool] = None) -> None: #True) -> None:
        """Validate a Volume object.

        Parameters
        ----------
        vol : Volume
            Volume object to validate.
        expectedSRate : float, optional
            Expected sampling rate (A/pix).
        expectedBoxSize : int, optional
            Expected cubic box size (pixels). If provided, dimensions
            are checked as (boxSize, boxSize, boxSize).
        hasCTF : bool, optional
            ``True``/``False`` to enforce CTF presence/absence.
        hasHalves : bool, optional
            ``True``/``False`` to enforce presence/absence of half-maps.
        expectedOriginShifts : list or tuple of float, optional
            Expected (x, y, z) origin shifts in A.
        sRateAngsPixTol : float
            Tolerance for sampling-rate comparisons.
        checkHeaderApix : bool
            If ``True`` also validate the sampling rate in the file header.
        """
        if not isinstance(vol, Volume):
            self.fail(f"Expected Volume, got {type(vol)}.")

        # Build dim tuple from box size
        dim = (expectedBoxSize, expectedBoxSize, expectedBoxSize) \
            if expectedBoxSize is not None else None

        # Origin tuple
        origin = tuple(expectedOriginShifts) \
            if expectedOriginShifts is not None else None

        self.assertEqual(vol.getClassId(), volumeId,
                         msg=f"Volume ID mismatch: "
                             f"expected {volumeId}, got {vol.getClassId()}")


        # Check inherited Image properties
        self.checkImage(vol, imageName=volumeName,
                        samplingRate=expectedSRate, dim=dim,
                        transformShifts=transformShifts,
                        hasCTF=hasCTF, origin=origin,
                        voltage=voltage,
                        sphericalAberration=sphericalAberration,
                        amplitudeContrast=amplitudeContrast,
                        magnification=magnification,
                        doseInitial=doseInitial,
                        dosePerFrame=dosePerFrame,
                        sRateAngsPixTol=sRateAngsPixTol,
                        checkHeaderApix=checkHeaderApix)

        # File existence
        fn = vol.getFileName()
        self.assertTrue(exists(fn),
                        msg=f"Volume file does not exist: {fn}")

        # Half-maps
        if hasHalves is not None:
            if hasHalves:
                self.assertTrue(vol.hasHalfMaps(),
                                msg="Volume should have half-maps registered "
                                    "but does not.")
                half1, half2 = vol.getHalfMaps().split(',')
                self.assertTrue(exists(half1),
                                msg=f"Volume first half does not exist: "
                                    f"{half1}")
                self.assertTrue(exists(half2),
                                msg=f"Volume second half does not exist: "
                                    f"{half2}")

    def checkCoordinate(self,
                        coord: Coordinate,
                        expectedX: Optional[int] = None,
                        expectedY: Optional[int] = None,
                        expectedMicId: Optional[int] = None) -> None:
        """Validate a Coordinate object.

        Parameters
        ----------
        coord : Coordinate
            Coordinate object to validate.
        expectedX : int, optional
            Expected X coordinate (pixels).
        expectedY : int, optional
            Expected Y coordinate (pixels).
        expectedMicId : int, optional
            Expected micrograph ID the coordinate belongs to.
        """
        if not isinstance(coord, Coordinate):
            self.fail(f"Expected Coordinate, got {type(coord)}.")
        if expectedX is not None:
            self.assertEqual(coord.getX(), expectedX,
                             msg=f"X coordinate mismatch: "
                                 f"expected {expectedX}, got {coord.getX()}")
        if expectedY is not None:
            self.assertEqual(coord.getY(), expectedY,
                             msg=f"Y coordinate mismatch: "
                                 f"expected {expectedY}, got {coord.getY()}")
        if expectedMicId is not None:
            self.assertEqual(coord.getMicId(), expectedMicId,
                             msg=f"Micrograph ID mismatch: "
                                 f"expected {expectedMicId}, "
                                 f"got {coord.getMicId()}")

    def checkAtomStruct(self,
                        atomStruct: AtomStruct,
                        hasVolume: bool = False,
                        pseudoatoms: bool = False) -> None:
        """Validate an AtomStruct (e.g. PDB file) object.

        Parameters
        ----------
        atomStruct : AtomStruct
            Object to validate.
        hasVolume : bool
            Whether an associated volume is expected.
        pseudoatoms : bool
            Whether the structure uses pseudo-atoms.
        """
        if not isinstance(atomStruct, AtomStruct):
            self.fail(f"Expected AtomStruct, got {type(atomStruct)}.")
        self.assertTrue(exists(atomStruct.getFileName()),
                        msg=f"AtomStruct file does not exist: "
                            f"{atomStruct.getFileName()}")
        self.assertEqual(atomStruct.hasVolume(), hasVolume,
                         msg=f"hasVolume mismatch: "
                             f"expected {hasVolume}, "
                             f"got {atomStruct.hasVolume()}")
        self.assertEqual(atomStruct.getPseudoAtoms(), pseudoatoms,
                         msg=f"pseudoatoms mismatch: "
                             f"expected {pseudoatoms}, "
                             f"got {atomStruct.getPseudoAtoms()}")

    def checkMask(self,
                  mask: Mask,
                  classId: Optional[int] = None,
                  micId: Optional[int] = None,
                  imageId: Optional[int] = None,
                  imageName: Optional[str] = None,
                  samplingRate: Optional[float] = None,
                  dim: Optional[Tuple[int, int, int]] = None,
                  hasCTF: Optional[bool] = None,
                  transformShifts: Optional[Tuple[float, float, float]] = None,
                  origin: Optional[Tuple[float, float, float]] = None,
                  voltage: Optional[float] = None,
                  sphericalAberration: Optional[float] = None,
                  amplitudeContrast: Optional[float] = None,
                  magnification: Optional[float] = None,
                  doseInitial: Optional[float] = None,
                  dosePerFrame: Optional[float] = None,
                  sRateAngsPixTol: float = 0.01,
                  corExpectedX: Optional[int] = None,
                  corExpectedY: Optional[int] = None,
                  corExpectedMicId: Optional[int] = None) -> None:
        """Validate a Mask object (delegates to ``checkParticle``)."""
        if not isinstance(mask, Mask):
            self.fail(f"Expected Mask, got {type(mask)}.")
        self.checkParticle(mask,
                           classId=classId, micId=micId,
                           particleId=imageId, imageName=imageName,
                           samplingRate=samplingRate, dim=dim, hasCTF=hasCTF,
                           transformShifts=transformShifts, origin=origin,
                           voltage=voltage,
                           sphericalAberration=sphericalAberration,
                           amplitudeContrast=amplitudeContrast,
                           magnification=magnification,
                           doseInitial=doseInitial,
                           dosePerFrame=dosePerFrame,
                           sRateAngsPixTol=sRateAngsPixTol,
                           corExpectedX=corExpectedX,
                           corExpectedY=corExpectedY,
                           corExpectedMicId=corExpectedMicId)

    def checkVolumeMask(self,
                        mask3d: VolumeMask,
                        expectedSRate: Optional[float] = None,
                        expectedBoxSize: Optional[int] = None,
                        hasCTF: Optional[bool] = None,
                        hasHalves: Optional[bool] = None,
                        expectedOriginShifts: Optional[Union[List[float],
                                                             Tuple[float, ...]]] = None,
                        volumeId: Optional[int] = None,
                        volumeName: Optional[str] = None,
                        transformShifts: Optional[Tuple[float, float, float]] = None,
                        origin: Optional[Tuple[float, float, float]] = None,
                        voltage: Optional[float] = None,
                        sphericalAberration: Optional[float] = None,
                        amplitudeContrast: Optional[float] = None,
                        magnification: Optional[float] = None,
                        doseInitial: Optional[float] = None,
                        dosePerFrame: Optional[float] = None,
                        checkHeaderApix: Optional[bool] = None,
                        sRateAngsPixTol: float = 0.01) -> None:
        """Validate a VolumeMask object (delegates to checkVolume)."""
        if not isinstance(mask3d, VolumeMask):
            self.fail(f"Expected VolumeMask, got {type(mask3d)}.")
        self.checkVolume(mask3d, expectedSRate=expectedSRate,
                         expectedBoxSize=expectedBoxSize,
                         hasCTF=hasCTF, hasHalves=hasHalves,
                         expectedOriginShifts=expectedOriginShifts,
                         volumeId=volumeId, volumeName=volumeName,
                         transformShifts=transformShifts, origin=origin,
                         voltage=voltage,
                         sphericalAberration=sphericalAberration,
                         amplitudeContrast=amplitudeContrast,
                         magnification=magnification,
                         doseInitial=doseInitial,
                         dosePerFrame=dosePerFrame,
                         checkHeaderApix=checkHeaderApix,
                         sRateAngsPixTol=sRateAngsPixTol)

    # =========================================================================
    # 3. SET CHECKS
    # =========================================================================
    def checkImageSet(self,
                      inImageSet: SetOfImages,
                      expectedSetSize: int,
                      expectedSRate: Optional[float] = None,
                      hasCtf: bool = False,
                      testAcqObj: Optional[Acquisition] = None,
                      streamState: Optional[int] = None,
                      sRateAngsPixTol: float = 0.01) -> None:
        self._checkImageSet(inImageSet, expectedSetSize, expectedSRate,
                            hasCtf, testAcqObj, streamState, sRateAngsPixTol)
        
        for img in inImageSet:
            imgId = img.getObjId()
            print(f'---> Checking volume (objId={imgId})')
            self.checkImage(img, expectedSRate=expectedSRate,
                             hasCTF=hasCtf,
                             sRateAngsPixTol=sRateAngsPixTol)

    def _checkImageSet(self,
                      inImageSet: SetOfImages,
                      expectedSetSize: int,
                      expectedSRate: Optional[float] = None,
                      hasCtf: bool = False,
                      testAcqObj: Optional[Acquisition] = None,
                      streamState: Optional[int] = None,
                      sRateAngsPixTol: float = 0.01) -> None:
        """Validate a SetOfImages (set-level properties only).
        Parameters
        ----------
        inImageSet : SetOfImages
            Set to validate.
        expectedSetSize : int
            Expected number of items.
        expectedSRate : float, optional
            Expected sampling rate (A/pix).
        hasCtf : bool
            Expected CTF presence flag on the set.
        testAcqObj : Acquisition, optional
            Expected acquisition values (set-level).
        streamState : int, optional
            Expected stream state.
        sRateAngsPixTol : float
            Tolerance for sampling-rate comparisons.
        """
        if not isinstance(inImageSet, SetOfImages):
            self.fail(f"Expected SetOfImages, got {type(inImageSet)}.")
        self.checkSetGeneralProps(inImageSet,
                                  expectedSetSize=expectedSetSize,
                                  expectedSRate=expectedSRate,
                                  streamState=streamState,
                                  sRateAngsPixTol=sRateAngsPixTol)
        self.assertEqual(inImageSet.hasCTF(), hasCtf,
                         msg=f"SetOfImages CTF flag mismatch: "
                             f"expected {hasCtf}, got {inImageSet.hasCTF()}")
        if testAcqObj is not None and inImageSet.hasAcquisition():
            acq = inImageSet.getAcquisition()
            self.checkAcquisition(acq,
                                  voltage=testAcqObj.getVoltage(),
                                  sphericalAberration=testAcqObj.getSphericalAberration(),
                                  amplitudeContrast=testAcqObj.getAmplitudeContrast(),
                                  magnification=testAcqObj.getMagnification())

    def checkVolumeSet(self,
                       inVolumeSet: SetOfVolumes,
                       expectedSetSize: int,
                       expectedSRate: Optional[float] = None,
                       expectedBoxSize: Optional[int] = None,
                       hasCtf: bool = False,
                       hasHalves: bool = False,
                       testAcqObj: Optional[Acquisition] = None,
                       checkHeaderApix: bool = True,
                       streamState: Optional[int] = None,
                       expectedOriginShifts: Optional[Union[List[float],
                                                                Tuple[float, ...]]] = None,
                       voltage: Optional[float] = None,
                       sphericalAberration: Optional[float] = None,
                       amplitudeContrast: Optional[float] = None,
                       magnification: Optional[float] = None,
                       doseInitial: Optional[float] = None,
                       dosePerFrame: Optional[float] = None,
                       sRateAngsPixTol: float = 0.01) -> None:
        """Validate a SetOfVolumes: set-level props + each item via checkVolume.

        Parameters
        ----------
        inVolumeSet : SetOfVolumes
            Set to validate.
        expectedSetSize : int
            Expected number of volumes.
        expectedSRate : float
            Expected sampling rate (A/pix).
        expectedBoxSize : int, optional
            Expected cubic box size (pixels). Passed to ``checkVolume``.
        hasCtf : bool
            Expected CTF flag on the set.
        hasHalves : bool
            Expected half-maps flag on the set.
        testAcqObj : Acquisition, optional
            Expected acquisition values (set-level).
        checkHeaderApix : bool
            Whether to check the voxel size in the MRC header.
        streamState : int, optional
            Expected stream state.
        expectedOriginShifts : list or tuple of float, optional
            Expected origin shifts for each volume.
        sRateAngsPixTol : float
            Tolerance for sampling-rate comparisons.
        """
        if not isinstance(inVolumeSet, SetOfVolumes):
            self.fail(f"Expected SetOfVolumes, got {type(inVolumeSet)}.")
        self._checkImageSet(inVolumeSet,
                            expectedSetSize=expectedSetSize,
                            expectedSRate=expectedSRate,
                            hasCtf=hasCtf,
                            testAcqObj=testAcqObj,
                            streamState=streamState,
                            sRateAngsPixTol=sRateAngsPixTol)
        self.assertEqual(inVolumeSet.hasCTF(), hasCtf,
                         msg=f"SetOfVolumes CTF flag mismatch: "
                             f"expected {hasCtf}, got {inVolumeSet.hasCTF()}")
        if testAcqObj is not None and inVolumeSet.hasAcquisition():
            self.checkAcquisition(inVolumeSet.getAcquisition(),
                                  voltage=testAcqObj.getVoltage(),
                                  sphericalAberration=testAcqObj.getSphericalAberration(),
                                  amplitudeContrast=testAcqObj.getAmplitudeContrast(),
                                  magnification=testAcqObj.getMagnification())
        for vol in inVolumeSet:
            print(f'---> Checking volume (objId={vol.getObjId()})')
            self.checkVolume(vol,
                             expectedSRate=expectedSRate,
                             expectedBoxSize=expectedBoxSize,
                             hasCTF=hasCtf,
                             hasHalves=hasHalves,
                             expectedOriginShifts=expectedOriginShifts,
                             voltage=voltage,
                             sphericalAberration=sphericalAberration,
                             amplitudeContrast=amplitudeContrast,
                             magnification=magnification,
                             doseInitial=doseInitial,
                             dosePerFrame=dosePerFrame,
                             sRateAngsPixTol=sRateAngsPixTol,
                             checkHeaderApix=checkHeaderApix)

    def checkSetOfMicrographs(self,
                              inMicSet: SetOfMicrographs,
                              expectedSetSize: int,
                              expectedSRate: float = None,
                              hasCtf: bool = False,
                              testAcqObj: Optional[Acquisition] = None,
                              streamState: Optional[int] = None,
                              sRateAngsPixTol: float = 0.01) -> None:
        """Validate a SetOfMicrographs: set-level props + each item.

        Parameters
        ----------
        inMicSet : SetOfMicrographs
            Set to validate.
        expectedSetSize : int
            Expected number of micrographs.
        expectedSRate : float
            Expected sampling rate (A/pix).
        hasCtf : bool
            Expected CTF flag.
        testAcqObj : Acquisition, optional
            Set-level expected acquisition.
        streamState : int, optional
            Expected stream state.
        sRateAngsPixTol : float
            Tolerance for sampling-rate comparisons.
        """
        if not isinstance(inMicSet, SetOfMicrographs):
            self.fail(f"Expected SetOfMicrographs, got {type(inMicSet)}.")
        self.checkSetGeneralProps(inMicSet,
                                  expectedSetSize=expectedSetSize,
                                  expectedSRate=expectedSRate,
                                  streamState=streamState,
                                  sRateAngsPixTol=sRateAngsPixTol)
        self.assertEqual(inMicSet.hasCTF(), hasCtf,
                         msg=f"SetOfMicrographs CTF flag mismatch: "
                             f"expected {hasCtf}, got {inMicSet.hasCTF()}")
        if testAcqObj is not None and inMicSet.hasAcquisition():
            self.checkAcquisition(inMicSet.getAcquisition(),
                                  voltage=testAcqObj.getVoltage(),
                                  sphericalAberration=testAcqObj.getSphericalAberration(),
                                  amplitudeContrast=testAcqObj.getAmplitudeContrast(),
                                  magnification=testAcqObj.getMagnification())
        for mic in inMicSet:
            print(f'---> Checking micrograph (objId={mic.getObjId()})')
            self.checkMicrograph(mic, samplingRate=expectedSRate,
                                 sRateAngsPixTol=sRateAngsPixTol)

    def checkSetOfMovies(self,
                         movieSet: SetOfMovies,
                         expectedSetSize: int,
                         movieIds: Optional[List[int]] = None,
                         movieNames: Optional[List[str]] = None,
                         dim: Optional[Tuple[int, int, int]] = None,
                         expectedSRate: Optional[float] = None,
                         expectedGain: Optional[str] = None,
                         expectedDark: Optional[str] = None,
                         streamState: Optional[int] = None,
                         sRateAngsPixTol: float = 0.01,
                         voltage: Optional[float] = None,
                         sphericalAberration: Optional[float] = None,
                         amplitudeContrast: Optional[float] = None,
                         magnification: Optional[float] = None,
                         doseInitial: Optional[float] = None,
                         dosePerFrame: Optional[float] = None) -> None:
        """Validate a SetOfMovies: set-level props + each item.

        Parameters
        ----------
        movieSet : SetOfMovies
            Set to validate.
        expectedSetSize : int
            Expected number of movies.
        movieIds : List[int], optional
            Expected movie IDs.
        movieNames : List[str], optional
            Expected movie names.
        expectedSRate : float, optional
            Expected sampling rate (A/pix).
        expectedGain : str, optional
            Expected gain reference filename.
        expectedDark : str, optional
            Expected dark reference filename.
        streamState : int, optional
            Expected stream state.
        sRateAngsPixTol : float
            Tolerance for sampling-rate comparisons.
        """
        if not isinstance(movieSet, SetOfMovies):
            self.fail(f"Expected SetOfMovies, got {type(movieSet)}.")
        self.checkSetGeneralProps(movieSet,
                                  expectedSetSize=expectedSetSize,
                                  expectedSRate=expectedSRate,
                                  streamState=streamState,
                                  sRateAngsPixTol=sRateAngsPixTol)
        if expectedGain is not None:
            self.assertEqual(movieSet.getGain(), expectedGain,
                             msg=f"Gain reference mismatch: "
                                 f"expected {expectedGain}, "
                                 f"got {movieSet.getGain()}")
        if expectedDark is not None:
            self.assertEqual(movieSet.getDark(), expectedDark,
                             msg=f"Dark reference mismatch: "
                                 f"expected {expectedDark}, "
                                 f"got {movieSet.getDark()}")
        for mov in movieSet:
            print(f'---> Checking movie (objId={mov.getObjId()})')
            self.checkMovie(mov, dim=dim, samplingRate=expectedSRate,
                            sRateAngsPixTol=sRateAngsPixTol,
                            voltage=voltage,
                            sphericalAberration=sphericalAberration,
                            amplitudeContrast=amplitudeContrast,
                            magnification=magnification,
                            doseInitial=doseInitial,
                            dosePerFrame=dosePerFrame)

    def checkSetOfParticles(self,
                            inParticleSet: SetOfParticles,
                            expectedSetSize: int,
                            expectedSRate: float,
                            isSubparticle: bool = False,
                            hasCtf: bool = False,
                            testAcqObj: Optional[Acquisition] = None,
                            streamState: Optional[int] = None,
                            sRateAngsPixTol: float = 0.01) -> None:
        """Validate a SetOfParticles: set-level props + each item.

        Parameters
        ----------
        inParticleSet : SetOfParticles
            Set to validate.
        expectedSetSize : int
            Expected number of particles.
        expectedSRate : float
            Expected sampling rate (A/pix).
        isSubparticle : bool
            Expected subparticle flag on the set.
        hasCtf : bool
            Expected CTF flag.
        testAcqObj : Acquisition, optional
            Set-level expected acquisition.
        streamState : int, optional
            Expected stream state.
        sRateAngsPixTol : float
            Tolerance for sampling-rate comparisons.
        """
        if not isinstance(inParticleSet, SetOfParticles):
            self.fail(f"Expected SetOfParticles, got {type(inParticleSet)}.")
        self.checkSetGeneralProps(inParticleSet,
                                  expectedSetSize=expectedSetSize,
                                  expectedSRate=expectedSRate,
                                  streamState=streamState,
                                  sRateAngsPixTol=sRateAngsPixTol)
        self.assertEqual(inParticleSet.getIsSubparticles(), isSubparticle,
                         msg=f"IsSubparticle flag mismatch: "
                             f"expected {isSubparticle}, "
                             f"got {inParticleSet.getIsSubparticles()}")
        self.assertEqual(inParticleSet.hasCTF(), hasCtf,
                         msg=f"SetOfParticles CTF flag mismatch: "
                             f"expected {hasCtf}, "
                             f"got {inParticleSet.hasCTF()}")
        if testAcqObj is not None and inParticleSet.hasAcquisition():
            self.checkAcquisition(inParticleSet.getAcquisition(),
                                  voltage=testAcqObj.getVoltage(),
                                  sphericalAberration=testAcqObj.getSphericalAberration(),
                                  amplitudeContrast=testAcqObj.getAmplitudeContrast(),
                                  magnification=testAcqObj.getMagnification())
        for particle in inParticleSet:
            print(f'---> Checking particle (objId={particle.getObjId()})')
            self.checkParticle(particle,
                               samplingRate=expectedSRate,
                               sRateAngsPixTol=sRateAngsPixTol)

    def checkSetOfCoordinates(self,
                              coordSet: SetOfCoordinates,
                              expectedSize: int,
                              expectedBoxSize: Optional[int] = None,
                              expectedMicSet: Optional[SetOfMicrographs] = None,
                              expectedCoordList: Optional[List[dict]] = None
                              ) -> None:
        """Validate a SetOfCoordinates: set-level props + each item.

        Parameters
        ----------
        coordSet : SetOfCoordinates
            Set to validate.
        expectedSize : int
            Expected number of coordinates.
        expectedBoxSize : int, optional
            Expected box size (pixels).
        expectedMicSet : SetOfMicrographs, optional
            If provided, the micrograph set linked to the coordinates
            will be verified.
        expectedCoordList : list of dict, optional
            Optional per-item expected values. Each dict may have keys
            ``expectedX``, ``expectedY``, ``expectedMicId``.
        """
        if not isinstance(coordSet, SetOfCoordinates):
            self.fail(f"Expected SetOfCoordinates, got {type(coordSet)}.")
        self.assertSetSize(coordSet, expectedSize)
        if expectedBoxSize is not None:
            self.assertEqual(coordSet.getBoxSize(), expectedBoxSize,
                             msg=f"Box size mismatch: "
                                 f"expected {expectedBoxSize}, "
                                 f"got {coordSet.getBoxSize()}")
        if expectedMicSet is not None:
            self.assertEqual(coordSet.getMicrographs().getObjId(),
                             expectedMicSet.getObjId(),
                             msg="The micrograph set linked to the "
                                 "coordinates does not match.")
        if expectedCoordList:
            self.assertEqual(len(expectedCoordList), expectedSize,
                             msg=f"expectedCoordList length "
                                 f"({len(expectedCoordList)}) must match "
                                 f"expectedSize ({expectedSize}).")
            for coord, expected in zip(coordSet.iterCoordinates(),
                                       expectedCoordList):
                self.checkCoordinate(coord,
                                     expectedX=expected.get('expectedX'),
                                     expectedY=expected.get('expectedY'),
                                     expectedMicId=expected.get('expectedMicId'))
        else:
            for coord in coordSet.iterCoordinates():
                self.checkCoordinate(coord)

    def checkSetOfClasses(self,
                          classesSet: SetOfClasses,
                          expectedSize: int,
                          hasRepresentatives: bool = True) -> None:
        """Validate a SetOfClasses: set-level props + representatives.

        Parameters
        ----------
        classesSet : SetOfClasses
            Set to validate.
        expectedSize : int
            Expected number of classes.
        hasRepresentatives : bool
            Whether representatives are expected.
        """
        if not isinstance(classesSet, SetOfClasses):
            self.fail(f"Expected SetOfClasses, got {type(classesSet)}.")
        self.assertSetSize(classesSet, expectedSize)
        self.assertEqual(classesSet.hasRepresentatives(), hasRepresentatives,
                         msg=f"hasRepresentatives mismatch: "
                             f"expected {hasRepresentatives}, "
                             f"got {classesSet.hasRepresentatives()}")
        if hasRepresentatives:
            repFiles = []
            for cls in classesSet:
                rep = cls.getRepresentative()
                repFn = rep.getFileName()
                repFiles.append(repFn)
                self.assertTrue(exists(repFn),
                                msg=f"Representative file does not exist: "
                                    f"{repFn}")
            self.assertEqual(len(set(repFiles)), expectedSize,
                             msg="At least one representative filename is "
                                 "repeated, which should not be possible.")

    def checkSetOfCTF(self,
                      ctfSet: SetOfCTF,
                      expectedSize: int,
                      expectedMicSet: Optional[SetOfMicrographs] = None,
                      expectedCtfList: Optional[List[dict]] = None) -> None:
        """Validate a SetOfCTF: set-level props + each CTF.

        Parameters
        ----------
        ctfSet : SetOfCTF
            Set to validate.
        expectedSize : int
            Expected number of CTF models.
        expectedMicSet : SetOfMicrographs, optional
            If provided, the micrograph set linked to the CTFs is verified.
        expectedCtfList : list of dict, optional
            Optional per-item expected values. Each dict may have keys
            ``defocusU``, ``defocusV``, ``defocusAngle``, ``resolution``,
            ``phaseShift``.
        """
        if not isinstance(ctfSet, SetOfCTF):
            self.fail(f"Expected SetOfCTF, got {type(ctfSet)}.")
        self.assertSetSize(ctfSet, expectedSize)
        micPtr = ctfSet.getMicrographs()
        self.assertIsNotNone(micPtr,
                             msg="SetOfCTF has no associated micrographs.")
        if expectedMicSet is not None:
            self.assertEqual(micPtr.getObjId(), expectedMicSet.getObjId(),
                             msg="The micrograph set linked to the CTFs "
                                 "does not match.")
        for ctf in ctfSet:
            self.assertNotEqual(ctf.getDefocusU(), -999,
                                msg=f"CTF id={ctf.getObjId()} has wrong "
                                    f"defocus marker (-999) -- parsing "
                                    f"likely failed.")
        if expectedCtfList:
            self.assertEqual(len(expectedCtfList), expectedSize,
                             msg=f"expectedCtfList length "
                                 f"({len(expectedCtfList)}) must match "
                                 f"expectedSize ({expectedSize}).")
            for ctf, expected in zip(ctfSet, expectedCtfList):
                self.checkCTF(ctf,
                              defocusU=expected['defocusU'],
                              defocusV=expected['defocusV'],
                              defocusAngle=expected['defocusAngle'],
                              resolution=expected.get('resolution'),
                              phaseShift=expected.get('phaseShift'))

    # =========================================================================
    # 4. STATIC HELPERS
    # =========================================================================

    @staticmethod
    def getMinAndMaxCoordValuesFromSet(inSet):
        """Return extreme coordinate values from a set.

        Works with ``SetOfCoordinates`` or ``SetOfParticles`` (in which
        case ``getCoordinates()`` is called internally).

        Returns
        -------
        numpy.ndarray
            ``[x_min, x_max, y_min, y_max]``.
        """
        if not isinstance(inSet, SetOfCoordinates):
            inSet = inSet.getCoordinates()
        dataDict = inSet.aggregate(['MAX'], '_micId', ['_x', '_y'])
        xcoords, ycoords = zip(*[(d['_x'], d['_y']) for d in dataDict])
        return np.array([min(xcoords), max(xcoords),
                         min(ycoords), max(ycoords)])
