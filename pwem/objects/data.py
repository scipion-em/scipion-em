# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This modules contains basic hierarchy
for EM data objects like: Image, SetOfImage and others
"""
import logging

logger = logging.getLogger(__name__)

import os
import json
import numpy as np

import pyworkflow.utils as pwutils
from pyworkflow.object import (Object, Float, Integer, String,
                               OrderedDict, CsvList, Boolean, Set, Pointer,
                               Scalar)
from pwem.constants import (NO_INDEX, ALIGN_NONE, ALIGN_2D, ALIGN_3D,
                            ALIGN_PROJ, ALIGNMENTS, LoopActions)


class EMObject(Object):
    """Base object for all EM classes"""

    def __str__(self):
        return self.getClassName()

    def getFiles(self):
        """ Get all filePaths """
        return None


class Acquisition(EMObject):
    """Acquisition information"""

    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._magnification = Float(kwargs.get('magnification', None))
        # Microscope voltage in kV
        self._voltage = Float(kwargs.get('voltage', None))
        # Spherical aberration in mm
        self._sphericalAberration = Float(kwargs.get('sphericalAberration',
                                                     None))
        self._amplitudeContrast = Float(kwargs.get('amplitudeContrast', None))
        self._doseInitial = Float(kwargs.get('doseInitial', 0))
        self._dosePerFrame = Float(kwargs.get('dosePerFrame', None))

        # Generic string to store information about optics groups
        self.opticsGroupInfo = String()

    def copyInfo(self, other):
        self.copy(other, copyId=False)

    def getMagnification(self):
        return self._magnification.get()

    def setMagnification(self, value):
        self._magnification.set(value)

    def getVoltage(self):
        return self._voltage.get()

    def setVoltage(self, value):
        self._voltage.set(value)

    def getSphericalAberration(self):
        return self._sphericalAberration.get()

    def setSphericalAberration(self, value):
        self._sphericalAberration.set(value)

    def getAmplitudeContrast(self):
        return self._amplitudeContrast.get()

    def setAmplitudeContrast(self, value):
        self._amplitudeContrast.set(value)

    def getDoseInitial(self):
        return self._doseInitial.get()

    def setDoseInitial(self, value):
        self._doseInitial.set(value)

    def getDosePerFrame(self):
        return self._dosePerFrame.get()

    def setDosePerFrame(self, value):
        self._dosePerFrame.set(value)

    def __str__(self):
        return "\n    mag=%s\n    volt= %s\n    Cs=%s\n    Q0=%s\n\n" % \
            (self._magnification.get(),
             self._voltage.get(),
             self._sphericalAberration.get(),
             self._amplitudeContrast.get())


class Transform(EMObject):
    """ This class will contain a transformation matrix
    that can be applied to 2D/3D objects like images and volumes.
    It should contain information about euler angles, translation(or shift)
    and mirroring.
    Shifts are stored in pixels as treated in extract coordinates, or assign angles,...
    """

    # Basic Transformation factory
    ROT_X_90_CLOCKWISE = 'rotX90c'
    ROT_Y_90_CLOCKWISE = 'rotY90c'
    ROT_Z_90_CLOCKWISE = 'rotZ90c'
    ROT_X_90_COUNTERCLOCKWISE = 'rotX90cc'
    ROT_Y_90_COUNTERCLOCKWISE = 'rotY90cc'
    ROT_Z_90_COUNTERCLOCKWISE = 'rotZ90cc'

    def __init__(self, matrix=None, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._matrix = Matrix()
        if matrix is not None:
            self.setMatrix(matrix)

    def getMatrix(self):
        return self._matrix.getMatrix()

    def getRotationMatrix(self):
        M = self.getMatrix()
        return M[:3, :3]

    def getEulerAngles(self):
        from pwem.convert.transformations import euler_from_matrix
        rotation = self.getRotationMatrix()
        return euler_from_matrix(rotation)

    def getShifts(self):
        M = self.getMatrix()
        return M[1, 4], M[2, 4], M[3, 4]

    def getMatrixAsList(self):
        """ Return the values of the Matrix as a list. """
        return self._matrix.getMatrix().flatten().tolist()

    def setMatrix(self, matrix):
        self._matrix.setMatrix(matrix)

    def __str__(self):
        return str(self._matrix)

    def scale(self, factor):
        m = self.getMatrix()
        m *= factor
        m[3, 3] = 1.

    def scaleShifts(self, factor):
        # By default Scipion uses a coordinate system associated with the volume rather than the projection
        m = self.getMatrix()
        m[0, 3] *= factor
        m[1, 3] *= factor
        m[2, 3] *= factor

    def invert(self):
        # Local import to avoid loop pwem --> data --> convert --> Plugin (at pwem)
        from pwem.convert.transformations import inverse_matrix

        self._matrix.setMatrix(inverse_matrix(self._matrix.getMatrix()))

        return self._matrix

    def getShifts(self):
        m = self.getMatrix()
        return m[0, 3], m[1, 3], m[2, 3]

    def setShifts(self, x, y, z):
        m = self.getMatrix()
        m[0, 3] = x
        m[1, 3] = y
        m[2, 3] = z

    def setShiftsTuple(self, shifts):
        self.setShifts(shifts[0], shifts[1], shifts[2])

    def composeTransform(self, matrix):
        """Apply a transformation matrix to the current matrix """
        new_matrix = np.matmul(matrix, self.getMatrix())
        # new_matrix = matrix * self.getMatrix()
        self._matrix.setMatrix(new_matrix)

    @classmethod
    def create(cls, type):
        if type == cls.ROT_X_90_CLOCKWISE:
            return Transform(matrix=np.array([
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, -1, 0, 0],
                [0, 0, 0, 1]]))
        elif type == cls.ROT_X_90_COUNTERCLOCKWISE:
            return Transform(matrix=np.array([
                [1, 0, 0, 0],
                [0, 0, -1, 0],
                [0, 1, 0, 0],
                [0, 0, 0, 1]]))
        elif type == cls.ROT_Y_90_CLOCKWISE:
            return Transform(matrix=np.array([
                [1, 0, -1, 0],
                [0, 1, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 0, 1]]))
        elif type == cls.ROT_Y_90_COUNTERCLOCKWISE:
            return Transform(matrix=np.array([
                [1, 0, 1, 0],
                [0, 1, 0, 0],
                [-1, 0, 0, 0],
                [0, 0, 0, 1]]))
        elif type == cls.ROT_Z_90_CLOCKWISE:
            return Transform(matrix=np.array([
                [0, 1, 0, 0],
                [-1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]]))
        elif type == cls.ROT_Z_90_COUNTERCLOCKWISE:
            return Transform(matrix=np.array([
                [0, -1, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]]))
        else:
            TRANSFORMATION_FACTORY_TYPES = [
                cls.ROT_X_90_CLOCKWISE,
                cls.ROT_Y_90_CLOCKWISE,
                cls.ROT_Z_90_CLOCKWISE,
                cls.ROT_X_90_COUNTERCLOCKWISE,
                cls.ROT_Y_90_COUNTERCLOCKWISE,
                cls.ROT_Z_90_COUNTERCLOCKWISE
            ]
            raise Exception('Introduced Transformation type is not recognized.\nAdmitted values are\n'
                            '%s' % ' '.join(TRANSFORMATION_FACTORY_TYPES))


class CTFModel(EMObject):
    """ Represents a generic CTF model. """

    DEFOCUS_V_MINIMUM_VALUE = 0.1
    DEFOCUS_RATIO_ERROR_VALUE = -1

    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._defocusU = Float(kwargs.get('defocusU', None))
        self._defocusV = Float(kwargs.get('defocusV', None))
        self._defocusAngle = Float(kwargs.get('defocusAngle', None))
        self._defocusRatio = Float()
        self._phaseShift = Float(kwargs['phaseShift']) if 'phaseShift' in kwargs else None
        self._psdFile = String()
        self._micObj = None
        self._resolution = Float()
        self._fitQuality = Float()

    def __str__(self):

        def strEx(value):
            return "None" if value is None else "%0.2f" % value

        ctfStr = "defU={}, defV={}, ast={}, " \
                 "psh={}, res={}, fit={}".format(
            strEx(self._defocusU.get()),
            strEx(self._defocusV.get()),
            strEx(self._defocusU.get(0) - self._defocusV.get(0)),
            strEx(self.getPhaseShift()),
            strEx(self._resolution.get()),
            strEx(self._fitQuality.get())
        )

        return ctfStr

    def getPhaseShift(self):
        if self._phaseShift is not None:
            return self._phaseShift.get()
        else:
            return None

    def setPhaseShift(self, value):
        self._phaseShift = Float(value)

    def hasPhaseShift(self):
        return False if self.getPhaseShift() is None else True

    def getResolution(self):
        return self._resolution.get()

    def hasResolution(self):
        return self._resolution.hasValue()

    def setResolution(self, value):
        self._resolution.set(value)

    def getFitQuality(self):
        return self._fitQuality.get()

    def setFitQuality(self, value):
        self._fitQuality.set(value)

    def getDefocusU(self):
        return self._defocusU.get()

    def setDefocusU(self, value):
        self._defocusU.set(value)

    def getDefocusV(self):
        return self._defocusV.get()

    def setDefocusV(self, value):
        self._defocusV.set(value)

    def getDefocusAngle(self):
        return self._defocusAngle.get()

    def setDefocusAngle(self, value):
        self._defocusAngle.set(value)

    def getDefocusRatio(self):
        return self._defocusRatio.get()

    def copyInfo(self, other):
        self.copyAttributes(other, '_defocusU', '_defocusV', '_defocusAngle',
                            '_defocusRatio', '_psdFile', '_micFile',
                            '_resolution', '_fitQuality')
        if other.hasPhaseShift():
            self.setPhaseShift(other.getPhaseShift())

    def getPsdFile(self):
        return self._psdFile.get()

    def setPsdFile(self, value):
        self._psdFile.set(value)

    def getMicrograph(self):
        self._micObj.copyObjId(self)
        return self._micObj

    def setMicrograph(self, mic):
        self._micObj = mic
        self.copyObjId(mic)

    def getDefocus(self):
        """ Returns defocusU, defocusV and defocusAngle. """
        return (self._defocusU.get(),
                self._defocusV.get(),
                self._defocusAngle.get())

    def setStandardDefocus(self, defocusU, defocusV, defocusAngle):
        """ Set defocus values following emx conventions.
        See _standardize function."""
        self._defocusU.set(defocusU)
        self._defocusV.set(defocusV)
        self._defocusAngle.set(defocusAngle)
        self.standardize()

    def standardize(self):
        """ Modify defocusU, defocusV and defocusAngle to conform
        the EMX standard: defocusU > defocusV, 0 <= defocusAngle < 180
        and the defocusAnges is between x-axis and defocusU. Also
        determine the defocusRatio(defocusU/defocusV).
        For more details see:
        http://i2pc.cnb.csic.es/emx/LoadDictionaryFormat.htm?type=Convention#ctf
        """
        if self._defocusV > self._defocusU:
            self._defocusV.swap(self._defocusU)  # exchange defocuU by defocuV
            self._defocusAngle.sum(90.)
        if self._defocusAngle >= 180.:
            self._defocusAngle.sum(-180.)
        elif self._defocusAngle < 0.:
            self._defocusAngle.sum(180.)
        # At this point defocusU is always greater than defocusV
        # following the EMX standard

        if self.getDefocusV() > self.DEFOCUS_V_MINIMUM_VALUE:
            self._defocusRatio.set(self.getDefocusU() / self.getDefocusV())
        else:
            self._defocusRatio.set(self.DEFOCUS_RATIO_ERROR_VALUE)

    def equalAttributes(self, other, ignore=[], verbose=False):
        """ Override default behaviour to compare two
        CTF objects, now ignoring the psdFile.
        """
        return (self._defocusU == other._defocusU and
                self._defocusV == other._defocusV and
                self._defocusAngle == other._defocusAngle
                )

    def setWrongDefocus(self):
        """ Set parameters if results parsing has failed.
        :param ctfModel: the model to be updated
        """
        self.setDefocusU(-999)
        self.setDefocusV(-1)
        self.setDefocusAngle(-999)


class DefocusGroup(EMObject):
    """ Groups CTFs by defocus"""

    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._defocusMin = Float()
        self._defocusMax = Float()
        self._defocusSum = Float(0)
        self._size = Integer(0)

    def getDefocusMin(self):
        return self._defocusMin.get()

    def getDefocusMax(self):
        return self._defocusMax.get()

    def getDefocusAvg(self):
        return self._defocusSum.get() / self.getSize()

    def getSize(self):
        return self._size.get()

    def addCTF(self, ctf):
        """ Add a new CTF to the group.
        Update values like min, max, avg and size.
        """
        self._size.increment()
        defocusU = ctf.getDefocusU()
        self._defocusMin.set(min(defocusU, self._defocusMin.get()))
        self._defocusMax.set(max(defocusU, self._defocusMax.get()))
        self._defocusSum.set(self._defocusSum.get() + defocusU)

    def containsCTF(self, ctf):
        """ Return True if a CTF is inside the group defocus range. """
        defocusU = ctf.getDefocusU()
        return self.getDefocusMax() >= defocusU >= self.getDefocusMin()


class SetOfDefocusGroups:
    """ Store a set of several defocus groups."""

    def __init__(self, inputSet,
                 groupRange=1000,
                 groupMinSize=1,
                 **kwargs):
        """ Create necessary defocus groups.
        Params:
            inputSet: input particles or micrographs with CTF information.
            groupRange: maximum defocus range allowed in one group.
            groupMinSize: impose a minimum number of particles per group.
        """
        self._groups = OrderedDict()
        self.__createGroups(inputSet, groupRange, groupMinSize)

    def __createGroups(self, inputSet, groupRange, groupMinSize):
        iterSet = iter(inputSet.iterItems(orderBy=['_ctfModel._defocusU',
                                                   'id'],
                                          direction='ASC'))
        first = next(iterSet)
        self.__addNewGroup(first.getCTF())

        for item in iterSet:
            ctf = item.getCTF()
            # Keep adding ctf to the current group
            # if the groupMinSize is not reached or
            # if the particle is inside the groupRange
            if (self._lastGroup.getSize() <
                    groupMinSize or
                    ctf.getDefocusU() - self._lastDefocus < groupRange):
                self._lastGroup.addCTF(ctf)
            else:
                self.__addNewGroup(ctf)

    def __addNewGroup(self, ctf):
        group = DefocusGroup()
        group.addCTF(ctf)
        count = len(self._groups) + 1
        defocusU = ctf.getDefocusU()
        groupName = 'ctfgroup_%06d_%05d' % (defocusU, count)
        self._groups[groupName] = group

        self._lastDefocus = defocusU
        self._lastGroup = group

    def getGroupByName(self, groupName):
        """ Return the Group for a given Name.
        If not exists, return None.
        """
        return self._groups.get(groupName, None)

    def getGroupByCTF(self, ctf):
        """ Get the CTF group for this ctf. """
        for group in self._groups.values():
            if group.containsCTF(ctf):
                return group
        return None

    def getDefocusMax(self):
        return self._defocusMax.get()

    def getDefocusAvg(self):
        return self._defocusAvg.get() / self.getSize()

    def getSize(self):
        return len(self._groups)


class ImageDim(CsvList):
    """ Just a wrapper to a CsvList to store image dimensions
    as X, Y and Z.
    """

    def __init__(self, x=None, y=None, z=None):
        CsvList.__init__(self, pType=int)
        if x is not None and y is not None:
            self.append(x)
            self.append(y)
            if z is not None:
                self.append(z)

    def getX(self):
        return self[0]

    def getY(self):
        return self[1]

    def getZ(self):
        return self[2]

    def __str__(self):
        if self.isEmpty():
            s = 'No-Dim'
        else:
            s = '%d x %d' % (self.getX(), self.getY())
            if self.getZ() > 1:
                s += ' x %d' % self.getZ()
        return s


class Image(EMObject):
    """Represents an EM Image object"""
    IMG_FILENAME_ATTR = '_filename'

    def __init__(self, location=None, **kwargs):
        """
         Params:
        :param location: Could be a valid location: (index, filename)
        or  filename
        """
        EMObject.__init__(self, **kwargs)
        # Image location is composed by an index and a filename
        self._index = Integer(0)
        self._filename = String()
        self._samplingRate = Float()
        self._ctfModel = None
        self._acquisition = None
        # _transform property will store the transformation matrix
        # this matrix can be used for 2D/3D alignment or
        # to represent projection directions
        self._transform = None
        # default origin by default is box center =
        # (Xdim/2, Ydim/2,Zdim/2)*sampling
        # origin stores a matrix that using as input the point (0,0,0)
        # provides  the position of the actual origin in the system of
        # coordinates of the default origin.
        # _origin is an object of the class Transform shifts
        # units are A.
        # Origin coordinates follow the MRC convention
        self._origin = None
        if location:
            self.setLocation(location)

    def getSamplingRate(self):
        """ Return image sampling rate. (A/pix) """
        return self._samplingRate.get()

    def setSamplingRate(self, sampling):
        self._samplingRate.set(sampling)

    def getFormat(self):
        pass

    def getDataType(self):
        pass

    def getDimensions(self):
        """getDimensions is redundant here but not in setOfVolumes
         create it makes easier to create protocols for both images
         and sets of images
        """
        return self.getDim()

    def getDim(self):
        """Return image dimensions as tuple: (Xdim, Ydim, Zdim)"""
        from pwem.emlib.image import ImageHandler
        x, y, z, n = ImageHandler().getDimensions(self)
        return None if x is None else (x, y, z)

    def getImage(self):
        """ Returns the actual image this objects represents"""
        from pwem.emlib.image import ImageHandler
        ih = ImageHandler()
        return ih.read(self)

    def getXDim(self):
        dim = self.getDim()
        return dim[0] if dim is not None else 0

    def getYDim(self):
        dim = self.getDim()
        return dim[1] if dim is not None else 0

    def getIndex(self):
        return self._index.get()

    def setIndex(self, index):
        self._index.set(index)

    def getFileName(self):
        """ Use the _objValue attribute to store filename. """
        return self._filename.get()

    def setFileName(self, filename):
        """ Use the _objValue attribute to store filename. """
        self._filename.set(filename)

    def getLocation(self):
        """ This function return the image index and filename.
        It will only differs from getFileName, when the image
        is contained in a stack and the index make sense.
        """
        return self.getIndex(), self.getFileName()

    def setLocation(self, *args):
        """ Set the image location, see getLocation.
        Params:
            First argument can be:
             1. a tuple with (index, filename)
             2. a index, this implies a second argument with filename
             3. a filename, this implies index=NO_INDEX
        """
        first = args[0]
        t = type(first)
        if t == tuple:
            index, filename = first
        elif t == int:
            index, filename = first, args[1]
        elif t == str:
            index, filename = NO_INDEX, first
        else:
            raise Exception('setLocation: unsupported type %s as input.' % t)

        self.setIndex(index)
        self.setFileName(filename)

    def getBaseName(self):
        return os.path.basename(self.getFileName())

    def copyInfo(self, other):
        """ Copy basic information """
        if type(self) is type(other):
            self.copy(other, copyId=False)
        else:
            self.copyAttributes(other, '_samplingRate')

    def copyLocation(self, other):
        """ Copy location index and filename from other image. """
        self.setIndex(other.getIndex())
        self.setFileName(other.getFileName())

    def hasCTF(self):
        return self._ctfModel is not None

    def getCTF(self):
        """ Return the CTF model """
        return self._ctfModel

    def setCTF(self, newCTF):
        self._ctfModel = newCTF

    def hasAcquisition(self):
        return (self._acquisition is not None and
                self._acquisition.getVoltage() is not None and
                self._acquisition.getMagnification() is not None
                )

    def getAcquisition(self):
        return self._acquisition

    def setAcquisition(self, acquisition):
        self._acquisition = acquisition

    def hasTransform(self):
        return self._transform is not None

    def getTransform(self) -> Transform:
        return self._transform

    def setTransform(self, newTransform):
        self._transform = newTransform

    def hasOrigin(self):
        return self._origin is not None

    def getOrigin(self, force=False):
        """shifts in A"""
        if self.hasOrigin():
            return self._origin
        else:
            if force:
                return self._getDefaultOrigin()
            else:
                return None

    def _getDefaultOrigin(self):
        sampling = self.getSamplingRate()
        t = Transform()
        x, y, z = self.getDim()
        if z > 1:
            z = z / -2.
        t.setShifts(x / -2. * sampling, y / -2. * sampling, z * sampling)
        return t  # The identity matrix

    def getShiftsFromOrigin(self):
        origin = self.getOrigin(force=True).getShifts()
        x = origin[0]
        y = origin[1]
        z = origin[2]
        return x, y, z
        # x, y, z are floats in Angstroms

    def setShiftsInOrigin(self, x, y, z):
        origin = self.getOrigin(force=True)
        origin.setShifts(x, y, z)

    def setOrigin(self, newOrigin=None):
        """If None, default origin will be set.
        Note: shifts are in Angstroms"""
        if newOrigin:
            self._origin = newOrigin
        else:
            self._origin = self._getDefaultOrigin()

    def originResampled(self, originNotResampled, oldSampling):
        factor = self.getSamplingRate() / oldSampling
        shifts = originNotResampled.getShifts()
        origin = self.getOrigin(force=True)
        origin.setShifts(shifts[0] * factor,
                         shifts[1] * factor,
                         shifts[2] * factor)
        return origin

    def __str__(self):
        """ String representation of an Image. """
        dim = self.getDim()
        dimStr = str(ImageDim(*dim)) if dim else 'No-Dim'
        return ("%s (%s, %0.2f Å/px)" % (self.getClassName(), dimStr,
                                         self.getSamplingRate() or 99999.))

    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        return filePaths

    def setMRCSamplingRate(self):
        """ Sets the sampling rate to the mrc file represented by this image"""
        from pwem.convert.headers import setMRCSamplingRate
        setMRCSamplingRate(self.getFileName(), self.getSamplingRate())


class Micrograph(Image):
    """ Represents an EM Micrograph object """
    MIC_NAME = '_micName'

    def __init__(self, location=None, **kwargs):
        Image.__init__(self, location, **kwargs)
        self._micName = String()

    def setMicName(self, micName):
        self._micName.set(micName)

    def getMicName(self):
        if self._micName.get():
            return self._micName.get()
        else:
            self.getFileName()

    def copyInfo(self, other):
        """ Copy basic information """
        Image.copyInfo(self, other)
        self.setMicName(other.getMicName())


class Particle(Image):
    """ Represents an EM Particle object """

    def __init__(self, location=None, **kwargs):
        Image.__init__(self, location, **kwargs)
        # This may be redundant, but make the Particle
        # object more independent for tracking coordinates
        self._coordinate = None
        self._micId = Integer()
        self._classId = Integer()

    def hasCoordinate(self):
        return self._coordinate is not None

    def setCoordinate(self, coordinate):
        self._coordinate = coordinate

    def getCoordinate(self):
        return self._coordinate

    def scaleCoordinate(self, factor):
        self.getCoordinate().scale(factor)

    def getMicId(self):
        """ Return the micrograph id if the coordinate is not None.
        or have set the _micId property.
        """
        if self._micId.hasValue():
            return self._micId.get()
        if self.hasCoordinate():
            return self.getCoordinate().getMicId()

        return None

    def setMicId(self, micId):
        self._micId.set(micId)

    def hasMicId(self):
        return self.getMicId() is not None

    def getClassId(self):
        return self._classId.get()

    def setClassId(self, classId):
        self._classId.set(classId)

    def hasClassId(self):
        return self._classId.hasValue()


class Mask(Particle):
    """ Represent a mask. """
    pass


class Volume(Image):
    """ Represents an EM Volume object """

    def __init__(self, location=None, **kwargs):
        Image.__init__(self, location, **kwargs)
        self._classId = Integer()
        self._halfMapFilenames = CsvList(pType=str)

    def getDim(self):
        """Return image dimensions as tuple: (Xdim, Ydim, Zdim)"""
        from pwem.emlib.image import ImageHandler

        fn = self.getFileName()
        if fn is not None and os.path.exists(fn.replace(':mrc', '')):
            x, y, z, n = ImageHandler().getDimensions(self)

            # Some volumes in mrc format can have the z dimension
            # as n dimension, so we need to consider this case.
            if z > 1:
                return x, y, z
            else:
                return x, y, n
        return None

    def getClassId(self):
        return self._classId.get()

    def setClassId(self, classId):
        self._classId.set(classId)

    def hasClassId(self):
        return self._classId.hasValue()

    def hasHalfMaps(self):
        return not self._halfMapFilenames.isEmpty()

    def getHalfMaps(self, asList=False):
        if asList:
            return self._halfMapFilenames
        else:
            return self._halfMapFilenames.get()

    def setHalfMaps(self, listFileNames):
        return self._halfMapFilenames.set(listFileNames)

    def fixMRCVolume(self, setSamplingRate=False):
        """ Fixes the header of the mrc file pointed by this object

        :param setSamplingRate: if true, it will set the header's sampling rate of the MRC file it refers

        """
        from pwem.convert.headers import fixVolume, setMRCSamplingRate
        fixVolume(self.getFileName())

        if setSamplingRate:
            setMRCSamplingRate(self.getFileName(), self.getSamplingRate())

    def __str__(self):
        """ returns string representation adding halves info to base image.__str__"""
        imgStr = super().__str__()

        if self.hasHalfMaps():
            imgStr += " - w/halves"

        return imgStr


class VolumeMask(Volume):
    """ A 3D mask to be used with volumes. """
    pass


class EMFile(EMObject):
    """ Class to link usually to text files. """

    def __init__(self, filename=None, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._filename = String(filename)

    def getFileName(self):
        """ Use the _objValue attribute to store filename. """
        return self._filename.get()

    def setFileName(self, filename):
        """ Use the _objValue attribute to store filename. """
        self._filename.set(filename)


class Alphabet():
    """ class with a dictionary of all valid alphabets"""
    # sequence types
    AMINOACIDS = 0
    NUCLEOTIDES = 1

    SEQ_TYPE = ['aminoacids', 'nucleotides']

    # alphabets for proteins
    PROTEIN_ALPHABET = 0
    EXTENDED_PROTEIN_ALPHABET = 1

    # alphabets for nucleotides
    AMBIGOUS_DNA_ALPHABET = 2
    UNAMBIGOUS_DNA_ALPHABET = 3
    EXTENDED_DNA_ALPHABET = 4
    AMBIGOUS_RNA_ALPHABET = 5
    UNAMBIGOUS_RNA_ALPHABET = 6
    NUCLEOTIDES_ALPHABET = 7

    # dummy alphabet
    DUMMY_ALPHABET = 8

    # dictionary with all alphabets
    alphabets = {};
    alphabetsLabels = {}

    alphabets[PROTEIN_ALPHABET] = 'ACDEFGHIKLMNPQRSTVWY'
    alphabets[EXTENDED_PROTEIN_ALPHABET] = 'ACDEFGHIKLMNPQRSTVWYBJOUXZ'
    alphabets[AMBIGOUS_DNA_ALPHABET] = 'GATCRYWSMKHBVDN'
    alphabets[UNAMBIGOUS_DNA_ALPHABET] = 'GATC'
    alphabets[EXTENDED_DNA_ALPHABET] = 'GATCBDSW'
    alphabets[AMBIGOUS_RNA_ALPHABET] = 'GAUCRYWSMKHBVDN'
    alphabets[UNAMBIGOUS_RNA_ALPHABET] = 'GAUC'
    alphabets[NUCLEOTIDES_ALPHABET] = 'GAC'
    alphabets[DUMMY_ALPHABET] = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    # dictionary with all alphabets labels
    alphabetsLabels[PROTEIN_ALPHABET] = 'Protein'
    alphabetsLabels[EXTENDED_PROTEIN_ALPHABET] = 'Extended Protein'
    alphabetsLabels[AMBIGOUS_DNA_ALPHABET] = 'Ambigous DNA'
    alphabetsLabels[UNAMBIGOUS_DNA_ALPHABET] = 'Unambigous DNA'
    alphabetsLabels[EXTENDED_DNA_ALPHABET] = 'Extended DNA'
    alphabetsLabels[AMBIGOUS_RNA_ALPHABET] = 'Ambigous RNA'
    alphabetsLabels[UNAMBIGOUS_RNA_ALPHABET] = 'Unambigous RNA'


class Sequence(EMObject):
    """Class containing a sequence of aminoacids/nucleotides
       Attribute names follow the biopython default ones
            param: name: name of the sequence
            param: sequence: string with the sequence
            param: alphabet: integer with the alphabet to be used
            param: isAminoacids: boolean indicating if the sequence is an aminoacid sequence
            param: id: string with the sequence id
            param: description: string with the sequence description
           
    """

    def __init__(self, name=None, sequence=None,
                 alphabet=None, isAminoacids=True, id=None, description=None,
                 **kwargs):
        EMObject.__init__(self, **kwargs)
        # sequence Id, usually from a database. E.g: P12345
        self._id = String(id)
        # Descriptive alias provided by the user. E.g: CicloxigenasaB
        self._name = String(name)
        # Id this a aminoacid or nucleotide sequence
        self._isAminoacids = Boolean(isAminoacids)
        # So far we just use _description when creating 'fasta' sequence files
        self._description = String(description)
        # sequence stores a string of residues (one character per residue)
        self._sequence = String(sequence)
        # alphabet is used to describe de convention followed to
        # store the _sequence. We follow biopython criteria
        self._alphabet = Integer(alphabet)

    def getId(self):
        return self._id.get()

    def setId(self, id):
        self._id.set(id)

    # Note that the natural name for the next two functions are
    # getName and  setName but
    # setName is defined in Object for another purpose
    def getSeqName(self):
        return self._name.get()

    def setSeqName(self, name):
        self._name.set(name)

    def getSequence(self):
        return self._sequence.get()

    def setSequence(self, sequence):
        self._sequence.set(sequence)

    def getDescription(self):
        return self._description.get()

    def setDescription(self, description):
        self._description.set(description)

    # Note: Alphabet is set when the sequence object is created
    # after that it makes no sense to change the alphabet
    def getAlphabet(self):
        return self._alphabet.get()

    def getIsAminoacids(self):
        return self._isAminoacids

    def importFromFile(self, seqFileName, isAmino=True):
        '''Fills the sequence attributes from the FIRST sequence found in the specified file'''
        import pwem.convert as emconv
        seqHandler = emconv.SequenceHandler()
        seqDic = seqHandler.readSequenceFromFile(seqFileName, type=None, isAmino=isAmino)
        self.setSequence(seqDic['sequence']), self.setId(seqDic['seqID'])
        self.setName(seqDic['name']), self.setDescription(seqDic['description'])
        self._alphabet.set(Integer(seqDic['alphabet']))
        self._isAminoacids.set(Boolean(seqDic['isAminoacids']))

    def exportToFile(self, seqFileName, doClean=True):
        '''Exports the sequence to the specified file'''
        import pwem.convert as emconv
        seqHandler = emconv.SequenceHandler(self.getSequence(),
                                            self._alphabet.get(), doClean)
        # retrieving  args from scipion object
        seqID = self.getId() if self.getId() is not None else 'seqID'
        seqName = self.getSeqName() if self.getSeqName() is not None else 'seqName'
        seqDescription = self.getDescription() if self.getDescription() is not None else ''
        seqHandler.saveFile(seqFileName, seqID,
                            name=seqName, seqDescription=seqDescription,
                            type=None)
        #seqFiP12345 USER_SEQ
        # Aspartate aminotransferase, mitochondrial

    def appendToFile(self, seqFileName, doClean=True):
        '''Exports the sequence to the specified file. If it already exists,
        the sequence is appended to the ones in the file'''
        logger.info("Appending sequence to file: %s" % seqFileName)
        import pwem.convert as emconv
        seqHandler = emconv.SequenceHandler(self.getSequence(),
                                            Alphabet.DUMMY_ALPHABET, doClean)
        # retrieving  args from scipion object
        seqID = self.getId() if self.getId() is not None else 'seqID'
        seqName = self.getSeqName() if self.getSeqName() is not None else 'seqName'
        seqDescription = self.getDescription() if self.getDescription() is not None else ''
        seqHandler.appendFile(seqFileName, seqID,
                              name=seqName, seqDescription=seqDescription,
                              type=None)

    def __str__(self):
        return "Sequence (name = {})\n".format(self.getSeqName())


class AtomStruct(EMFile):
    """Represents an PDB file. """

    def __init__(self, filename=None, pseudoatoms=False, **kwargs):
        EMFile.__init__(self, filename, **kwargs)
        self._pseudoatoms = Boolean(pseudoatoms)
        self._volume = None
        # origin stores a matrix that using as input the point (0,0,0)
        # provides  the position of the actual origin in the system of
        # coordinates of the default origin.
        # _origin is an object of the class Transform shifts
        # units are Angstroms (in Image units are A)
        self._origin = None

    def getPseudoAtoms(self):
        return self._pseudoatoms.get()

    def setPseudoAtoms(self, value):
        self._pseudoatoms.set(value)

    def getVolume(self):
        return self._volume

    def hasVolume(self):
        return self._volume is not None

    def setVolume(self, volume):
        if issubclass(type(volume), Volume):
            self._volume = volume
        else:
            raise Exception('TypeError', 'ERROR: SetVolume, This is not a '
                                         'volume')

    def __str__(self):
        return ("%s (pseudoatoms=%s, volume=%s)" %
                (self.getClassName(), self.getPseudoAtoms(), self.hasVolume()))

    def hasOrigin(self):
        return self._origin is not None

    def getOrigin(self, force=False):
        if self.hasOrigin():
            return self._origin
        else:
            if force:
                t = Transform()
                t.setShifts(0., 0., 0.)
                return t  # The identity matrix
            else:
                return None

    def setOrigin(self, newOrigin):
        self._origin = newOrigin


class PdbFile(AtomStruct):
    def __init__(self, filename=None, pseudoatoms=False, **kwargs):
        AtomStruct.__init__(self, filename, pseudoatoms, **kwargs)


class EMSet(Set, EMObject):
    _classesDict = None

    def _loadClassesDict(self):

        if self._classesDict is None:
            from pwem import Domain
            self._classesDict = Domain.getObjects()
            self._classesDict.update(globals())

        return self._classesDict

    def copyInfo(self, other):
        """ Define a dummy copyInfo function to be used
        for some generic operations on sets.
        """
        pass

    def clone(self):
        """ Override the clone defined in Object
        to avoid copying _mapperPath property
        """
        pass

    def copyItems(self, otherSet,
                  updateItemCallback=None,
                  itemDataIterator=None,
                  copyDisabled=False,
                  doClone=True,
                  itemSelectedCallback=None,
                  rowFilter=None,
                  orderBy='id',
                  direction='ASC'
                  ):
        """ Copy items from another set, allowing to update items information
        based on another source of data, paired with each item.

        Params:
            otherSet: input set from where the items will be copied.
            updateItemCallback: if not None, this will be called for each item
                and each data row (if the itemDataIterator is not None). Apart
                from setting item values or new attributes, it is possible to
                set the special flag _appendItem to False, and then this item
                will not be added to the set.
            itemDataIterator: if not None, it must be an iterator that have one
                data item for each item in the set. Usually the data item is a
                data row, coming from a table stored in text files (e.g STAR)
            copyDisabled: By default, disabled items are not copied from the other
                set. If copyDisable=True, then the enabled property of the item
                will be ignored.
            doClone: By default, the new item that will be inserted is a "clone"
                of the input item. By using doClone=False, the same input item
                will be passed to the callback and added to the set. This will
                avoid the clone operation and the related overhead.
            itemSelectedCallback: Optional, callback receiving an item and
                returning true if it has to be copied
            orderBy: Attribute by which the items will be sorted before copying. Default is 'id'.
            direction: Sorting direction, either 'ASC' for ascending or 'DESC' for descending. Default is 'ASC'.
        """

        if itemSelectedCallback is None:
            itemSelectedCallback = EMSet.isItemEnabled

        if isinstance(otherSet, Set):
            itemIterator = otherSet.iterItems(rowFilter=rowFilter,
                                              orderBy=orderBy,
                                              direction=direction)
        else:
            itemIterator = otherSet

        for item in itemIterator:
            # copy items if enabled or copyDisabled=True
            if copyDisabled or itemSelectedCallback(item):
                newItem = item.clone() if doClone else item
                if updateItemCallback:
                    row = None if itemDataIterator is None else next(itemDataIterator)
                    updateItemCallback(newItem, row)
                # If updateCallBack function returns attribute
                # _appendItem to False do not append the item
                if getattr(newItem, "_appendItem", True):
                    self.append(newItem)
            else:
                if itemDataIterator is not None:
                    next(itemDataIterator)  # just skip disabled data row

    @classmethod
    def create(cls, outputPath,
               prefix=None, suffix=None, ext=None,
               **kwargs):
        """ Create an empty set from the current Set class.
         Params:
            outputPath: where the output file will be written.
            prefix: prefix of the created file, if None, it will be deduced
                from the ClassName.
            suffix: additional suffix that will be added to the prefix with a
                "_" in between.
            ext: extension of the output file, be default will use .sqlite
        """
        fn = prefix or cls.__name__.lower().replace('setof', '')

        if suffix:
            fn += '_%s' % suffix

        if ext and ext[0] == '.':
            ext = ext[1:]
        fn += '.%s' % (ext or 'sqlite')

        setPath = os.path.join(outputPath, fn)
        pwutils.cleanPath(setPath)

        return cls(filename=setPath, **kwargs)

    def createCopy(self, outputPath,
                   prefix=None, suffix=None, ext=None,
                   copyInfo=False, copyItems=False,
                   itemSelectedCallback=None,
                   rowFilter=None):
        """ Make a copy of the current set to another location (e.g file).
        Params:
            outputPath: where the output file will be written.

            prefix: prefix of the created file, if None, it will be deduced
                from the ClassName.

            suffix: additional suffix that will be added to the prefix with a
                "_" in between.

            ext: extension of the output file, be default will use the same
                extension of this set filename.

            copyInfo: if True the copyInfo will be called after set creation.

            copyItems: if True the copyItems function will be called.

            itemSelectedCallback: Optional, callback receiving an item and returning
                true if it has to be copied
        """
        setObj = self.create(outputPath,
                             prefix=prefix,
                             suffix=suffix,
                             ext=ext or pwutils.getExt(self.getFileName()))

        if copyInfo:
            setObj.copyInfo(self)

        if copyItems:
            setObj.copyItems(self, itemSelectedCallback=itemSelectedCallback, rowFilter=rowFilter)

        return setObj

    def getFiles(self):
        return Set.getFiles(self)


class SetOfImages(EMSet):
    """ Represents a set of Images """
    ITEM_TYPE = Image
    _compatibilityDict = {'sampling rates': 'getSamplingRate',
                          'dimensions': 'getDimensions'}

    def __init__(self, **kwargs):
        EMSet.__init__(self, **kwargs)
        self._samplingRate = Float()
        self._hasCtf = Boolean(kwargs.get('ctf', False))
        self._alignment = String(ALIGN_NONE)
        self._isPhaseFlipped = Boolean(False)
        self._isAmplitudeCorrected = Boolean(False)
        self._acquisition = Acquisition()
        self._firstDim = ImageDim()  # Dimensions of the first image

    def getAcquisition(self):
        return self._acquisition

    def setAcquisition(self, acquisition):
        self._acquisition = acquisition

    def hasAcquisition(self):
        return self._acquisition.getMagnification() is not None

    def hasCTF(self):
        """Return True if the SetOfImages has associated a CTF model"""
        return self._hasCtf.get()

    def setHasCTF(self, value):
        self._hasCtf.set(value)

    def hasAlignment(self):
        return self._alignment != ALIGN_NONE

    def hasAlignment2D(self):
        return self._alignment == ALIGN_2D

    def hasAlignment3D(self):
        return self._alignment == ALIGN_3D

    def hasAlignmentProj(self):
        return self._alignment == ALIGN_PROJ

    def getAlignment(self):
        return self._alignment.get()

    def setAlignment(self, value):
        if value not in ALIGNMENTS:
            raise Exception('Invalid alignment value: "%s"' % value)
        self._alignment.set(value)

    def setAlignment2D(self):
        self.setAlignment(ALIGN_2D)

    def setAlignment3D(self):
        self.setAlignment(ALIGN_3D)

    def setAlignmentProj(self):
        self.setAlignment(ALIGN_PROJ)

    def isPhaseFlipped(self):
        return self._isPhaseFlipped.get()

    def setIsPhaseFlipped(self, value):
        self._isPhaseFlipped.set(value)

    def isAmplitudeCorrected(self):
        return self._isAmplitudeCorrected.get()

    def setIsAmplitudeCorrected(self, value):
        self._isAmplitudeCorrected.set(value)

    def append(self, image):
        """ Add a image to the set. """
        # If the sampling rate was set before, the same value
        # will be set for each image added to the set
        if self.getSamplingRate() or not image.getSamplingRate():
            image.setSamplingRate(self.getSamplingRate())
        # Copy the acquisition from the set to images
        # only override image acquisition if setofImages acquisition
        # is not none
        if self.hasAcquisition():
            # TODO: image acquisition should not be overwritten
            if not image.hasAcquisition():
                image.setAcquisition(self.getAcquisition())
        # Store the dimensions of the first image, just to
        # avoid reading image files for further queries to dimensions
        # only check this for first time append is called
        if self.isEmpty():
            self._setFirstDim(image)

        EMSet.append(self, image)

    def _setFirstDim(self, image):
        """ Store dimensions when the first image is found.
        This function should be called only once, to avoid reading
        dimension from image file. """
        logger.info("Getting the dimensions for the first item: %s" % image.getFileName())
        self._firstDim.set(image.getDim())

    def copyInfo(self, other):
        """ Copy basic information (sampling rate and ctf)
        from other set of images to current one"""
        if type(self) is type(other):
            self.copy(other, copyId=False)
        else:
            self.copyAttributes(other, '_samplingRate', '_isPhaseFlipped',
                                '_isAmplitudeCorrected', '_alignment')
            self._acquisition.copyInfo(other._acquisition)

    def getFiles(self):
        filePaths = set()
        uniqueFiles = self.aggregate(['count'], '_filename', ['_filename'])

        for row in uniqueFiles:
            filePaths.add(row['_filename'])
        return filePaths

    def setDownsample(self, downFactor):
        """ Update the values of samplingRate and scannedPixelSize
        after applying a downsampling factor of downFactor.
        """
        self.setSamplingRate(self.getSamplingRate() * downFactor)

    def setSamplingRate(self, samplingRate):
        """ Set the sampling rate and adjust the scannedPixelSize. """
        self._samplingRate.set(samplingRate)

    def getSamplingRate(self):
        return self._samplingRate.get()

    def writeStack(self, fnStack, orderBy='id', direction='ASC',
                   applyTransform=False):
        # TODO create empty file to improve efficiency
        from pwem.emlib.image import ImageHandler
        ih = ImageHandler()
        applyTransform = applyTransform and self.hasAlignment2D()

        for i, img in enumerate(self.iterItems(orderBy=orderBy,
                                               direction=direction)):
            transform = img.getTransform() if applyTransform else None
            ih.convert(img, (i + 1, fnStack), transform=transform)

    # TODO: Check whether this function can be used.
    # for example: protocol_apply_mask
    def readStack(self, fnStack, postprocessImage=None):
        """ Populate the set with the images in the stack """
        from pwem.emlib.image import ImageHandler

        _, _, _, ndim = ImageHandler().getDimensions(fnStack)
        img = self.ITEM_TYPE()
        for i in range(1, ndim + 1):
            img.setObjId(None)
            img.setLocation(i, fnStack)
            if postprocessImage is not None:
                postprocessImage(img)
            self.append(img)

    def getDim(self):
        """ Return the dimensions of the first image in the set. """
        if self._firstDim.isEmpty():
            return None
        x, y, z = self._firstDim
        return x, y, z

    def setDim(self, newDim):
        self._firstDim.set(newDim)

    def getXDim(self):
        return self.getDim()[0] if self.getDim() is not None else 0

    def isOddX(self):
        """ Return True if the first item x dimension is odd. """
        return self.getXDim() % 2 == 1

    def getDimensions(self):
        """Return first image dimensions as a tuple: (xdim, ydim, zdim)"""
        return self.getFirstItem().getDim()

    def getClassNameStr(self):
        return self.getClassName().replace("SetOf", '')

    def __str__(self):
        """ String representation of a set of images. """
        try:
            s = "%s (%d items, %s, %s%s)" % \
                (self.getClassNameStr(), self.getSize(),
                 self._dimStr(), self._samplingRateStr(), self._appendStreamState())
        except Exception as e:
            s = "Couldn't convert the set to text."
            logger.error(s, exc_info=e)

        return s

    def _samplingRateStr(self):
        """ Returns how the sampling rate is presented in a 'str' context."""
        sampling = self.getSamplingRate()

        if not sampling:
            logger.error("FATAL ERROR: Object %s has no sampling rate!!!"
                         % self.getName())
            sampling = -999.0

        return "%0.2f Å/px" % sampling

    def _dimStr(self):
        """ Return the string representing the dimensions. """
        return str(self._firstDim)

    def iterItems(self, orderBy='id', direction='ASC', where='1', limit=None, iterate=True, rowFilter=None):
        """ Redefine iteration to set the acquisition to images. """
        for img in Set.iterItems(self, orderBy=orderBy, direction=direction,
                                 where=where, limit=limit, iterate=iterate, rowFilter=rowFilter):

            # Sometimes the images items in the set could
            # have the acquisition info per data row, and we
            # don't want to override with the set acquisition for this case
            if not img.hasAcquisition():
                img.setAcquisition(self.getAcquisition())
            yield img

    def appendFromImages(self, imagesSet, itemSelectedCallback=None, rowFilter=None):
        """ Iterate over the images and append
        every image that is enabled.

        :param imagesSet: Set to go copy items from
        :param itemSelectedCallback: Optional, callback receiving an item and returning true if it has to be added

        """

        if itemSelectedCallback is None:
            itemSelectedCallback = SetOfImages.isItemEnabled

        for img in imagesSet.iterItems(rowFilter=rowFilter):
            if itemSelectedCallback(img):
                self.append(img)

    def appendFromClasses(self, classesSet, filterClassFunc=None):
        """ Iterate over the classes and the element inside each
        class and append to the set all that are enabled.
        """

        if filterClassFunc is None:
            filterClassFunc = SetOfImages.isItemEnabled

        for cls in classesSet.iterItems():
            if filterClassFunc(cls) and cls.getSize() > 0:
                for img in cls:
                    if img.isEnabled():
                        self.append(img)


class SetOfMicrographsBase(SetOfImages):
    """ Create a base class for both Micrographs and Movies,
    but avoid to select Movies when Micrographs are required.
    """

    def __init__(self, **kwargs):
        SetOfImages.__init__(self, **kwargs)
        self._scannedPixelSize = Float()

    def copyInfo(self, other):
        """ Copy basic information (voltage, spherical aberration and
        sampling rate) from other set of micrographs to current one.
        """
        SetOfImages.copyInfo(self, other)
        self._scannedPixelSize.set(other.getScannedPixelSize())

    def setSamplingRate(self, samplingRate):
        """ Set the sampling rate and adjust the scannedPixelSize. """
        self._samplingRate.set(samplingRate)
        mag = self._acquisition.getMagnification()
        if mag is None:
            self._scannedPixelSize.set(None)
        else:
            self._scannedPixelSize.set(1e-4 * samplingRate * mag)

    def getScannedPixelSize(self):
        return self._scannedPixelSize.get()

    def setScannedPixelSize(self, scannedPixelSize):
        """ Set scannedPixelSize and update samplingRate. """
        mag = self._acquisition.getMagnification()
        if mag is None:
            raise Exception("SetOfMicrographs: cannot set scanned pixel size "
                            "if Magnification is not set.")
        self._scannedPixelSize.set(scannedPixelSize)
        self._samplingRate.set((1e+4 * scannedPixelSize) / mag)


class SetOfMicrographs(SetOfMicrographsBase):
    """Represents a set of Micrographs"""
    ITEM_TYPE = Micrograph


class SetOfParticles(SetOfImages):
    """ Represents a set of Particles.
    The purpose of this class is to separate the
    concepts of Micrographs and Particles, even if
    both are considered Images
    """
    ITEM_TYPE = Particle
    REP_TYPE = Particle

    def __init__(self, **kwargs):
        SetOfImages.__init__(self, **kwargs)
        self._coordsPointer = Pointer()
        self._isSubParticles = Boolean(False)

    def getClassNameStr(self):
        if self._isSubParticles.get():
            return 'SubParticles'
        else:
            return super().getClassNameStr()

    def hasCoordinates(self):
        return self._coordsPointer.hasValue()

    def getCoordinates(self):
        """ Returns the SetOfCoordinates associated with
        this SetOfParticles"""
        return self._coordsPointer.get()

    def getIsSubparticles(self):
        return self._isSubParticles

    def setIsSubparticles(self, value):
        self._isSubParticles = Boolean(value)

    def setCoordinates(self, coordinates):
        """ Set the SetOfCoordinates associated with
        this set of particles.
         """
        if coordinates.isPointer():
            self._coordsPointer.copy(coordinates)
        else:
            self._coordsPointer.set(coordinates)

        if not self._coordsPointer.hasExtended():
            logger.warning("FOR DEVELOPERS: Direct pointers to objects should be avoided. "
                           "They are problematic in complex streaming scenarios. "
                           "Pass a pointer to a protocol with extended "
                           "(e.g.: input param are this kind of pointers. Without get()!)")

    def copyInfo(self, other):
        """ Copy basic information (voltage, spherical aberration and
        sampling rate) from other set of micrographs to current one.
        """
        SetOfImages.copyInfo(self, other)
        if hasattr(other, '_coordsPointer'):
            self.copyAttributes(other, "_coordsPointer")
        self.setHasCTF(other.hasCTF())


class SetOfAverages(SetOfParticles):
    """Represents a set of Averages.
    It is a SetOfParticles but it is useful to differentiate outputs."""

    def __init__(self, **kwargs):
        SetOfParticles.__init__(self, **kwargs)


class SetOfVolumes(SetOfImages):
    """Represents a set of Volumes"""
    ITEM_TYPE = Volume
    REP_TYPE = Volume

    # Hint to GUI components to expose internal items for direct selection
    EXPOSE_ITEMS = True

    def __init__(self, **kwargs):
        SetOfImages.__init__(self, **kwargs)


class SetOfMorphing(SetOfVolumes):
    def __init__(self, **kwargs):
        SetOfVolumes.__init__(self, **kwargs)


class SetOfCTF(EMSet):
    """ Contains a set of CTF models estimated for a set of images."""
    ITEM_TYPE = CTFModel

    def __init__(self, **kwargs):
        EMSet.__init__(self, **kwargs)
        self._micrographsPointer = Pointer()

    def getMicrographs(self):
        """ Return the SetOFImages used to create the SetOfClasses. """
        return self._micrographsPointer.get()

    def setMicrographs(self, micrographs):
        """ Set the micrographs from which this CTFs were estimated.
        Params:
            micrographs: Either a SetOfMicrographs object or a pointer to it.
        """
        if micrographs.isPointer():
            self._micrographsPointer.copy(micrographs)
        else:
            self._micrographsPointer.set(micrographs)


class SetOfDefocusGroup(EMSet):
    """ Contains a set of DefocusGroup.
        if min/max/avg exists the corresponding flag must be
        set to true.
    """
    ITEM_TYPE = DefocusGroup

    def __init__(self, **kwargs):
        EMSet.__init__(self, **kwargs)
        self._minSet = False
        self._maxSet = False
        self._avgSet = False

    def getMinSet(self):
        return self._minSet

    def setMinSet(self, value):
        self._minSet = value

    def getMaxSet(self):
        return self._maxSet

    def setMaxSet(self, value):
        self._maxSet = value

    def getAvgSet(self):
        return self._avgSet

    def setAvgSet(self, value):
        self._avgSet = value


class SetOfAtomStructs(EMSet):
    """ Set containing PDB items. """
    ITEM_TYPE = AtomStruct
    # Hint to GUI components to expose internal items for direct selection
    EXPOSE_ITEMS = True

    def getFiles(self):
        files = []
        for atomStruct in self.iterItems():
            files.append(atomStruct.getFileName())

        return files


class SetOfPDBs(SetOfAtomStructs):
    """ Set containing PDB items. """

    def __init__(self):
        SetOfAtomStructs.__init__(self)
        logger.warning("FOR DEVELOPERS: SetOfPDBs class has been renamed to SetOfAtomStructs. "
                       "Please update your code.")


class SetOfSequences(EMSet):
    """Set containing Sequence items."""
    ITEM_TYPE = Sequence

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.aligned = Boolean(kwargs.get('aligned', False))

    def exportToFile(self, seqFileName):
        '''Writes the sequences in the set in the specified file'''
        for sequence in self:
            sequence.appendToFile(seqFileName)

    def importFromFile(self, seqFileName, isAmino=True, type=None):
        '''Appends elements to the set from sequences found in the specified file'''
        import pwem.convert as emconv
        seqHandler = emconv.SequenceHandler()
        seqsDic = seqHandler.readSequencesFromFile(seqFileName, type=type, isAmino=isAmino)
        for seqDic in seqsDic:
            newSeq = Sequence(sequence=seqDic['sequence'], id=seqDic['seqID'],
                              name=seqDic['name'], description=seqDic['description'],
                              alphabet=seqDic['alphabet'], isAminoacids=seqDic['isAminoacids'])
            self.append(newSeq)


class Coordinate(EMObject):
    """This class holds the (x,y) position and other information
    associated with a coordinate"""

    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._micrographPointer = Pointer(objDoStore=False)
        self._x = Integer(kwargs.get('x', None))
        self._y = Integer(kwargs.get('y', None))
        self._micId = Integer()
        self._micName = String()

    def getX(self):
        return self._x.get()

    def setX(self, x):
        self._x.set(x)

    def shiftX(self, shiftX):
        self._x.sum(shiftX)

    def getY(self):
        return self._y.get()

    def setY(self, y):
        self._y.set(y)

    def shiftY(self, shiftY):
        self._y.sum(shiftY)

    def scale(self, factor):
        """ Scale both x and y coordinates by a given factor.
        """
        self._x.multiply(factor)
        self._y.multiply(factor)

    def getPosition(self):
        """ Return the position of the coordinate as a (x, y) tuple.
        mode: select if the position is the center of the box
        or in the top left corner.
        """
        return self.getX(), self.getY()

    def setPosition(self, x, y):
        self.setX(x)
        self.setY(y)

    def getMicrograph(self):
        """ Return the micrograph object to which
        this coordinate is associated.
        """
        return self._micrographPointer.get()

    def setMicrograph(self, micrograph):
        """ Set the micrograph to which this coordinate belongs. """
        self._micrographPointer.set(micrograph)
        self._micId.set(micrograph.getObjId())
        if isinstance(micrograph, Micrograph):
            self._micName.set(micrograph.getMicName())

    def copyInfo(self, coord):
        """ Copy information from other coordinate. """
        self.setPosition(*coord.getPosition())
        self.setObjId(coord.getObjId())

    def getMicId(self):
        return self._micId.get()

    def setMicId(self, micId):
        self._micId.set(micId)

    def invertY(self):
        if not self.getMicrograph() is None:
            dims = self.getMicrograph().getDim()
            height = dims[1]
            self.setY(height - self.getY())
        # else: error TODO

    def setMicName(self, micName):
        self._micName.set(micName)

    def getMicName(self):
        return self._micName.get()


class SetOfCoordinates(EMSet):
    """ Encapsulate the logic of a set of particles coordinates.
    Each coordinate has a (x,y) position and is related to a Micrograph
    The SetOfCoordinates can also have information about TiltPairs.
    """
    ITEM_TYPE = Coordinate

    def __init__(self, **kwargs):
        EMSet.__init__(self, **kwargs)
        self._micrographsPointer = Pointer()
        self._boxSize = Integer()
        self._micrographs = None

    def getBoxSize(self):
        """ Return the box size of the particles.
        """
        return self._boxSize.get()

    def setBoxSize(self, boxSize):
        """ Set the box size of the particles. """
        self._boxSize.set(boxSize)

    def iterMicrographs(self):
        """ Iterate over the micrographs set associated with this
        set of coordinates.
        """
        return self.getMicrographs()

    def iterMicrographCoordinates(self, micrograph):
        """ Iterates over the set of coordinates belonging to that micrograph.
        """
        pass

    def iterCoordinates(self, micrograph=None):
        """ Iterate over the coordinates associated with a micrograph.
        If micrograph=None, the iteration is performed over the whole
        set of coordinates.
        """
        if micrograph is None:
            micId = None
        elif isinstance(micrograph, int):
            micId = micrograph
        elif isinstance(micrograph, Micrograph):
            micId = micrograph.getObjId()
        else:
            raise Exception('Invalid input micrograph of type %s'
                            % type(micrograph))

        # Iterate over all coordinates if micId is None,
        # otherwise use micId to filter the where selection
        coordWhere = '1' if micId is None else '_micId=%d' % micId

        for coord in self.iterItems(where=coordWhere):
            # Associate the micrograph
            self._associateMicrograph(coord)
            yield coord

    def _associateMicrograph(self, coord):
        coord.setMicrograph(self._getMicrograph(coord.getMicId()))

    def _getMicrograph(self, micId):
        """ Returns  the micrograph from a micId"""
        micrographs = self._getMicrographs()

        if micId not in micrographs.keys():
            mic = self.getMicrographs()[micId]
            self._micrographs[micId] = mic
            return mic
        else:
            return micrographs[micId]

    def initMicrographs(self):
        """ Initialize internal _micrographs to a dictionary if not done already"""
        if self._micrographs is None:
            self._micrographs = dict()

    def _getMicrographs(self):
        if self._micrographs is None:
            self.initMicrographs()

        return self._micrographs

    def getMicrographs(self, asPointer=False):
        """ Returns the SetOfMicrographs associated with
        this SetOfCoordinates"""

        if asPointer:
            return self._micrographsPointer
        else:
            return self._micrographsPointer.get()

    def setMicrographs(self, micrographs):
        """ Set the micrographs associated with this set of coordinates.
        Params:
            micrographs: Either a SetOfMicrographs object or a pointer to it.
        """
        if micrographs.isPointer():
            self._micrographsPointer.copy(micrographs)
        else:
            self._micrographsPointer.set(micrographs)

        if not self._micrographsPointer.hasExtended():
            logger.warning("FOR DEVELOPERS: Direct pointers to objects should be avoided. "
                           "They are problematic in complex streaming scenarios. "
                           "Pass a pointer to a protocol with extended "
                           "(e.g.: input param are this kind of pointers. Without get()!)")

    def getFiles(self):
        filePaths = set()
        filePaths.add(self.getFileName())
        return filePaths

    def __str__(self):
        """ String representation of a set of coordinates. """
        if self._boxSize.hasValue():
            boxSize = self._boxSize.get()
            boxStr = ' %d x %d' % (boxSize, boxSize)
        else:
            boxStr = 'No-Box'
        s = "%s (%d items, %s%s)" % (self.getClassName(), self.getSize(),
                                     boxStr, self._appendStreamState())

        return s

    def copyInfo(self, other):
        """ Copy basic information (boxsize)
                from other set of coordinates to current one"""
        self.copyAttributes(other, '_boxSize')

        # TODO: we might what here to copy the mics too, same as done with
        # acquisition in SetOfImages
        if isinstance(other, SetOfCoordinates):
            self.setMicrographs(other.getMicrographs(asPointer=True))


class Matrix(Scalar):
    def __init__(self, **kwargs):
        Scalar.__init__(self, **kwargs)
        self._matrix = np.eye(4)

    def _convertValue(self, value):
        """Value should be a str with comma separated values
        or a list.
        """
        self._matrix = np.array(json.loads(value))

    def getObjValue(self):
        self._objValue = json.dumps(self._matrix.tolist())
        return self._objValue

    def setValue(self, i, j, value):
        self._matrix[i, j] = value

    def getMatrix(self):
        """ Return internal numpy matrix. """
        return self._matrix

    def setMatrix(self, matrix):
        """ Override internal numpy matrix. """
        self._matrix = matrix

    def __str__(self):
        return np.array_str(self._matrix)

    def _copy(self, other, copyDict, copyId, level=1, ignoreAttrs=[], copyEnable=False):
        """ Override the default behaviour of copy
        to also copy array data.
        """
        self.setMatrix(np.copy(other.getMatrix()))
        self._objValue = other._objValue


class Class2D(SetOfParticles):
    """ Represent a Class that groups Particles objects.
    Usually the representative of the class is another Particle
    (some kind of average particle from the particles assigned
    to the class)
    """

    def copyInfo(self, other):
        """ Copy basic information (id and other properties) but
        not _mapperPath or _size from other set of micrographs to current one.
        """
        self.copy(other, copyId=False, ignoreAttrs=['_mapperPath', '_size'])

    def clone(self):
        clone = self.getClass()()
        clone.copy(self, ignoreAttrs=['_mapperPath', '_size'])
        return clone

    def close(self):
        # Do nothing on close, since the db will be closed by SetOfClasses
        pass


class Class3D(SetOfParticles):
    """ Represent a Class that groups Particles objects.
    Usually the representative of the class is a Volume
    reconstructed from the particles assigned to the class.
    """
    REP_TYPE = Volume

    def copyInfo(self, other):
        """ Copy basic information (id and other properties) but not
        _mapperPath or _size from other set of micrographs to current one.
        """
        self.copy(other, copyId=False, ignoreAttrs=['_mapperPath', '_size'])

    def clone(self):
        clone = self.getClass()()
        clone.copy(self, ignoreAttrs=['_mapperPath', '_size'])
        return clone

    def close(self):
        # Do nothing on close, since the db will be closed by SetOfClasses
        pass


class ClassVol(SetOfVolumes):
    """ Represent a Class that groups Volume objects.
    Usually the representative of the class is another Volume.
    """

    def close(self):
        # Do nothing on close, since the db will be closed by SetOfClasses
        pass


class SetOfClasses(EMSet):
    """ Store results from a classification. """
    ITEM_TYPE = None  # type of classes stored in the set
    REP_TYPE = None  # type of the representatives of each class

    def __init__(self, **kwargs):
        EMSet.__init__(self, **kwargs)
        # Store the average images of each class(SetOfParticles)
        self._representatives = Boolean(False)
        self._imagesPointer = Pointer()

    def copyInfo(self, other):
        """ Copy basic information from other set of classes to current one"""
        self.copyAttributes(other, '_representatives', '_imagesPointer')

    def iterClassItems(self, iterDisabled=False):
        """ Iterate over the images of a class.
        Params:
            iterDisabled: If True, also include the disabled items. """
        for cls in self.iterItems():
            if iterDisabled or cls.isEnabled():
                for img in cls:
                    if iterDisabled or img.isEnabled():
                        yield img

    def hasRepresentatives(self):
        return self._representatives.get()

    def getImages(self):
        """ Return the SetOfImages used to create the SetOfClasses. """
        return self._imagesPointer.get()

    def getImagesPointer(self):
        """" Return the pointer to the SetOfImages used to create the SetOfClasses. """
        return self._imagesPointer

    def setImages(self, images):
        """ Set the images (particles) 2d associated with this set of classes.
        Params:
            images: An indirect pointer (with extended) to a set of images.
        """
        if images.isPointer():
            self._imagesPointer.copy(images)
        else:
            self._imagesPointer.set(images)

        if not self._imagesPointer.hasExtended():
            logger.warning("FOR DEVELOPERS: Direct pointers to objects should be avoided. "
                           "They are problematic in complex streaming scenarios. "
                           "Pass a pointer to a protocol with extended "
                           "(e.g.: input param are this kind of pointers. Without get()!)")

    def getDimensions(self):
        """Return first image dimensions as a tuple: (xdim, ydim, zdim)"""
        if self.hasRepresentatives():
            return self.getFirstItem().getRepresentative().getDim()
        return None

    def _setItemMapperPath(self, classItem):
        """ Set the mapper path of this class according to the mapper
        path of the SetOfClasses and also the prefix according to class id
        """
        classPrefix = 'Class%03d' % classItem.getObjId()
        classItem._mapperPath.set('%s,%s' % (self.getFileName(), classPrefix))
        classItem._mapperPath.setStore(False)
        classItem.load()

    def _insertItem(self, classItem):
        """ Create the SetOfImages assigned to a class.
        If the file exists, it will load the Set.
        """
        if classItem.hasRepresentative():
            self._representatives.set(True)

        self._setItemMapperPath(classItem)
        EMSet._insertItem(self, classItem)
        classItem.write(properties=False)  # Set.write(self)

    def __getitem__(self, itemId):
        """ Setup the mapper classes before returning the item. """
        classItem = EMSet.__getitem__(self, itemId)
        self._setItemMapperPath(classItem)
        return classItem

    def getFirstItem(self):
        classItem = EMSet.getFirstItem(self)
        self._setItemMapperPath(classItem)
        return classItem

    def iterItems(self, orderBy='id', direction='ASC', rowFilter=None):
        for classItem in EMSet.iterItems(self, orderBy=orderBy,
                                         direction=direction,
                                         rowFilter=rowFilter):
            self._setItemMapperPath(classItem)
            yield classItem

    def iterRepresentatives(self, orderBy='id', direction='ASC'):
        for classItem in self.iterItems(orderBy, direction):
            if classItem.hasRepresentative():
                rep = classItem.getRepresentative()
                rep.setObjId(classItem.getObjId())
                rep.setSamplingRate(classItem.getSamplingRate())

                yield rep

    def getFiles(self):

        files = []
        for rep in self.iterRepresentatives():
            files.append(rep.getFileName())
        return files

    def getSamplingRate(self):
        return self.getImages().getSamplingRate()

    def appendFromClasses(self, classesSet, filterClassFunc=None, updateClassCallback=None):
        """ Iterate over the classes and the elements inside each
        class and append classes and items that are enabled.

        :param classesSet: Set of classes to copy items from
        :param filterClassFunc: Extra callback to exclude classes. Receives a class item, should return a boolean
        """
        if filterClassFunc is None:
            filterClassFunc = lambda cls: True

        for cls in classesSet.iterItems():
            if cls.isEnabled() and filterClassFunc(cls):
                newCls = self.ITEM_TYPE()
                newCls.copyInfo(cls)
                newCls.setObjId(cls.getObjId())
                if updateClassCallback:
                    updateClassCallback(newCls)
                self.append(newCls)
                for img in cls:
                    if img.isEnabled():
                        newCls.append(img)
                self.update(newCls)

    def classifyItems(self,
                      updateItemCallback=None,
                      updateClassCallback=None,
                      itemDataIterator=None,
                      classifyDisabled=False,
                      iterParams=None,
                      doClone=True,
                      raiseOnNextFailure=True,
                      cancelNextWhenAppendIsFalse=False):
        """ Classify items from the self.getImages() and add the needed classes.
        This function iterates over each item in the images and call
        the updateItemCallback to register the information coming from
        the iterator in itemDataIterator. The callback function should
        set the classId of the image that will be used to classify it.
        It is also possible to pass a callback to update the class properties.

        :param updateItemCallback: callback to be invoked on each item's loop (e.g.: 2d image in a 2d classification)
        :param updateClassCallback: callback to be invoked when a item.getClassId() changes
        :param itemDataIterator: an iterator (usually on metadata files, star, xmd,...) that will be called on each loop.
        usually has that same lines as items and iteration is kept in sync
        :param classifyDisabled: classify disabled items. By default, they are skipped.
        :param iterParams: Parameters for self.getImages() to leave oot images/filter
        :param doClone: Make a clone of the item (defaults to true)
        :param raiseOnNextFailure: A boolean flag indicating whether to raise an exception if there is a failure while
                                   attempting to retrieve the next element from itemDataIterator. If set to True(default),
                                   an exception will be raised; if set to False, the loop will be terminated in case
                                   of failure. Pass False when itemDataIterator is a subset of the inputSet.
                                   Important: Pass the right iterParams to make sure the iteration on the inputSet
                                              matches the iteration of the itemDataIterator
        :param cancelNextWhenAppendIsFalse: A boolean flag that determines whether the classification process should skip
                                    the next item when the `_appendItem` attribute of the current item is set to `False`.
                                    If set to `True`, after encountering an item where `_appendItem` is `False`, the
                                    next item in the iteration will be skipped. If set to `False` (default),
                                    classification continues normally without skipping the next item.

        """

        # Dictionary to store the {classId: class} pairs
        clsDict = self._getExistingItems()

        inputSet = self.getImages()
        iterParams = iterParams or {}
        cancelNext = False

        # For each item in the input set: Particles tipically (which will contribute to the main items here: class2d or 3d).
        for item in inputSet.iterItems(**iterParams):
            # copy items if enabled or copyDisabled=True
            if classifyDisabled or item.isEnabled():

                # get the new item (cloned or not)
                newItem = item.clone() if doClone else item

                # If we need have a callback to update the item
                if updateItemCallback:

                    action, cancelNext = self._updateItem(cancelNext, cancelNextWhenAppendIsFalse, itemDataIterator, newItem,
                                     raiseOnNextFailure, updateItemCallback)

                    if action == LoopActions.BREAK:
                        break
                    elif action == LoopActions.CONTINUE:
                        continue

                # Get the class id (reference)
                ref = newItem.getClassId()
                if ref is None:
                    raise Exception('Particle classId is None!!!')
                if ref == 0:
                    continue

                # Get the class the newItem belongs to.
                classItem = self._get_or_create_class(clsDict, ref, updateClassCallback)
                classItem.append(newItem)
                # cancel next() cancelation --> Enable next()
                cancelNext = False
            else:
                if itemDataIterator is not None:
                    next(itemDataIterator)  # just skip disabled data row

        for classItem in clsDict.values():
            self.update(classItem)

    def _updateItem(self, cancelNext, cancelNextWhenAppendIsFalse, itemDataIterator, newItem, raiseOnNextFailure,
                    updateItemCallback)->LoopActions:
        # Declare row
        row = None

        # If next is not canceled and have and iterator to do the next()
        if not cancelNext and itemDataIterator is not None:
            try:
                # ... do the next
                row = next(itemDataIterator)
            except Exception as ex:
                # ... if need to raise an exception
                if raiseOnNextFailure:
                    raise ex
                else:
                    # tolerate next execeptions
                    return LoopActions.BREAK, cancelNext
        # Update items using callback
        try:
            updateItemCallback(newItem, row)
        except Exception as ex:
            logger.error("There was an error updating the particle %s (row: %s): %s" % (
                newItem.getObjId(), str(row), str(ex)), exc_info=ex)
            raise ex
        # If updateCallBack function returns attribute
        # _appendItem to False do not append the item
        if not getattr(newItem, "_appendItem", True):
            cancelNext = cancelNextWhenAppendIsFalse
            return LoopActions.CONTINUE, cancelNext

        return LoopActions.NONE, cancelNext

    def _get_or_create_class(self, clsDict, ref, updateClassCallback):
        # Register a new class set if the ref was not found.
        # if not ref in clsDict:

        inputSet = self.getImages()

        if ref not in clsDict:
            classItem = self.ITEM_TYPE(objId=ref)
            rep = self.REP_TYPE()
            classItem.setRepresentative(rep)
            clsDict[ref] = classItem
            classItem.copyInfo(inputSet)
            classItem.setAcquisition(inputSet.getAcquisition())
            if updateClassCallback is not None:
                updateClassCallback(classItem)
            self.append(classItem)
        else:
            classItem = clsDict[ref]
        return classItem

    def _getExistingItems(self):
        """ Returns a dictionary with the class id as key and the class as value asa representation
        of what already is in the set"""
        clsDict = dict()

        if not self.isEmpty():
            for item in self.iterItems():
                # clone with a param to clone also mapper path?
                clone = item.clone()
                self._setItemMapperPath(clone)
                # Maybe enableAppend of class based on enableAppend of set?
                clone.enableAppend()
                clsDict[item.getObjId()] = clone

        return clsDict

class SetOfClasses2D(SetOfClasses):
    """ Store results from a 2D classification of Particles. """
    ITEM_TYPE = Class2D
    REP_TYPE = Particle

    def writeStack(self, fnStack):
        """ Write an stack with the classes averages. """
        from pwem.emlib.image import ImageHandler

        if not self.hasRepresentatives():
            raise Exception('Could not write Averages stack '
                            'if not hasRepresentatives!!!')
        ih = ImageHandler()

        for i, class2D in enumerate(self):
            img = class2D.getRepresentative()
            ih.convert(img, (i + 1, fnStack))


class SetOfClasses3D(SetOfClasses):
    """ Store results from a 3D classification of Particles. """
    ITEM_TYPE = Class3D
    REP_TYPE = Volume

    pass


class SetOfClassesVol(SetOfClasses3D):
    """ Store results from a classification of Volumes. """
    ITEM_TYPE = ClassVol

    pass


class NormalMode(EMObject):
    """ Store normal mode information. """

    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._modeFile = String(kwargs.get('modeFile', None))
        self._collectivity = Float(kwargs.get('collectivity', None))
        self._score = Float(kwargs.get('score', None))

    def getModeFile(self):
        return self._modeFile.get()

    def setModeFile(self, value):
        self._modeFile.set(value)

    def getCollectivity(self):
        return self._collectivity.get()

    def setCollectivity(self, value):
        self._collectivity.set(value)

    def getScore(self):
        return self._score.get()

    def setScore(self, value):
        self._score.set(value)


class SetOfNormalModes(EMSet):
    """ Set containing NormalMode items. """
    ITEM_TYPE = NormalMode

    def __init__(self, **kwargs):
        EMSet.__init__(self, **kwargs)
        # Store a pointer to the PdbFile object
        # from which this normal modes where computed.
        self._pdbPointer = Pointer()

    def setPdb(self, pdbFile):
        self._pdbPointer.set(pdbFile)

    def getPdb(self):
        return self._pdbPointer.get()

    def copyInfo(self, other):
        self._pdbPointer.copy(other._pdbPointer, copyId=False)


class SetOfPrincipalComponents(SetOfNormalModes):
    "Set of Principal Components is a SetOfNormalModes with a new name"
    pass


class FramesRange(CsvList):
    """ Store first and last frames in a movie. The last element will be
     the index in the stack of the first frame.
    """

    def __init__(self, firstFrame=1, lastFrame=0, firstFrameIndex=1):
        """
        Frames numbering always refer to the initial movies imported into
        the project. Counting goes from 1 to the number of frames.

        The underlying stack file could have all original frames or just
        a subset of them, so we need to keep which is the index of the
        first frame in the stack, again, the first image in a stack is 1
        """
        CsvList.__init__(self, pType=int)
        self.set([firstFrame, lastFrame, firstFrameIndex])

    def getFirstFrame(self):
        return self[0]

    def setFirstFrame(self, value):
        self[0] = value

    def getLastFrame(self):
        return self[1]

    def setLastFrame(self, value):
        self[1] = value

    def getRange(self):
        """ Return the frames range. """
        return self[0], self[1]

    def getFirstFrameIndex(self):
        return self[2]

    def setFirstFrameIndex(self, value):
        self[2] = value

    def rangeStr(self):
        return '[%d-%d]' % (self[0], self[1])

    def __str__(self):
        return '%s starts: %d' % (self.rangeStr(), self[2])


class Movie(Micrograph):
    """ Represent a set of frames of micrographs.
    """

    def __init__(self, location=None, **kwargs):
        Micrograph.__init__(self, location, **kwargs)

        # A movie could could have alignment information after some protocol
        self._alignment = None
        self._framesRange = FramesRange()

    def isCompressed(self):
        fn = self.getFileName()
        return fn.endswith('bz2') or fn.endswith('tbz')

    def getDim(self):
        """Return image dimensions as tuple: (Xdim, Ydim, Zdim)
        Consider compressed Movie files"""
        from pwem.emlib.image import ImageHandler

        if not self.isCompressed():
            x, y, z, n = ImageHandler().getDimensions(self)
            if x is not None:
                return x, y, max(z, n)
        return None

    def getNumberOfFrames(self):
        """ Return the number of frames of this movie
        """
        first, last, _ = self._framesRange
        if last > 0:
            return last - first + 1

        if not self.isCompressed():
            from pwem.emlib.image import ImageHandler

            x, y, z, n = ImageHandler().getDimensions(self)
            if x is not None:
                return max(z, n)  # Protect against evil mrc files
        return None

    def hasAlignment(self):
        return self._alignment is not None

    def getAlignment(self):
        return self._alignment

    def setAlignment(self, alignment):
        """Alignment are stored as a vector
        containing x and y coordinates. In this way 1 2 3 4
        are the data related with 2 frames with shifts (1,2)
        and (3,4)
        """
        self._alignment = alignment

    def getFramesRange(self):
        return self._framesRange

    def setFramesRange(self, value):
        self._framesRange.set(value)


class MovieAlignment(EMObject):
    """ Store the alignment between the different Movie frames.
    Also store the first and last frames used for alignment.
    """

    def __init__(self, first=-1, last=-1, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._first = Integer(first)
        self._last = Integer(last)
        self._xshifts = CsvList(pType=float)
        self._yshifts = CsvList(pType=float)
        self._xshifts.set(kwargs.get('xshifts', []))
        self._yshifts.set(kwargs.get('yshifts', []))
        # This list contain the coordinate where you begin the crop (x, y),
        # the width and height of the frames.
        # The order is: x,y, width and height.
        # For width and height, 0 means the entire frame.
        self._roi = CsvList(pType=int)

    def getRange(self):
        """ Return the first and last frames used for alignment.
        The first frame in a movie stack is numbered 1 and the maximum value
        for the last one is N,  where N is the total number of frames.
        """
        return self._first.get(), self._last.get()

    def getShifts(self):
        """ Return the list of alignment between one frame
        to another, from first to last frame used.
        """
        return self._xshifts, self._yshifts

    def setRoi(self, roiList):
        self._roi.set(roiList)

    def getRoi(self):
        """ Return the size used to align the movie
        """
        return self._roi


class SetOfMovies(SetOfMicrographsBase):
    """ Represents a set of Movies. """
    ITEM_TYPE = Movie

    def __init__(self, **kwargs):
        SetOfMicrographsBase.__init__(self, **kwargs)
        self._gainFile = String()
        self._darkFile = String()
        # Store the frames range to avoid loading the items
        self._firstFramesRange = FramesRange()

    def setGain(self, gain):
        self._gainFile.set(gain)

    def getGain(self):
        return self._gainFile.get()

    def setDark(self, dark):
        self._darkFile.set(dark)

    def getDark(self):
        return self._darkFile.get()

    def getFramesRange(self):
        return self._firstFramesRange

    def setFramesRange(self, value):
        self._firstFramesRange.set(value)

    def _setFirstDim(self, image):
        if self._firstDim.isEmpty():
            self._firstDim.set(image.getDim())
            self._firstFramesRange.set(image.getFramesRange())

    def _dimStr(self):
        """ Return the string representing the dimensions. """
        if self._firstDim.isEmpty():
            dimStr = 'No-Dim'
        else:
            try:
                x, y, z = self._firstDim
            except:
                return str(self._firstDim)
            first, last, i = self._firstFramesRange
            if last == 0:
                last = z
            dimStr = str(ImageDim(x, y, last - first + 1))

        return '%s %s' % (dimStr, self._firstFramesRange.rangeStr())

    def copyInfo(self, other):
        """ Copy SoM specific information plus inherited """
        SetOfMicrographsBase.copyInfo(self, other)
        self._gainFile.set(other.getGain())
        self._darkFile.set(other.getDark())
        self._firstFramesRange.set(other.getFramesRange())


class MovieParticle(Particle):
    def __init__(self, **kwargs):
        Particle.__init__(self, **kwargs)
        self._particleId = Integer()
        self._frameId = Integer()

    def getParticleId(self):
        return self._particleId.get()

    def setParticleId(self, partId):
        self._particleId.set(partId)

    def getFrameId(self):
        return self._frameId.get()

    def setFrameId(self, frameId):
        self._frameId.set(frameId)


class SetOfMovieParticles(SetOfParticles):
    """ This is just to distinguish the special case
    when the particles have been extracted from a set of movies.
    """
    ITEM_TYPE = MovieParticle


class FSC(EMObject):
    """Store a Fourier Shell Correlation"""

    def __init__(self, **kwargs):
        EMObject.__init__(self, **kwargs)
        self._x = CsvList(pType=float)
        self._y = CsvList(pType=float)

    def loadFromMd(self, mdObj, labelX, labelY):
        """
        Fill the x and y data of the FSC from a metadata.
        Params:
            mdObj: either a metadata object of a filename
            labelX: label used for frequency
            labelY: label used for FSC values
        """
        # iterate through x and y and create csvLists
        import pwem.emlib.metadata as md
        self._x.clear()
        self._y.clear()

        for row in md.iterRows(mdObj):
            self._x.append(row.getValue(labelX))
            self._y.append(row.getValue(labelY))

    def getData(self):
        return self._x, self._y

    def setData(self, xData, yData):
        self._x.set(xData)
        self._y.set(yData)

    def calculateResolution(self, threshold=0.143):
        """
        Calculate the FSC resolution value
        """
        dataLength = len(self._x)
        for i in range(dataLength):
            if float(self._y[i]) < threshold or i == dataLength - 1:
                above_res = float(self._x[i - 1])
                above_fsc = float(self._y[i - 1])
                below_res = float(self._x[i])
                below_fsc = float(self._y[i])
                break
        resolution = below_res - ((threshold - below_fsc) / (above_fsc - below_fsc) * (below_res - above_res))
        return "{0:.1f}".format(1 / resolution)


class SetOfFSCs(EMSet):
    """Represents a set of FSCs"""
    ITEM_TYPE = FSC


class SetOfData(EMSet):
    """Represents lines of data elements"""
    ITEM_TYPE = Object


class SetOfStats(SetOfData):
    """Replaces by SetOfData. Here for backwards compatibility. """
    pass
