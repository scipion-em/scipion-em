# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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


import numpy as np

from pwem.objects import (SetOfParticles, Particle, SetOfClasses, Volume, Transform, SetOfImages,
                          AtomStruct, EMSet, SetOfAtomStructs, EMObject, NO_INDEX)

from pyworkflow.object import String, CsvList, Float, Integer


class ParticleFlex(Particle):
    """Particle with flexibility information stored"""

    def __init__(self, progName="", **kwargs):
        Particle.__init__(self, **kwargs)
        self._flexInfo = FlexInfo(progName)
        self._zFlex = CsvList()
        self._zRed = CsvList()
        self._transform = Transform()

    def getFlexInfo(self):
        return self._flexInfo

    def setFlexInfo(self, flexInfo):
        self._flexInfo = flexInfo

    def getZFlex(self):
        return np.fromstring(self._zFlex.get(), sep=",")

    def setZFlex(self, zFlex):
        csvZFlex = CsvList()
        for c in zFlex:
            csvZFlex.append(c)
        self._zFlex = csvZFlex

    def getZRed(self):
        return np.fromstring(self._zRed.get(), sep=",")

    def setZRed(self, zRed):
        csvZRed = CsvList()
        for c in zRed:
            csvZRed.append(c)
        self._zRed = csvZRed

    def copyInfo(self, other):
        self.copy(other, copyId=False)


class SetOfParticlesFlex(SetOfParticles):
    """SetOfParticles with flexibility information stored"""
    ITEM_TYPE = ParticleFlex

    def __init__(self, progName="", **kwargs):
        SetOfParticles.__init__(self, **kwargs)
        self._flexInfo = FlexInfo(progName)

    def getFlexInfo(self):
        return self._flexInfo

    def setFlexInfo(self, flexInfo):
        self._flexInfo = flexInfo

    def copyInfo(self, other):
        super(SetOfParticles, self).copyInfo(other)
        if hasattr(other, "_flexInfo"):
            self._flexInfo.copyInfo(other._flexInfo)


class FlexInfo(EMObject):
    """ Object storing any information needed by other protocols/plugins to work properly such us
    neural network paths, basis degrees..."""

    def __init__(self, progName="", **kwargs):
        EMObject.__init__(self, **kwargs)
        self._progName = String(progName)

    def copyInfo(self, other):
        self.copy(other, copyId=False)

    def getProgName(self):
        return self._progName.get()

    def setProgName(self, progName):
        self._progName = String(progName)

    def setAttr(self, attrName, value):
        if isinstance(value, float):
            setattr(self, attrName, Float(value))
        elif isinstance(value, int):
            setattr(self, attrName, Integer(value))
        elif isinstance(value, str):
            setattr(self, attrName, String(value))
        elif isinstance(value, list) or isinstance(value, np.ndarray):
            csv = CsvList()
            for c in value:
                csv.append(c)
            setattr(self, attrName, csv)
        else:
            setattr(self, attrName, value)

    def getAttr(self, attrName):
        return getattr(self, attrName).get()


class VolumeFlex(Volume):
    """Volume with flexibility information stored"""

    def __init__(self, progName="", **kwargs):
        Volume.__init__(self, **kwargs)
        self._flexInfo = FlexInfo(progName=progName)
        self._zFlex = CsvList()
        self._zRed = CsvList()

    def getFlexInfo(self):
        return self._flexInfo

    def setFlexInfo(self, flexInfo):
        self._flexInfo = flexInfo

    def getZFlex(self):
        return np.fromstring(self._zFlex.get(), sep=",")

    def setZFlex(self, zFlex):
        csvZFlex = CsvList()
        for c in zFlex:
            csvZFlex.append(c)
        self._zFlex = csvZFlex

    def getZRed(self):
        return np.fromstring(self._zRed.get(), sep=",")

    def setZRed(self, zRed):
        csvZRed = CsvList()
        for c in zRed:
            csvZRed.append(c)
        self._zRed = csvZRed

    def copyInfo(self, other):
        self.copy(other, copyId=False)


class SetOfVolumesFlex(SetOfImages):
    """Represents a set of Volumes with flexibility information stored"""
    ITEM_TYPE = VolumeFlex
    REP_TYPE = VolumeFlex

    # Hint to GUI components to expose internal items for direct selection
    EXPOSE_ITEMS = True

    def __init__(self, progName="", **kwargs):
        SetOfImages.__init__(self, **kwargs)
        self._flexInfo = FlexInfo(progName=progName)

    def getFlexInfo(self):
        return self._flexInfo

    def setFlexInfo(self, flexInfo):
        self._flexInfo = flexInfo


class ClassFlex(SetOfParticlesFlex):
    """Class3D with flexibility information stored"""
    REP_TYPE = VolumeFlex
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


class SetOfClassesFlex(SetOfClasses):
    """ SetOfClasses3D with flexibility information stored"""
    ITEM_TYPE = ClassFlex
    REP_TYPE = VolumeFlex
    REP_SET_TYPE = SetOfVolumesFlex

    pass


class AtomStructFlex(AtomStruct):
    """ AtomStruct with flexibility information stored"""

    def __init__(self, progName="", **kwargs):
        AtomStruct.__init__(self, **kwargs)
        self._flexInfo = FlexInfo(progName)
        self._zFlex = CsvList()
        self._zRed = CsvList()

    def getFlexInfo(self):
        return self._flexInfo

    def setFlexInfo(self, flexInfo):
        self._flexInfo = flexInfo

    def getZFlex(self):
        return np.fromstring(self._zFlex.get(), sep=",")

    def setZFlex(self, zFlex):
        csvZFlex = CsvList()
        for c in zFlex:
            csvZFlex.append(c)
        self._zFlex = csvZFlex

    def getZRed(self):
        return np.fromstring(self._zRed.get(), sep=",")

    def setZRed(self, zRed):
        csvZRed = CsvList()
        for c in zRed:
            csvZRed.append(c)
        self._zRed = csvZRed

    def copyInfo(self, other):
        self.copy(other, copyId=False)

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
            _, filename = first
        elif t == int:
            _, filename = first, args[1]
        elif t == str:
            _, filename = NO_INDEX, first
        else:
            raise Exception('setLocation: unsupported type %s as input.' % t)

        self.setFileName(filename)

    def getSamplingRate(self):
        # This is only a way to use the default way to generate subsets
        return None

    def setSamplingRate(self, sr):
        # This is only a way to use the default way to generate subsets
        pass

    def setAlignment(self, sr):
        # This is only a way to use the default way to generate subsets
        pass


class SetOfAtomStructFlex(SetOfAtomStructs):
    """ Set containing AtomStructFlex items. """
    ITEM_TYPE = AtomStructFlex
    EXPOSE_ITEMS = True

    def __init__(self, progName="", **kwargs):
        EMSet.__init__(self, **kwargs)
        self._flexInfo = FlexInfo(progName)

    def getFlexInfo(self):
        return self._flexInfo

    def setFlexInfo(self, flexInfo):
        self._flexInfo = flexInfo

    def copyInfo(self, other):
        super(SetOfAtomStructFlex, self).copyInfo(other)
        if hasattr(other, "_flexInfo"):
            self._flexInfo.copyInfo(other._flexInfo)

    def getSamplingRate(self):
        # This is only a way to use the default way to generate subsets
        return None

    def setSamplingRate(self, sr):
        # This is only a way to use the default way to generate subsets
        pass

    def setAlignment(self, sr):
        # This is only a way to use the default way to generate subsets
        pass


class ClassStructFlex(SetOfParticlesFlex):
    """Class3D with flexibility information stored for structures"""
    REP_TYPE = AtomStructFlex
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

class SetOfClassesStructFlex(SetOfClasses):
    """ SetOfClasses3D with flexibility information stored for structures"""
    ITEM_TYPE = ClassStructFlex
    REP_TYPE = AtomStructFlex
    REP_SET_TYPE = SetOfAtomStructFlex

    pass