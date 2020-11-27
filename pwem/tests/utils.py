from unittest.mock import MagicMock

from pwem.objects import CTFModel, Micrograph, Particle, Coordinate


def getMicList(start=1, end=3):
    """ Mocks a mic dict"""

    newMicDict = {}

    for index in range(start, end + 1):
        mic = Micrograph()
        mic.setMicName("mic%s" % index)
        newMicDict[index] = mic

    return newMicDict


def getSoMMock(start=1, end=3, name=None):
    """ Mocks a set of micrographs"""
    som = MagicMock(name=name)
    micList = getMicList(start, end)
    som.__iter__.return_value = micList.values()
    som.iterItems.return_value = list(micList.values())
    som.items.side_effect = micList.items
    som.__getitem__.side_effect = lambda key: micList.__getitem__(key)
    return som


def getCTFList(start=1, end=3):
    """ Mocks a CTF dict"""
    newCtfDict = {}

    for index in range(start, end + 1):
        ctf = CTFModel()
        mic = Micrograph()
        mic.setMicName("mic%s" % index)
        ctf.setMicrograph(mic)
        newCtfDict[mic.getMicName()] = ctf

    return newCtfDict


def getSoCTFsMock(start=1, end=3):
    """ Mocks a set of CTFs"""
    som = MagicMock()
    som.__iter__.return_value = getCTFList(start, end).values()
    return som


def getParticleList(start=1, end=3):
    """ Mocks a Particle dict"""
    newPartDict = {}

    for index in range(start, end + 1):
        part = Particle()
        coord = Coordinate()
        coord.setMicId(index)
        coord.setMicName(getMicNameFromId(index))
        part.setCoordinate(coord)
        newPartDict[index] = part

    return newPartDict


def getMicNameFromId(id):
    return "mic%s" % id


def getSoPartMock(start=1, end=3, name=None):
    """ Mocks a set of Particles"""
    som = MagicMock(name=name)
    som.__iter__.return_value = getParticleList(start, end).values()
    return som


def get2DCoordList(start=1, end=3):
    """ Mocks a 2d coordinates dict"""

    newCoordDict = {}

    for index in range(start, end + 1):
        coord = Coordinate()
        newCoordDict[index] = coord

    return newCoordDict


def getSo2DCMock(start=1, end=3):
    """ Mocks a set of 2D coordinates"""
    som = MagicMock()
    som.__iter__.return_value = get2DCoordList(start, end).values()
    return som
