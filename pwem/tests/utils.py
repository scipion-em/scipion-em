from unittest.mock import MagicMock

from pwem.objects import CTFModel, Micrograph

def getMicList(start=1, end = 3):
    """ Mocks a mic dict"""

    newMicDict = {}

    for index in range(start, end+1):
        mic = Micrograph()
        mic.setMicName("mic%s" % index)
        newMicDict[mic.getMicName()]= mic

    return newMicDict

def getSoMMock(start=1, end=3):
    """ Mocks a set of micrographs"""
    som = MagicMock()
    som.__iter__.return_value = getMicList(start, end).values()
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