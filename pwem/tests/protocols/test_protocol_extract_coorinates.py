import unittest
from unittest.mock import patch
from pwem.protocols import ProtExtractCoords
from pwem.tests.utils import getSoPartMock, getSoMMock, getMicNameFromId


class TestExtractCoordinates(unittest.TestCase):
    """ Tests extract coordinates protocol mocking some behaviour"""

    def test_extractCoordinatesById(self):
        """ Tests ProtExtractCoords.extractCoordinates method using mic id for matching"""

        output = self._extractCoordinatesMocker(sop=getSoPartMock(name="inParts"),
                                                som=getSoMMock(name="inMics"))

        self.assertEqual(3, output.getSize(), "Wrong coordinates extraction")

        # Test missing mics
        output = self._extractCoordinatesMocker(sop=getSoPartMock(name="inPartsMissing"),
                                                som=getSoMMock(start=2, end=4, name="inMicsMissing"))

        self.assertEqual(2, output.getSize(), "Wrong coordinates extraction when missing mic ids")

    def test_extractCoordinatesByMicName(self):
        """ Tests ProtExtractCoords.extractCoordinates method using mic name for matching"""

        # Mic id wil range from 5-7
        som = getSoMMock(5, 7, "noIdButName")

        # Set micnames from 1-3
        for key, mic in som.items():
            mic.setMicName(getMicNameFromId(key - 4))

        output = self._extractCoordinatesMocker(sop=getSoPartMock(name="inPartsNoIdButNames"),
                                                som=som)

        self.assertEqual(3, output.getSize(), "Wrong coordinates extraction, matching by micname fails")

    @staticmethod
    def _extractCoordinatesMocker(sop, som):
        extractionProt = ProtExtractCoords(workingDir="/tmp")

        # Patch getInputParticles
        with patch.object(ProtExtractCoords, 'getInputParticles',
                          return_value=sop) as mock_method:
            # Patch getInputMics
            with patch.object(ProtExtractCoords, 'getInputMicrographs',
                              return_value=som) as mock_method2:
                # Not in streaming
                extractionProt.streamingModeOn = False

                return extractionProt.extractCoordinates()


if __name__ == '__main__':
    unittest.main()
