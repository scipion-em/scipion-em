import unittest
from unittest.mock import Mock, patch, MagicMock
from pwem.protocols.protocol_particles import ProtExtractParticles
from pwem.tests.utils import getSoCTFsMock, getSoMMock


class TestExtractParticles(unittest.TestCase):

    def test_simpleLoadInputList(self):

        self.assertLoadInputList((1, 3),
                                 3, "Simple extraction of mics does not work.")

    def test_LoadInputListWithMicOthers(self):

        self.assertLoadInputList((1, 3),
                                 2, "Partial extraction of mics and others does not work.",
                                 micOthers=(2, 3))

        self.assertLoadInputList((1, 3),
                                 3, "Partial extraction of mics and others does not work.",
                                 micOthers=(1, 3))

        self.assertLoadInputList((1, 3),
                                 0, "Non overlapping extraction of mics and others does not work.",
                                 micOthers=(4, 9))

    def test_LoadInputListWithCTFs(self):

        self.assertLoadInputList((1, 3),
                                 2, "Partial extraction of mics and ctf does not work.",
                                 ctfs=(2, 3))

        self.assertLoadInputList((1, 3),
                                 3, "Extraction of same mics and ctfs does not work.",
                                 ctfs=(1, 3))

        self.assertLoadInputList((1, 3),
                                 0, "Extraction of non overlapping mics and ctfs does "
                                    "return mics and it shouldn't.",
                                 ctfs=(4, 8))

    def test_LoadInputListWithAll(self):

        # **** micOthers == mics
        #  ctf ==
        self.assertLoadInputList((1, 3),
                                 3, "Extraction of mics = others = ctfs does not work.",
                                 micOthers=(1, 3),
                                 ctfs=(1, 3))

        # ctf partial
        self.assertLoadInputList((1, 3),
                                 1, "Extraction of overlapping mics, others and partial ctfs does not work.",
                                 micOthers=(1, 3),
                                 ctfs=(3, 4))

        # ctf excluded
        self.assertLoadInputList((1, 3),
                                 0, "Extraction of overlapping mics, others and excluded ctfs does not work.",
                                 micOthers=(1, 3),
                                 ctfs=(4, 6))

        # **** micOthers partial overlap
        # ctf == mics
        self.assertLoadInputList((1, 3),
                                 1, "Extraction of mics with partial others and full ctfs does not work.",
                                 micOthers=(3, 6),
                                 ctfs=(1, 3))

        # ctf == overlaps mics
        self.assertLoadInputList((1, 3),
                                 1, "Extraction of mics, overlapping others and overlapping ctfs does not work.",
                                 micOthers=(3, 8),
                                 ctfs=(2, 5))

        # ctf != mics
        self.assertLoadInputList((1, 3),
                                 0, "Extraction of mics, overlapping others and non overlapping ctfs does not work.",
                                 micOthers=(3, 8),
                                 ctfs=(4, 5))

        # **** micOthers non overlap
        # ctf == mics
        self.assertLoadInputList((1, 3),
                                 0, "Extraction of mics, excluded others and full ctfs does not work.",
                                 micOthers=(4, 8),
                                 ctfs=(1, 3))

        # ctf overlaps mics
        self.assertLoadInputList((1, 3),
                                 0, "Extraction of mics, excluded others and full ctfs does not work.",
                                 micOthers=(4, 8),
                                 ctfs=(3, 8))

        # ctf != mics
        self.assertLoadInputList((1, 3),
                                 0, "Extraction of mics, != others and != ctfs does not work.",
                                 micOthers=(4, 8),
                                 ctfs=(9, 12))

    def assertLoadInputList(self, mics, expectedMicSize, msg, micOthers=None, ctfs=None):
        """ Asserts particle extraction using ctfs
            :parameter mics tuple with mic indexes  --> start and end
            :parameter ctfs tuple with ctf indexes  --> start and end
            """

        with patch('pwem.objects.SetOfCTF') as  setOfCtfPatcher:
            with patch('pwem.objects.SetOfMicrographs') as setOfMicPatcher:

                # If ctfs are passed
                if ctfs:
                    ctfs = getSoCTFsMock(ctfs[0], ctfs[1])

                    #  Make SoCtf return items
                    side_effects = [ctfs]

                    # Side effects are a way to mock what is instantiated, order is important
                    setOfCtfPatcher.side_effect = side_effects
                else:
                    setOfCtfPatcher = None

                if micOthers:
                    micOthers = getSoMMock(micOthers[0], micOthers[1])

                micDict = self.loadInputListRunner(setOfMicPatcher,
                                                   getSoMMock(mics[0], mics[1]),
                                                   othersMics=micOthers,
                                                   ctfsPatcher=setOfCtfPatcher)

                self.assertEqual(expectedMicSize, len(micDict), msg)

    def loadInputListRunner(self, setOfMicPatcher, coordSoMics, othersMics=None, ctfsPatcher=None, processedMics={}):
        """ :parameter setOfMicPatcher: Mocker object that patches SetOfMicrographs instantiation
            :parameter coordSoMics: SoM mock for the mics coming from the set of coordinates
            :parameter othersMics: Optional, SoM for testing extraction on other micrographs
            :parameter ctfsPatcher: Optional SetOfCTFs patcher to test extraction with ctf
            :parameter processedMics: already processed micDictionary
            """

        # Make SoM return items
        side_effects = [coordSoMics]

        # Side effects are a way to mock what is instantiated, order is important
        setOfMicPatcher.side_effect = side_effects

        # Mock : self._loadInputCoords
        # Luckily, _loadInputCoords receives micDict and fills micDict.
        # So, here we are making _loadInputCoords mock return the same param it receives
        def returnSame(value):
            return value

        with patch.object(ProtExtractParticles, '_loadInputCoords', side_effect=returnSame) as mock_method:

            # Instantiate
            extractParticles = ProtExtractParticles()

            # if testing extracting on others
            if othersMics:
                side_effects.append(othersMics)

                # Activate if _micOthers
                extractParticles._micsOther = Mock

                # Mock the pointer
                inputMicrographs = MagicMock()
                inputMicrographs.get.return_value = othersMics
                extractParticles.inputMicrographs = inputMicrographs

            # if testing ctf extraction
            if ctfsPatcher:
                # Activate if _micOthers
                extractParticles._useCTF = Mock

                # Mock the pointer
                ctfRelations = MagicMock()
                ctfRelations.get.return_value = MagicMock()
                extractParticles.ctfRelations = ctfRelations

            # Mock : self.inputCoordinates.get().getMicrographs()
            coordsMics = Mock()
            coordsMics.get.return_value.getMicrographs.return_value = coordSoMics
            extractParticles.inputCoordinates = coordsMics

            # Mock: _isStreamClosed
            extractParticles._isStreamClosed = Mock()
            # Mics already extracted (empty
            extractParticles.micDict = processedMics

            return extractParticles._loadInputList()
