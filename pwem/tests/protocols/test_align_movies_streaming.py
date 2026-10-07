import unittest
from unittest.mock import MagicMock, Mock

from pwem.protocols.protocol_align_movies import ProtAlignMovies, OUT_MOVIES


class TestAlignMoviesStreaming(unittest.TestCase):

    def test_FirstEmptyStreamCloseDoesNotRequireSummaryMovie(self):
        protocol = ProtAlignMovies()
        protocol._firstTimeOutput = True
        protocol.getAttributeValue = Mock(return_value=False)

        movieSet = MagicMock()
        movieSet.getIdSet.return_value = set()

        protocol._loadOutputSet = Mock(return_value=movieSet)
        protocol._updateOutputSet = Mock()
        protocol._storeSummary = Mock()
        protocol._defineTransformRelation = Mock()

        protocol.inputMovies = MagicMock()
        protocol.inputMovies.get.return_value.getDim.return_value = (64, 64, 7)

        streamMode = 1
        protocol._updateOutputMovieSet([], streamMode)

        protocol._storeSummary.assert_not_called()
        protocol._updateOutputSet.assert_called_once_with(
            OUT_MOVIES, movieSet, streamMode
        )
        protocol._defineTransformRelation.assert_called_once_with(
            protocol.inputMovies, movieSet
        )


if __name__ == "__main__":
    unittest.main()
