import unittest
from unittest.mock import patch

from pwem import Plugin, findFolderWithPattern, CUDA_LIB_VAR, NO_VERSION_FOUND_STR


class TestPluginCudaGuess(unittest.TestCase):
    """ Tests code that tries to guess what is the CUDA version and functionality related"""

    def test_cuda_internal_methods(self):

        # findFolderWithPattern
        self.assertIsNone(findFolderWithPattern("/one/two/three", "four"))
        self.assertEqual(findFolderWithPattern("/one/two/three", "tw"), "two")

        v = Plugin.getVersionFromPath("/usr/local/cuda-10.2")
        self.assertEqual(v.major, 10, "Major not parsed from path with version")
        self.assertEqual(v.minor, 2, "Minor not parsed from path with version")

        v = Plugin.getVersionFromPath("/usr/local/cuda10.2")
        self.assertEqual(str(v), NO_VERSION_FOUND_STR, "Default value not working")

    def test_cuda_version_in_variable(self):

        with patch("os.path.realpath") as realpath_faked:
            realpath_faked.return_value = "/usr/local/cuda-13.8/lib-12.4"

            with patch.object(Plugin, "getVar") as getVarFaked:

                # Mock getVar return value
                getVarFaked.return_value = "/usr/local/useless-11.6"

                v = Plugin.getVersionFromVariable(CUDA_LIB_VAR)
                self.assertEqual(v.major, 12, "Major not parsed from link path")
                self.assertEqual(v.minor, 4, "Minor not parsed from link path")

                v = Plugin.getVersionFromVariable(CUDA_LIB_VAR, pattern="cuda")
                self.assertEqual(v.major, 13, "Major not parsed from link path with a pattern")
                self.assertEqual(v.minor, 8, "Minor not parsed from link path with a pattern")

                v = Plugin.guessCudaVersion(CUDA_LIB_VAR)
                self.assertEqual(v.major, 13, "Major not parsed from link path with a pattern")
                self.assertEqual(v.minor, 8, "Minor not parsed from link path with a pattern")

                # Mock getVar return value. Grigory's case ;-)
                realpath_faked.return_value = "/usr/lib/x86_64-linux-gnu"

                v = Plugin.guessCudaVersion(CUDA_LIB_VAR)
                self.assertEqual(str(v), "10.1", "Guess cuda version does not return default=10.1")
