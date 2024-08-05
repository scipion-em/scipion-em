import logging
from functools import lru_cache
from pwem.emlib.image.image_handler import ImageReadersRegistry
from metadataviewer.dao.numpy_dao import NumpyDao

logger = logging.getLogger(__name__)

import metadataviewer
from metadataviewer.model import Table
from metadataviewer.model.renderers import ImageReader
from pwem.viewers.mdviewer.star_dao import StarFile
from pwem.viewers.mdviewer.sqlite_dao import ScipionSetsDAO


class ScipionImageReader(ImageReader):
    @classmethod
    def getCompatibleFileTypes(cls) -> list:
        return list(ImageReadersRegistry._readers.keys())

    @classmethod
    @lru_cache
    def open(cls, path):
        ext = path.split('.')[-1]
        imageReader = ImageReadersRegistry._readers[ext]
        return imageReader.open(path)

    @classmethod
    def write(cls, pilImages: list, fileName: str, sr=1.0, isStack=False) -> None:
        """Generate a stack of images from a list of PIL images."""
        ext = fileName.split('.')[-1]
        imageReader = ImageReadersRegistry._readers[ext]
        imageReader.write(pilImages, fileName, sr, isStack)


def extendMDViewer(om: metadataviewer.model.ObjectManager):
    """ Function to extend the object manager with DAOs and readers"""
    om.registerDAO(ScipionSetsDAO)
    om.registerReader(ScipionImageReader)
    NumpyDao.addCompatibleFileType('cs')
    om.registerDAO(StarFile)
