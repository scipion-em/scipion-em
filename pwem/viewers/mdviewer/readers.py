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


def extendMDViewer(om: metadataviewer.model.ObjectManager):
    """ Function to extend the object manager with DAOs and readers"""
    om.registerDAO(ScipionSetsDAO)
    om.registerReader(ScipionImageReader)
    NumpyDao.addCompatibleFileType('cs')
    om.registerDAO(StarFile)
