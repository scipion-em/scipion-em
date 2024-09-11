import logging
from functools import lru_cache

import numpy
from PIL import Image

from pwem.emlib.image.image_handler import ImageReadersRegistry
from metadataviewer.dao.numpy_dao import NumpyDao

logger = logging.getLogger(__name__)

import metadataviewer
from metadataviewer.model.renderers import ImageReader
from pwem.viewers.mdviewer.star_dao import StarFile
from pwem.viewers.mdviewer.sqlite_dao import ScipionSetsDAO


class ScipionImageReader(ImageReader):
    @classmethod
    def getCompatibleFileTypes(cls) -> list:
        return ImageReadersRegistry.getAvailableExtensions()

    @classmethod
    @lru_cache
    def open(cls, path):

        imgStack = ImageReadersRegistry.open(path)
        # So far readers are not returning the whole stack...even if they are tomograms or volumes.
        # They return in this case the central slice. So we always (for now) have a sinlge image.
        return cls._normalize(imgStack.getImage())

    @classmethod
    def _normalize(cls, npImage):
        iMax = npImage.max()
        iMin = npImage.min()
        im255 = ((npImage - iMin) / (iMax - iMin) * 255).astype(numpy.uint8)
        return Image.fromarray(im255)


def extendMDViewer(om: metadataviewer.model.ObjectManager):
    """ Function to extend the object manager with DAOs and readers"""
    om.registerDAO(ScipionSetsDAO)
    om.registerReader(ScipionImageReader)
    NumpyDao.addCompatibleFileType('cs')
    om.registerDAO(StarFile)
