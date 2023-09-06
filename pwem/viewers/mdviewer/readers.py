
import logging
from functools import lru_cache

logger = logging.getLogger(__name__)

from metadataviewer.model.renderers import ImageReader
from PIL import Image
import mrcfile
import numpy as np


class MRCImageReader(ImageReader):
    @classmethod
    def open(cls, path):
        path = path.replace(":mrc", "")
        if not "@" in path:
            path ="1@"+path
        filePath = path.split('@')

        index = int(filePath[0])
        fileName = filePath[-1]
        mrcImg = cls.getMrcImage(fileName)
        if mrcImg.is_volume():
            dim = mrcImg.data.shape
            x = int(dim[0] /2)
            imfloat = mrcImg.data[x]
        elif mrcImg.is_image_stack():
            imfloat = mrcImg.data[index-1]
        else:
            imfloat = mrcImg.data

        iMax = imfloat.max()
        iMin = imfloat.min()
        im255 = ((imfloat - iMin) / (iMax - iMin) * 255).astype(np.uint8)
        img = Image.fromarray(im255)
        return img


    @classmethod
    @lru_cache
    def getMrcImage(cls, fileName):
        logger.info("Reading %s" % fileName)
        return mrcfile.open(fileName, permissive=True)

    @classmethod
    def getCompatibleFileTypes(cls) -> list:
        return ['mrc', 'mrc:mrc' 'mrcs', 'em', 'rec', 'ali']


class STKImageReader(ImageReader):
    IMG_BYTES = None
    stk_handler = None
    header_info = None
    HEADER_OFFSET = 1024
    FLOAT32_BYTES = 4
    TYPE = None

    @classmethod
    def open(cls, path):
        stk = path.split('@')
        if len(stk) > 1:
            image = cls.read(stk[-1], int(stk[0]))
            return image

    @classmethod
    def read(cls, filename, id):
        """
        Reads a given image
           :param filename (str) --> Image to be read
        """
        cls.stk_handler = open(filename, "rb")
        cls.header_info = cls.readHeader()
        cls.IMG_BYTES = cls.FLOAT32_BYTES * cls.header_info["n_columns"] ** 2
        image = cls.readImage(id - 1)
        iMax = image.max()
        iMin = image.min()
        image = ((image - iMin) / (iMax - iMin) * 255).astype('uint8')
        image = Image.fromarray(image)
        return image

    @classmethod
    def readHeader(cls):
        """
        Reads the header of the current file as a dictionary
            :returns The current header as a dictionary
        """
        header = cls.readNumpy(0, cls.HEADER_OFFSET)

        header = dict(img_size=int(header[1]), n_images=int(header[25]),
                      offset=int(header[21]),
                      n_rows=int(header[1]), n_columns=int(header[11]),
                      n_slices=int(header[0]),
                      sr=float(header[20]))

        cls.TYPE = "stack" if header["n_images"] > 1 else "volume"

        return header

    @classmethod
    def readNumpy(cls, start, end):
        """
        Read bytes between start and end as a Numpy array
            :param start (int) --> Start byte
            :param end (int) --> End byte
            :returns decoded bytes as Numpy array
        """
        return np.frombuffer(cls.readBinary(start, end), dtype=np.float32)

    @classmethod
    def readBinary(cls, start, end):
        """
        Read bytes between start and end
            :param start (int) --> Start byte
            :param end (int) --> End byte
            :returns the bytes read
        """
        cls.seek(start)
        return cls.stk_handler.read(end)

    @classmethod
    def readImage(cls, iid):
        """
        Reads a given image in the stack according to its ID
            :param iid (int) --> Image id to be read
            :returns Image as Numpy array
        """

        if cls.TYPE == "stack":
            start = 2 * cls.header_info["offset"] + iid * (
                    cls.IMG_BYTES + cls.header_info["offset"])
        else:
            start = cls.header_info["offset"] + iid * cls.IMG_BYTES

        img_size = cls.header_info["n_columns"]
        return cls.readNumpy(start, cls.IMG_BYTES).reshape([img_size, img_size])

    @classmethod
    def seek(cls, pos):
        """
        Move file pointer to a given position
            :param pos (int) --> Byte to move the pointer to
        """
        cls.stk_handler.seek(pos)

    @classmethod
    def getCompatibleFileTypes(cls) -> list:
        return ['stk', 'vol']
