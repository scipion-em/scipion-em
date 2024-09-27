from functools import lru_cache

import numpy
from PIL import Image
from tifffile import TiffFile, imread, imwrite
import mrcfile

import pwem.constants as emcts
from .. import lib

import logging

logger = logging.getLogger(__name__)

# Classes to replace one day the functionality covered by ImageHandler... which uses xmipp binding
class ImageStack:
    """Class to hold image stacks. A single image is considered a stack of one image """
    def __init__(self, images=None, properties=None):
        """
        :param images: either None, an image as returned by the readers or a list of them.
             Images are numpy arrays
        :param properties: optional: dictionary of key value pairs for header information for those files tha may need it
        """
        if images is None:
            images = []
        elif isinstance(images, numpy.ndarray):
            images=[images]
        elif not isinstance(images, list):
            logger.warning("ImageStack initialized with an invalid type. Valid types are None, a singe numy array or a list of them. Current value is a %s. Continuing as an empty image." % type(images))
            images = []
            
        self._images = images
        self._properties = dict() if properties is None else properties

    def getImage(self, index=0, pilImage=False):
        if index >= len(self._images):
            raise IndexError("Image at %s position dos not exists. Current stack has %s images." % (index, len(self._images)))

        npImg = self._images[index]

        if pilImage:
            return self._normalize(npImg)
        else:
            return npImg

    def getImages(self):
        """Returns all the images"""
        return self._images

    def getProperty(self, property):
        """ Returns the property passed"""

        return self._properties.get(property, None)
    def getProperties(self):
        """Returns the properties dictionary"""
        return self._properties

    def append(self, imgStack):
        """ Appends to its local list of images the images inside the imgStack passed as parameter"""
        self._images.extend(imgStack.getImages())

    def _normalize(cls, npImage):
        iMax = npImage.max()
        iMin = npImage.min()
        im255 = ((npImage - iMin) / (iMax - iMin) * 255).astype(numpy.uint8)
        return Image.fromarray(im255)


class ImageReader:
    @staticmethod
    def getCompatibleExtensions() -> list:
        """ Returns a list of the compatible extensions the reader can handle"""
        pass

    @staticmethod
    def getDimensions(filePath):
        """ Returns the dimensions [X,Y,Z,N] of the file"""
        pass

    @staticmethod
    def write(images: ImageStack, fileName: str, isStack: bool) -> None:
        """ Generate a stack of images or a volume from a list of PIL images.
        :param images: An ImageStack instance with one or more images
        :param fileName: Path of the new stack
        :param isStack: Specifies whether to generate a volume or an image stack
        """
        logger.warning("write method not implemented. Cannot write %s" % fileName)


class ImageReadersRegistry:
    """ Class to register image readers to provide basic information about an image like dimensions or getting an image"""
    _readers = dict()  # Dictionary to hold the readers. The key is the extension

    @classmethod
    def addReader(cls, imageReader: ImageReader):

        for ext in imageReader.getCompatibleExtensions():
            ext_low = ext.lower()
            logger.debug("Adding %s as image reader for %s" % (imageReader, ext_low))
            cls._readers[ext_low] = imageReader

    @classmethod
    def getReader(cls, filePath) -> ImageReader:
        """ Returns the reader or None able to deal with filePath based on the extension."""
        ext = filePath.split(".")[-1].lower()
        logger.debug("Getting ImageReader for %s (%s)" % (filePath, ext))

        reader = cls._readers.get(ext, None)

        # Fall back to XmippImage reader
        if reader is None:
            logger.info("Reader not registered for %s files. Falling back to XmippReader" % ext)
            reader = XMIPPImageReader

        return reader

    @classmethod
    @lru_cache
    def open(cls, filePath):
        "Opens the file and returns and image stack"

        # Get the reader
        imageReader = cls.getReader(filePath)

        return ImageStack(imageReader.open(filePath))

    @classmethod
    def write(cls, imgStack:ImageStack, fileName: str, isStack=False) -> None:
        """Generate a stack of images from a list of PIL images."""

        imageWriter = cls.getReader(fileName)
        return imageWriter.write(imgStack, fileName, isStack)

    @classmethod
    def getAvailableExtensions(cls):
        """ Returns all the extensions it can handle"""
        return cls._readers.keys()



class PILImageReader(ImageReader):
    """ PIL image reader"""

    @staticmethod
    def getCompatibleExtensions() -> list:
        return ['png', 'jpg', 'jpeg']

    @staticmethod
    def getDimensions(filePath):
        im = Image.open(filePath)
        x, y = im.size  # (width,height) tuple
        return x, y, 1, 1

    @classmethod
    def open(cls,filePath:str):

        pilImg=Image.open(filePath)
        return numpy.array(pilImg)
    @classmethod
    def write(cls, imgStack: ImageStack, fileName: str, isStack=False) -> None:

        # So far write the first image in teh stack
        np_img = imgStack.getImage()
        im = Image.fromarray(numpy.uint8(np_img))
        im.save(fileName)

        return True

class TiffImageReader(ImageReader):
    """ Tiff image reader"""

    @staticmethod
    def getCompatibleExtensions() -> list:
        return ['tif', 'tiff', 'gain', 'eer']

    @staticmethod
    def getDimensions(filePath):
        tif = TiffFile(filePath)
        frames = len(tif.pages)  # number of pages in the file
        page = tif.pages[0]  # get shape and dtype of the image in the first page
        x, y = page.imagewidth, page.imagelength  # IMPORTANT: to match xmipp convention

        return x, y, 1, frames

    @classmethod
    def open(cls, path: str):

        key = 0
        if "@" in path:
            key, path=path.split("@")

        npImg = imread(path, key=key)
        return npImg

    @classmethod
    def write(cls, imgStack: ImageStack, fileName: str, isStack=False) -> None:

        npImg = imgStack.getImage().astype("uint8")
        imwrite(fileName, npImg)
        return True

class EMANImageReader(ImageReader):
    """ Image reader for eman file formats"""

    @staticmethod
    def getCompatibleExtensions() -> list:
        return ["img"]

    @staticmethod
    def getDimensions(filePath):
        from pwem import Domain
        getImageDimensions = Domain.importFromPlugin(
            'eman2.convert', 'getImageDimensions',
            doRaise=True)
        return getImageDimensions(filePath)  # we are ignoring index here


class XMIPPImageReader(ImageReader):
    @staticmethod
    def getCompatibleExtensions():
        return emcts.ALL_MRC_EXTENSIONS + emcts.ALL_TIF_EXTENSIONS + ["hdf5", "dm4", "stk", "spi", "vol", "tif", "em", "map"]

    @staticmethod
    def getDimensions(filePath):
        img = lib.Image()
        img.read(filePath, lib.HEADER)
        return img.getDimensions()


class MRCImageReader(ImageReader):
    """ Image reader for MRC files"""

    @staticmethod
    def getCompatibleExtensions() -> list:
        return emcts.ALL_MRC_EXTENSIONS

    @staticmethod
    def getDimensions(filePath):
        from pwem.convert import headers
        header = headers.Ccp4Header(filePath, readHeader=True)
        return header.getXYZN()

    @classmethod
    def open(cls, path: str):
        isVol = path.endswith(":mrc")
        forceIndex = False
        path = path.replace(":mrc", "")
        if not "@" in path:
            path = "1@" + path
            forceIndex = True
        filePath = path.split('@')

        index = int(filePath[0])
        fileName = filePath[-1]
        mrcImg = cls.getMrcImage(fileName)
        if mrcImg.is_volume() or isVol:
            dim = mrcImg.data.shape
            x = int(dim[0] / 2) if forceIndex or index == 0 else index
            imfloat = mrcImg.data[x-1, :, :]
        elif mrcImg.is_image_stack():
            imfloat = mrcImg.data[index - 1]
        else:
            imfloat = mrcImg.data

        return imfloat

    @classmethod
    @lru_cache
    def getMrcImage(cls, fileName):
        logger.info("Reading %s" % fileName)
        return mrcfile.mmap(fileName, mode='r+', permissive=True)

    @classmethod
    def getArray(cls, filename):
        filename = filename.split('@')[-1]
        with mrcfile.open(filename, permissive=True) as mrc:
            return numpy.array(mrc.data)

    @classmethod
    def write(cls, imageStack: ImageStack, fileName: str, isStack=False) -> None:
        """Generate a stack of images or a volume from a list of PIL images."""
        sr = imageStack.getProperties().get("sr", 1.0)
        stack = numpy.stack(imageStack.getImages(), axis=0)

        with mrcfile.new(fileName, overwrite=True) as mrc:
            mrc.set_data(stack.astype(numpy.float32))
            if isStack:
                mrc.header.ispg = 0
            mrc.update_header_from_data()
            mrc.voxel_size = sr
        return True


class STKImageReader(ImageReader):
    IMG_BYTES = None
    stk_handler = None
    header_info = None
    HEADER_OFFSET = 1024
    FLOAT32_BYTES = 4
    TYPE = None

    @classmethod
    def __init__(cls, fileName):
        cls.stk_handler = open(fileName, "rb")
        cls.header_info = cls.readHeader()
        cls.IMG_BYTES = cls.FLOAT32_BYTES * cls.header_info["n_columns"] ** 2

    @classmethod
    def open(cls, path):
        stk = path.split('@')
        id = int(stk[0]) if len(stk) > 1 else 1
        image = cls.read(stk[-1], id)
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
        return image

    @staticmethod
    def getDimensions(filePath):

        STKImageReader.stk_handler = open(filePath, "rb")
        STKImageReader.header_info = STKImageReader.readHeader()
        STKImageReader.IMG_BYTES = STKImageReader.FLOAT32_BYTES * STKImageReader.header_info["n_columns"] ** 2
        header = STKImageReader.header_info
        return (header['n_rows'], header['n_columns'], header['n_slices'],
                header['n_images'])

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
        return numpy.frombuffer(cls.readBinary(start, end), dtype=numpy.float32)

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
    def getCompatibleExtensions(cls) -> list:
        return ['stk', 'vol']

    @classmethod
    def getArray(cls, filename):
        filename = filename.split('@')[-1]
        dimX, dimY, dimZ, nImages = cls.getDimensions(filename)
        if cls.TYPE == 'volume':
            numpyStack = numpy.stack([cls.readImage(ii) for ii in range(0, dimZ, 1)])
        else:
            numpyStack = numpy.stack([cls.readImage(ii) for ii in range(0, nImages, 1)])

        return numpyStack



# Register reader in the registry. Latest registered will take priority.
ImageReadersRegistry.addReader(XMIPPImageReader)
ImageReadersRegistry.addReader(MRCImageReader)
ImageReadersRegistry.addReader(STKImageReader)
ImageReadersRegistry.addReader(EMANImageReader)
ImageReadersRegistry.addReader(PILImageReader)
ImageReadersRegistry.addReader(TiffImageReader)
