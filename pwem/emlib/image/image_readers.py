import enum
from functools import lru_cache
from typing import Union

import numpy
import numpy as np
from PIL import Image
from tifffile import TiffFile, imread, imwrite
import mrcfile

import pwem.constants as emcts
from .. import lib
from scipy.ndimage import rotate, shift
from skimage.transform import rescale


import logging

logger = logging.getLogger(__name__)

class ROT_MODE(enum.Enum):
    FIXED=1  # Final image will have exactly the same dims as the input image in the SAME orientation
    NATURAL=2 # Rotation will require a bigger image to avoid loosing informatión in the corners
    CONDITIONAL=3 # Similar to FIXED, but shifting x and y original dimensions in some cases to reduce
                  # information loss. 45<rot<135 and 225<rot<315

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
            shape = images.shape
            if len(shape) == 2:
                images = [images]

        elif not isinstance(images, list):
            logger.warning("ImageStack initialized with an invalid type. Valid types are None, a singe numpy "
                           "array or a list of them. Current value is a %s. Continuing as an empty image." % type(images))
            images = []
            
        self._images = images
        self._properties = dict() if properties is None else properties

    def __iter__(self):
        return iter(self._images)

    def __len__(self):
        return len(self._images)

    def getImage(self, index=0, pilImage=False):
        if index >= len(self._images):
            raise IndexError("Image at %s position dos not exists. Current stack has %s images." % (index, len(self._images)))

        npImg = self._images[index]

        if pilImage:
            return self.asPilImage(npImg)
        else:
            return npImg

    def getCentralImage(self, pilImage=False):
        """ Returns the central image"""
        size = len(self._images)
        if size == 0:
            raise FileNotFoundError("Cannot get a central image. It may not exist or is not yet opened.")

        midIndex = size - 1 if size == 1 else (size // 2)
        return self.getImage(midIndex, pilImage=pilImage)

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

        if isinstance(imgStack, ImageStack):
            self._images.extend(imgStack.getImages())
        else:  # Numpy slice, validate is a slice?
            self._images.append(imgStack)

    ######### operations section ############
    @classmethod
    def normalizeSlice(cls, npImage):
        iMax = npImage.max()
        iMin = npImage.min()
        im255 = ((npImage - iMin) / (iMax - iMin) * 255).astype(numpy.uint8)
        return im255

    @classmethod
    def asPilImage(cls, npArray, normalize=True):
        """ Returns the npArray a numpy image
        :param npArray: 2d numpy array (image)
        :param normalize: by default it has to be normalized. Cancel this is you are sure it hase been normalized before"""

        if normalize:
            npArray = cls.normalizeSlice(npArray)

        return Image.fromarray(npArray)

    @classmethod
    def _center_crop(cls, npArray, target_height, target_width):

        h, w = npArray.shape

        start_y = (h - target_height) // 2
        start_x = (w - target_width) // 2

        return npArray[start_y:start_y + target_height, start_x:start_x + target_width]

    @classmethod
    def rotateSlice(cls, npArray: numpy.ndarray, angle: float, mode=ROT_MODE.FIXED, bg=None) -> numpy.ndarray:
        """Rotates a numpy array"""

        bg = npArray.mean() if bg is None else bg # Get the mean value
        reshape = mode != ROT_MODE.FIXED  # Fixed mode should not reshape the array

        # Rotate the image
        rotated = rotate(npArray, angle, reshape=reshape, mode='constant', cval=bg)

        # If mode
        if mode == ROT_MODE.CONDITIONAL:
            angle = angle % 360  # negative angles should turn into its equivalent: -15 -> 345
            target_height, target_width = npArray.shape

            # If in the region to shift dimension
            if (45 <= angle <= 135) or (225 <= angle <= 315):
                # Crop the image
                target_height, target_width = target_width, target_height

            rotated = cls._center_crop(rotated, target_height, target_width)


        return rotated

    @classmethod
    def shiftSlice(cls, image: numpy.ndarray, shifts: float, bg=None) -> numpy.ndarray:
        """Shifts a numpy array
        :param shifts = float or sequence. If a sequence, first value should be X shift and second Y shift
        """

        bg = image.mean() if bg is  None else bg# Get the mean value

        if not isinstance(shifts, float):
            # Swap: shift expect first element to be y abd then x. We have opposite convention
            shifts = (shifts[1], shifts[0])

        # Rotate the image
        return shift(image, shifts,  mode='constant', cval=bg)


    @classmethod
    def transformSlice(cls, npImage:numpy.ndarray, shifts: float, angle: float, mode=ROT_MODE.FIXED, bg=None):
        """ Apply the rotation and the shift to the npImage passed"""

        bg = npImage.mean() if bg is None else bg

        return cls.shiftSlice(cls.rotateSlice(npImage, angle, mode=mode, bg=bg), shifts, bg=bg)

    @classmethod
    def scaleSlice(cls, npImage, factors, anti_aliasing=True):
        """ Scales the npImage by the factor/s
        :param npImage: 2d numpy array
        :param factors: float or sequence
            The zoom factor along the axes. If a float, `zoom` is the same for each
            axis. If a sequence, `zoom` should contain one value for each axis.
        :param anti_aliasing:

        """
        return rescale(npImage, factors, anti_aliasing=anti_aliasing)

    @classmethod
    def thumbnailSlice(cls, npImage, width, height, normalize=True):
        """ Calls PIl thumbnail. It is less precise but faster than scaleSlice"""

        img = cls.asPilImage(npImage, normalize=normalize)
        img.thumbnail((width,height))
        return np.array(img)



    @classmethod
    def flipSlice(cls, npImage: numpy.ndarray, vertically=True):

        mode = 0 if vertically else 1
        return numpy.flip(npImage, mode)

    def flip(self, vertically=True):
        """Flip all images of an ImageStack horizontally or vertically.
            Vertically is up-down, horizontally is left-right."""

        return self._applyOperation(self.flipSlice, vertically)

    def flipV(self):
        """ Flips this stack vertically: up to down"""
        return self.flip()

    def flipH(self):
        """ Flips this stack horizontally: left to right"""
        return self.flip(False)

    def shift(self, shifts):
        """ Shifts the whole stack x and y returning a new stack.

        :param shift: The shift along the axes. If a float, shift is the same for each axis. If a sequence, shift should contain one value for each axis.

        """
        return self._applyOperation(self.shiftSlice, shifts)

    def scale(self, factors, anti_aliasing=True):
        """ Scales the stack by the factors
        :param: factors: Scale factors for spatial dimensions.
        """

        return self._applyOperation(self.scaleSlice, factors, anti_aliasing)

    def rotate(self, angle, mode=ROT_MODE.FIXED, bg=None):
        """rotates all its images the angle (deg) passed and returns a new ImageStack rotated"""
        return self._applyOperation(self.rotateSlice, angle, mode=mode, bg=bg)

    def transform(self, shifts, angle, mode=ROT_MODE.FIXED):
        """rotates all its images the angle (deg) passed and returns a new ImageStack rotated"""
        return self._applyOperation(self.transformSlice, shifts, angle, mode=mode)

    def multiply(self, factor: float):
        """ Multiplies the image stack by a factor
        :param: factor: to multiply values by it
        """

        return self._applyOperation(lambda npImage, factor: npImage*factor, factor)
    def invert(self):

        return self.multiply(factor=-1)

    def normalize(self):
        return self._applyOperation(self.normalizeSlice)

    def _applyOperation(self, operation, *args, **kwargs):
        rotImg = ImageStack()

        for image in self._images:
            rot_slice = operation(image, *args, **kwargs)
            rotImg.append(rot_slice)

        return rotImg

    def write(self, path):
        ImageReadersRegistry.write(self, path)




class ImageReader:
    @classmethod
    def open(cls, path):
        """ Opens the image in the path and returns a numpy array with the whole file (all slices included)"""
        raise NotImplementedError("Image reader %s does not implement 'open' method." % cls.__name__)

    @classmethod
    def canOpenSlices(cls):
        """ Returns True if the reader can open slices optimally. Without loading the whole file.
        If so, openSlice method should be implemented
        """
        return False

    @classmethod
    def openSlice(cls, path, slice):
        """ Opens a specific slice"""
        raise NotImplementedError("Image reader %s does not implement 'openSlice' method." % cls.__name__)

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
    def open(cls, filePath) -> ImageStack:
        """"Opens the file and returns and ImageStack.
        Accepts formats like 1@path/to/my/image.ext.
        In this case there can be a performance penalty since readers are reading all the
        stack but we are only taking one slice. This may be unavoidable in cases when you want
        to read a single image but not optimal when you are going to go through all the stack.
        """

        parts = filePath.split("@")

        filePath = parts[-1]

        # Get the reader that deals with the file extension.
        imageReader = cls.getReader(filePath)

        # If requesting a slice 1@ppath/to/image.ext
        if len(parts) == 2:

            sliceIndex = int(parts[0])

            if imageReader.canOpenSlices():
                data = imageReader.openSlice(filePath, sliceIndex)

            else:
                logger.debug("Requesting slice %s from %s. Suboptimal?." % (sliceIndex, filePath))
                data = imageReader.open(filePath)
                data = data[sliceIndex - 1]
        else:
            # Get the numpy array
            data = imageReader.open(filePath)

        return ImageStack(data)

    @classmethod
    def write(cls, imgStack: ImageStack, fileName: str, isStack=False) -> None:
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
    def open(cls, filePath: str):
        pilImg = Image.open(filePath)
        return numpy.array(pilImg)

    @classmethod
    def write(cls, imgStack: ImageStack, fileName: str, isStack=False) -> None:
        # So far write the first image in the stack
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
            key, path = path.split("@")

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
    @classmethod
    def open(cls, path):
        img = lib.Image()
        img.read(path)
        return img.getData()

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
    def canOpenSlices(cls):
        return True

    @classmethod
    def openSlice(cls, path, slice):
        """
        Reads a given image
           :param path (str) --> Image to be read
        """
        npImg = cls.open(path)
        return npImg[slice-1]

    @classmethod
    def open(cls, path: str):
        path = path.replace(":mrc", "")
        if "@" in path:
            path = path.split('@')[-1]

        mrcImg = cls.getMrcImage(path)
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
    def write(cls,
              imageStack: ImageStack,
              fileName: str,
              isStack: bool = False,
              samplingRate: Union[float, None] = None) -> None:
        """Generate a stack of images or a volume from a list of images."""
        sr = samplingRate if samplingRate else imageStack.getProperties().get("sr", 1.0)
        stack = numpy.stack(imageStack.getImages(), axis=0)

        with mrcfile.new(fileName, overwrite=True) as mrc:
            mrc.set_data(stack.astype(numpy.float32))
            if isStack:
                mrc.header.ispg = 0
            mrc.update_header_from_data()
            mrc.voxel_size = sr
        return True

    @classmethod
    def isMrcVolume(cls, mrcImg):
        if mrcImg.is_volume():
            return True
        return False

    @classmethod
    def isMrcStack(cls, mrcImg):
        if mrcImg.is_image_stack():
            return True
        return False

      
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
        """ Opens the full stack"""
        return cls.readAll(path)

    @classmethod
    def canOpenSlices(cls):
        return True

    @classmethod
    def openSlice(cls, path, slice):
        """
        Reads a given image
           :param filename (str) --> Image to be read
        """
        cls.stk_handler = open(path, "rb")
        cls.header_info = cls.readHeader()
        cls.IMG_BYTES = cls.FLOAT32_BYTES * cls.header_info["n_columns"] ** 2
        image = cls.readImage(slice - 1)
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
        return cls.readAll(filename)

    @classmethod
    def readAll(cls, filename):
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
