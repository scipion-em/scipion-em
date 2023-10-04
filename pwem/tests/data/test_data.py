"""
Created on May 20, 2013

@author: laura
"""

from glob import iglob
import sqlite3
from unittest import TestCase

import numpy as np

from pyworkflow import SCIPION_DEBUG_NOCLEAN
from pyworkflow.tests import *
import pyworkflow.utils as pwutils

import pwem.objects as emobj
from pwem import emlib
from pwem.emlib import metadata as md
import pwem.protocols as emprot
import pyworkflow.tests as pwtests
from pwem.convert import SequenceHandler

# set to true if you want to check how fast is the access to
# the database
SPEEDTEST = True


def getIndex(fileName):
    """ Returns list of indexes for the table 'objects'
    in the database fileName"""
    conn = sqlite3.connect(fileName)
    cur = conn.cursor()
    cur.execute("PRAGMA index_list(objects);")
    rows = cur.fetchall()
    conn.close()
    return rows  # list with index names


def dropIndex(fileName, indexName):
    """ Drop index called indexName for table objects
    in database fileName"""
    conn = sqlite3.connect(fileName)
    cur = conn.cursor()
    cur.execute("DROP INDEX %s" % indexName)
    conn.close()


def createDummyProtocol(projName):
    """create a dummy protocol, returns protocol object"""
    proj = Manager().createProject(projName)
    os.chdir(proj.path)

    from pwem.protocols import EMProtocol

    prot = proj.newProtocol(EMProtocol)
    prot.setObjLabel('dummy protocol')
    try:
        proj.launchProtocol(prot)
    except Exception as e:
        logger.error("Can't launch EMProtocol", exc_info=e)
    return prot


class TestFSC(unittest.TestCase):
    _labels = [SMALL, WEEKLY]

    def testIO(self):
        """Test basic FSC object"""
        xList = [0.00, 0.05, 0.10, 0.15, 0.2]
        yList = [1.00, 0.95, 0.90, 0.85, 0.2]
        fsc = emobj.FSC()
        fsc.setData(xList, yList)
        x, y = fsc.getData()
        self.assertEqual(xList, x)
        self.assertEqual(yList, y)

    def testMd(self):
        """test create FSC from metadata"""
        xList = [0.00, 0.05, 0.10, 0.15, 0.2]
        yList = [1.00, 0.95, 0.90, 0.85, 0.2]
        md1 = emlib.MetaData()
        for freq, fscValue in zip(xList, yList):
            id = md1.addObject()
            md1.setValue(emlib.MDL_RESOLUTION_FREQ, freq, id)
            md1.setValue(emlib.MDL_RESOLUTION_FRC, fscValue, id)
        fsc = emobj.FSC()
        fsc.loadFromMd(md1, emlib.MDL_RESOLUTION_FREQ,
                       emlib.MDL_RESOLUTION_FRC)
        x, y = fsc.getData()
        self.assertEqual(xList, x)
        self.assertEqual(yList, y)


class TestImage(unittest.TestCase):
    _labels = [SMALL, WEEKLY]

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.mic1 = cls.dataset.getFile('mic1')

    def testLocation(self):
        fn = self.mic1
        mic = emobj.Micrograph()
        mic.setFileName(fn)

        # Check that setFileName-getFileName is working properly
        self.assertEqual(fn, mic.getFileName())

        # Check the location is accepted from constructor
        mic2 = emobj.Micrograph(fn)
        self.assertEqual(fn, mic2.getFileName())

        volStk = '/data/current/volumes/all_volumes.stk'
        vol1 = emobj.Volume((1, volStk))
        self.assertEqual(1, vol1.getIndex())
        self.assertEqual(volStk, vol1.getFileName())

        self.assertEqual('all_volumes.stk', vol1.getBaseName())

    def testOrigin(self):
        fn = self.mic1
        mic = emobj.Micrograph()
        mic.setFileName(fn)
        mic.setSamplingRate(1.)

        referenceT = emobj.Transform()
        referenceT.setShifts(-4608., -4720., 1)
        mreferenceT = referenceT.getMatrix()

        t = mic.getOrigin(True)
        mt = t.getMatrix()
        for i in range(4):
            for j in range(4):
                self.assertAlmostEquals(mt[i][j], mreferenceT[i][j], delta=0.5)


class TestImageHandler(unittest.TestCase):
    _labels = [SMALL, WEEKLY]

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.dsFormat = DataSet.getDataSet('movies')
        cls.dsRelion = DataSet.getDataSet('relion_tutorial')

    def testExistLocation(self):
        volFn = self.dataset.getFile('volumes/volume_1_iter_002.mrc')

        ih = emlib.image.ImageHandler()
        # Test the volume filename exists
        self.assertTrue(ih.existsLocation(volFn))
        # Test missing filename
        self.assertFalse(ih.existsLocation(volFn.replace('.mrc', '_fake.mrc')))
        # Test the :mrc is append when used as volume
        newFn = ih.getVolFileName(volFn)
        self.assertEqual(newFn, volFn + ":mrc")
        # Test that the new filename still exists even with the :mrc suffix
        self.assertTrue(ih.existsLocation(newFn))

    def testGetDimensions(self):
        volLabel = ':mrc'
        movieLabel = ':mrcs'
        ih = emlib.image.ImageHandler()

        # MICROGRAPH
        expectedSize_Mic = [9216, 9441, 1, 1]
        micFn = self.dataset.getFile('micrographs/BPV_1386.mrc')
        self.assertTrue(ih.existsLocation(micFn))
        X, Y, Z, N = ih.getDimensions(micFn)
        self.assertEqual([X, Y, Z, N], expectedSize_Mic)
        # Labelled as movie
        X, Y, Z, N = ih.getDimensions(micFn + movieLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Mic)
        # Labelled as volume
        X, Y, Z, N = ih.getDimensions(micFn + volLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Mic)

        # MOVIE .MRC
        expectedSize_Mov = [4096, 4096, 1, 7]
        expectedSize_Vol = [4096, 4096, 7, 1]
        movFn = self.dsFormat.getFile('qbeta/qbeta.mrc')  # Not labelled, so treated as volume
        self.assertTrue(ih.existsLocation(movFn))
        X, Y, Z, N = ih.getDimensions(movFn)
        self.assertEqual([X, Y, Z, N], expectedSize_Vol)
        # Labelled as movie
        X, Y, Z, N = ih.getDimensions(movFn + movieLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Mov)
        # Labelled as volume
        X, Y, Z, N = ih.getDimensions(movFn + volLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Vol)

        # MOVIE .MRCS
        expectedSize_Mov = [1950, 1950, 1, 16]
        expectedSize_Vol = [1950, 1950, 16, 1]
        movFn = self.dsFormat.getFile('Falcon_2012_06_12-16_55_40_0_movie.mrcs')
        self.assertTrue(ih.existsLocation(movFn))
        X, Y, Z, N = ih.getDimensions(movFn)
        self.assertEqual([X, Y, Z, N], expectedSize_Mov)
        # Labelled as movie
        X, Y, Z, N = ih.getDimensions(movFn + movieLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Mov)
        # Labelled as volume
        X, Y, Z, N = ih.getDimensions(movFn + volLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Vol)

        # VOLUME
        expectedSize_Mov = [64, 64, 1, 64]
        expectedSize_Vol = [64, 64, 64, 1]
        volFn = self.dataset.getFile('volumes/volume_1_iter_002.mrc')
        self.assertTrue(ih.existsLocation(volFn))
        X, Y, Z, N = ih.getDimensions(volFn)
        self.assertEqual([X, Y, Z, N], expectedSize_Vol)
        # Labelled as movie
        X, Y, Z, N = ih.getDimensions(volFn + movieLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Mov)
        # Labelled as volume
        X, Y, Z, N = ih.getDimensions(volFn + volLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Vol)

        # STACK OF VOLUMES
        expectedSize_Mov = [60, 60, 1, 180]
        expectedSize_Svol = [60, 60, 60, 3]
        sVolFn = self.dsRelion.getFile('import/case2/relion_volumes.mrc')
        self.assertTrue(ih.existsLocation(sVolFn))
        X, Y, Z, N = ih.getDimensions(sVolFn)
        self.assertEqual([X, Y, Z, N], expectedSize_Svol)
        # Labelled as movie
        X, Y, Z, N = ih.getDimensions(sVolFn + movieLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Mov)
        # Labelled as volume
        X, Y, Z, N = ih.getDimensions(sVolFn + volLabel)
        self.assertEqual([X, Y, Z, N], expectedSize_Svol)

        # eer format (In case we want to cancel the warning)
        import logging
        tifffileLogger = logging.getLogger("tifffile.tifffile")
        tifffileLogger.disabled = True

        eerFile = self.dsFormat.getFile("eer")
        logger.info("Reading dimension for %s" % eerFile )
        X, Y, Z, N = ih.getDimensions(eerFile)
        self.assertEqual([X, Y, Z, N], [4096, 4096, 1, 567])

    def test_convertMicrographs(self):
        """ Convert micrographs to different formats.
         EMAN2 required for .img
        """
        micFn = self.dataset.getFile('micrographs/BPV_1386.mrc')
        outSuffix = pwutils.replaceBaseExt(micFn, 'img')
        ih = emlib.image.ImageHandler()

        outFn = join('/tmp', outSuffix)
        logger.info("Converting: \n%s -> %s" % (micFn, outFn))

        ih.convert(micFn, outFn)

        self.assertTrue(os.path.exists(outFn))
        self.assertTrue(pwutils.getFileSize(outFn) > 0)

        pwutils.cleanPath(outFn)
        pwutils.cleanPath(outFn.replace('.img', '.hed'))

    def test_readDM4(self):
        """ Check we can read dm4 files (using EMAN)
        """
        micFn = self.dsFormat.getFile('SuperRef_c3-adp-se-xyz-0228_001.dm4')

        ih = emlib.image.ImageHandler()
        # Check that we can read the dimensions of the dm4 file:
        EXPECTED_SIZE = (7676, 7420, 1, 1)
        self.assertEqual(ih.getDimensions(micFn), EXPECTED_SIZE)

        # We could even convert to an mrc file:
        outSuffix = pwutils.replaceBaseExt(micFn, 'mrc')

        outFn = join('/tmp', outSuffix)
        logger.info("Converting: \n%s -> %s" % (micFn, outFn))

        ih.convert(micFn, outFn)

        self.assertTrue(os.path.exists(outFn))
        self.assertTrue(pwutils.getFileSize(outFn) > 0)
        # Check dimensions are still the same:
        self.assertEqual(ih.getDimensions(outFn), EXPECTED_SIZE)

        # Clean up tmp files
        pwutils.cleanPath(outFn)

    def test_readCompressedTIF(self):
        """ Check we can read tif files
        """
        micFn = self.dsFormat.getFile('c3-adp-se-xyz-0228_200.tif')

        logger.info("Reading and converting %s" % micFn)
        ih = emlib.image.ImageHandler()
        # Check that we can read the dimensions of the dm4 file:
        EXPECTED_SIZE = (7676, 7420, 1, 38)
        self.assertEqual(ih.getDimensions(micFn), EXPECTED_SIZE)

        # We could even convert to an mrc file:
        outSuffix = pwutils.replaceBaseExt(micFn, 'mrc')

        outFn = join('/tmp', outSuffix)
        logger.info("Converting: \n%s -> %s" % ((1, micFn), outFn))

        ih.convert((1, micFn), outFn)

        self.assertTrue(os.path.exists(outFn))
        self.assertTrue(pwutils.getFileSize(outFn) > 0)
        self.assertEqual(ih.getDimensions(outFn), (7676, 7420, 1, 1))

        # Clean up tmp files
        pwutils.cleanPath(outFn)

    def test_convertMovie(self):
        """Check movie conversion"""
        movFn = self.dsFormat.getFile('qbeta/qbeta.mrc') + ":mrcs"

        ih = emlib.image.ImageHandler()
        # Check that we can read the dimensions of the dm4 file:
        EXPECTED_SIZE = (4096, 4096, 1, 7)
        EXPECTED_DT = emlib.DT_USHORT

        self.assertEqual(ih.getDimensions(movFn), EXPECTED_SIZE)
        self.assertEqual(ih.getDataType(movFn), EXPECTED_DT)

        outFn = join('/tmp/qbeta_converted.mrcs')

        ih.convertStack(movFn, outFn, 2, 6)

        self.assertTrue(os.path.exists(outFn))
        self.assertTrue(pwutils.getFileSize(outFn) > 0)
        self.assertEqual(ih.getDimensions(outFn), (4096, 4096, 1, 5))
        self.assertEqual(ih.getDataType(outFn), EXPECTED_DT)

        if pwutils.envVarOn(SCIPION_DEBUG_NOCLEAN):
            logger.info("Not cleaning output movie: ", outFn)
        else:
            pwutils.cleanPath(outFn)

    def test_truncateMask(self):
        ih = emlib.image.ImageHandler()

        maskFn = self.dataset.getFile('masks/mask.vol')
        outFn = join('/tmp', 'mask.vol')

        ih.truncateMask(maskFn, outFn, newDim=128)
        EXPECTED_SIZE = (128, 128, 128, 1)

        self.assertTrue(os.path.exists(outFn))
        self.assertTrue(pwutils.getFileSize(outFn) > 0)
        self.assertEqual(ih.getDimensions(outFn), EXPECTED_SIZE)

        pwutils.cleanPath(outFn)

    def test_createEmptyImage(self):
        outFn = join('/tmp', 'empty.mrc')
        SIZE = (128, 128, 1, 1)
        ih = emlib.image.ImageHandler()
        DT = emlib.DT_FLOAT
        ih.createEmptyImage(outFn, SIZE[0], SIZE[1], dataType=DT)

        self.assertTrue(pwutils.getFileSize(outFn) > 0)
        self.assertEqual(ih.getDimensions(outFn), SIZE)
        self.assertEqual(ih.getDataType(outFn), DT)

    def test_scaleStack(self):
        particles = self.dataset.getFile("particles/BPV_1386.stk")
        outFn = join('/tmp', 'scaled.mrc')
        ih = emlib.image.ImageHandler()
        DT = emlib.DT_FLOAT

        # Scaled with a higher dimension using finalDimension parameter
        EXPECTED_SIZE = (256, 256, 131, 1)
        ih.scale2DStack(particles, outFn, finalDimension=EXPECTED_SIZE[0])

        self.assertTrue(os.path.exists(outFn))
        self.assertTrue(pwutils.getFileSize(outFn) > 0)
        self.assertEqual(ih.getDimensions(outFn), EXPECTED_SIZE)
        self.assertEqual(ih.getDataType(outFn), DT)

        # Scaled with a lower dimension using finalDimension parameter
        EXPECTED_SIZE = (64, 64, 131, 1)
        ih.scale2DStack(particles, outFn, finalDimension=EXPECTED_SIZE[0])

        self.assertTrue(os.path.exists(outFn))
        self.assertTrue(pwutils.getFileSize(outFn) > 0)
        self.assertEqual(ih.getDimensions(outFn), EXPECTED_SIZE)
        self.assertEqual(ih.getDataType(outFn), DT)

        # Scaled with a higher dimension using scaleFactor parameter
        scaleFactor = 1.5
        EXPECTED_SIZE = (210, 210, 131, 1)
        ih.scale2DStack(particles, outFn, scaleFactor=scaleFactor)

        self.assertTrue(os.path.exists(outFn))
        self.assertTrue(pwutils.getFileSize(outFn) > 0)
        self.assertEqual(ih.getDimensions(outFn), EXPECTED_SIZE)
        self.assertEqual(ih.getDataType(outFn), DT)


class TestSetOfMicrographs(BaseTest):
    _labels = [SMALL, WEEKLY]

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.dbGold = cls.dataset.getFile('micsGoldSqlite')
        cls.micsPattern = cls.dataset.getFile('allMics')
        cls.dbFn = cls.getOutputPath('micrographs.sqlite')
        cls.acquisition = emobj.Acquisition(magnification=60000, voltage=300,
                                            sphericalAberration=2,
                                            amplitudeContrast=0.1)

        # cls.mics = glob(cls.micsPattern)
        cls.mics = []
        for mic in iglob(cls.micsPattern):
            cls.mics.append(cls.getRelPath(cls.dataset.getPath(), mic))

        if len(cls.mics) == 0:
            raise Exception('There are not micrographs matching pattern')
        cls.mics.sort()

    def checkSet(self, micSet):
        idCount = 1

        for mic1, fn in zip(micSet, self.mics):
            micFn = mic1.getFileName()
            self.assertEqual(fn, micFn,
                             "micrograph NAME in the set is wrong, \n   expected: '%s'\n        got: '%s'"
                             % (fn, micFn))
            self.assertEqual(idCount, mic1.getObjId(),
                             "micrograph ID in the set is wrong, \n   expected: '%s'\n        got: '%s'"
                             % (idCount, mic1.getObjId()))
            mic2 = micSet[idCount]  # Test getitem
            self.assertEqual(mic1.getObjId(), mic2.getObjId(), "micrograph got from ID is wrong")
            idCount += 1

    def testCreate(self):
        cwd = os.getcwd()
        # Change to test path
        os.chdir(self.dataset.getPath())
        """ Create a SetOfMicrographs from a list of micrographs """
        micSet = emobj.SetOfMicrographs(filename=self.dbFn)

        micSet.setAcquisition(self.acquisition)
        micSet.setSamplingRate(1.2)
        for fn in self.mics:
            mic = emobj.Micrograph()
            mic.setFileName(fn)
            micSet.append(mic)
        micSet.write()

        self.checkSet(micSet)

        # Check copy info also copies the samplingRate
        micSet2 = emobj.SetOfMicrographs(filename=':memory:')
        micSet2.copyInfo(micSet)
        self.assertAlmostEqual(micSet.getSamplingRate(), micSet2.getSamplingRate())

        os.chdir(cwd)

    def testRead(self):
        """ Read micrographs from a an sqlite file.
        It should contains Acquisition info. """
        micFn = self.dataset.getFile('micsGoldSqlite2')
        logger.info(">>> Reading gold micrographs from %s " % micFn)

        micSet = emobj.SetOfMicrographs(filename=micFn)
        self.assertEqual(2, micSet.getSize())
        acquisition = emobj.Acquisition()
        acquisition.setMagnification(10000.)
        acquisition.setVoltage(200.)
        acquisition.setSphericalAberration(2.26)
        acquisition.setAmplitudeContrast(0.1)

        mic2 = emobj.Micrograph()
        mic2.setSamplingRate(2.8)
        mic2.setAcquisition(acquisition)

        fileNames = ['/home/roberto/Scipion/Test/Test2/down2_12585',
                     '/home/roberto/Scipion/Test/Test2/down2_12610']
        counter = 0
        for mic in micSet:
            mic2.setFileName(fileNames[counter])
            self.assertTrue(mic.equalAttributes(mic2))
            counter += 1

    def test_mapper(self):
        """ test that indexes are created when a
        setOfParticles is created """
        MICNUMBER = 10
        indexesNames = ['_index']
        prot = createDummyProtocol("dummy_protocol")

        # create set of micrographs
        micSet = prot._createSetOfMicrographs()
        mic = emobj.Micrograph()
        for i in range(MICNUMBER + 1):
            mic.setLocation(i, "stack.stk")
            micSet.append(mic)
            mic.cleanObjId()
            micSet.write()

        # check defined indexes
        setOfMicrographsFileName = prot._getPath("micrographs.sqlite")
        indexes = sorted([index[1] for index in
                          getIndex(setOfMicrographsFileName)])
        for index, indexName in zip(indexes, indexesNames):
            self.assertEqual(index, 'index_' + indexName)


class TestSetOfParticles(BaseTest):
    """ Check if the information of the images is copied to another image when
    a new SetOfParticles is created"""

    _labels = [SMALL, WEEKLY]

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')

    def test_orderBy(self):
        # create setofProjections sorted
        imgSet1 = emobj.SetOfParticles(filename=':memory:', prefix='set1')
        imgSet2 = emobj.SetOfParticles(filename=':memory:', prefix='set2')
        imgSet1.setSamplingRate(1.5)
        imgSet2.setSamplingRate(1.5)
        img = emobj.Particle()

        for i in range(1, 10):
            img.setLocation(i, 'mystack.stk')
            img.setMicId(10 - i)
            img.setClassId(i % 5)
            imgSet1.append(img)
            img.setMicId(i)
            imgSet2.append(img)
            img.cleanObjId()
        # orderby
        for item1, item2 in zip(imgSet1.iterItems(orderBy='_micId',
                                                  direction='ASC'),
                                imgSet2.iterItems(orderBy='_micId',
                                                  direction='ASC')):
            self.assertEquals(item1.getMicId(), item2.getMicId())

    def test_readStack(self):
        """ Read an stack of 29 particles from .hdf file.
        Particles should be of 500x500 pixels.
        """
        size = 29
        xdim = 500
        inStack = self.dataset.getFile('particles1')
        outFn = self.getOutputPath('particles.sqlite')

        imgSet = emobj.SetOfParticles(filename=outFn)
        imgSet.setSamplingRate(1.0)
        imgSet.readStack(inStack)  # This should add 29 new items to the set

        self.assertEquals(size, imgSet.getSize())  # Check same size
        self.assertEquals(xdim, imgSet.getDim()[0])  # Check same dimensions

        logger.info("writing particles to %s" % outFn)
        imgSet.write()

        imgSet2 = emobj.SetOfParticles(filename=':memory:')
        imgSet2.copyInfo(imgSet)
        self.assertAlmostEqual(imgSet.getSamplingRate(), imgSet2.getSamplingRate())

    def test_hugeSet(self):
        """ Create a set of a big number of particles to measure
        creation time with sqlite operations. 
        """
        # Allow what huge means to be defined with environment var
        n = int(os.environ.get('SCIPION_TEST_HUGE', 10000))
        logger.info(">>>> Creating a set of %d particles." % n)
        logger.info("     (set SCIPION_TEST_HUGE environment var to other value)")

        dbFn = self.getOutputPath('huge_set.sqlite')
        # dbFn = ':memory:'

        img = emobj.Particle()
        imgSet = emobj.SetOfParticles(filename=dbFn)
        imgSet.setSamplingRate(1.0)

        for i in range(1, n + 1):
            # Creating object inside the loop significantly
            # decrease performance
            # img = Particle()
            img.setLocation(i, "images.stk")

            imgSet.append(img)
            img.cleanObjId()
            # img.setObjId(None)

        imgSet.write()

    def test_hugeSetToMd(self):
        """ Just as a benchrmark comparing to test_hugeSet ."""
        # Allow what huge means to be defined with environment var
        n = int(os.environ.get('SCIPION_TEST_HUGE', 10000))
        logger.info(">>>> Creating a set of %d particles." % n)
        logger.info("     (set SCIPION_TEST_HUGE environment var to other value)")

        imgMd = md.MetaData()

        for i in range(1, n + 1):
            # Creating object inside the loop significantly
            # decrease performance
            # img = Particle()
            objId = imgMd.addObject()
            imgMd.setValue(md.MDL_IMAGE, '%06d@images.stk' % (i + 1), objId)

        mdFn = self.getOutputPath('huge_set.xmd')
        logger.info("Writing metadata to %s." % mdFn)
        imgMd.write(mdFn)

    def test_hugeSetToText(self):
        """ Just as a benchrmark comparing to test_hugeSet ."""
        # Allow what huge means to be defined with environment var
        n = int(os.environ.get('SCIPION_TEST_HUGE', 10000))
        logger.info(">>>> Creating a set of %d particles." % n)
        logger.info("     (set SCIPION_TEST_HUGE environment var to other value)")

        textFn = self.getOutputPath('huge_set.txt')
        logger.info("Writing to text file %s " % textFn)
        f = open(textFn, 'w')

        for i in range(1, n + 1):
            string = '%06d@images.stk' % i
            f.write(string)

        f.close()

    def test_getFiles(self):
        # create setofImages
        dbFn = self.getOutputPath('multistack_set.sqlite')
        # dbFn = ':memory:'
        n = 10
        m = 3
        img = emobj.Particle()
        imgSet = emobj.SetOfParticles(filename=dbFn)
        imgSet.setSamplingRate(1.0)
        goldStacks = set()

        for j in range(1, m + 1):
            stackName = 'stack%02d.stk' % j
            goldStacks.add(stackName)
            for i in range(1, n + 1):
                img.setLocation(i, stackName)
                imgSet.append(img)
                img.cleanObjId()

        self.assertEquals(goldStacks, imgSet.getFiles())

        imgSet.close()

        # It should automatically call load
        # before accessing items db        
        imgSet.getFirstItem()

    def test_mapper(self):
        """ test that indexes are created when a
        setOfParticles is created """
        PARTNUMBER = 10
        indexesNames = ['_classId', '_micId']
        prot = createDummyProtocol("dummy_protocol")

        # create set of particles
        partSet = prot._createSetOfParticles()
        part = emobj.Particle()
        for i in range(PARTNUMBER + 1):
            part.setLocation(i, "stack.vol")
            partSet.append(part)
            part.cleanObjId()
        partSet.write()

        # check defined indexes
        setOfPArticleFileName = prot._getPath("particles.sqlite")
        indexes = sorted([index[1] for index in
                          getIndex(setOfPArticleFileName)])
        for index, indexName in zip(indexes, indexesNames):
            self.assertEqual(index, 'index_' + indexName)


class TestSetOfCoordinates(BaseTest):
    # TODO: A proper test for setOfCoordinates is missing
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def test_mapper(self):
        """ test that indexes are created when a
        setOfCoordinates is created """
        PARTNUMBER = 10
        MICNUMBER = 600
        NUMBERCOORDINATES = PARTNUMBER * MICNUMBER
        indexesNames = ['_micId']
        prot = createDummyProtocol("dummy_protocol")

        # create set of micrographs
        micSet = emobj.SetOfMicrographs(filename=":memory:")
        mic = emobj.Micrograph()
        for i in range(NUMBERCOORDINATES):
            mic.setLocation(i, "mic_%06d.mrc" % i)
            micSet.append(mic)
            mic.cleanObjId()
        micSet.write()

        # create a set of particles
        coordSet = prot._createSetOfCoordinates(micSet)
        coord = emobj.Coordinate()
        for i in range(NUMBERCOORDINATES):
            coordSet.append(coord)
            coord.cleanObjId()
        coordSet.write()

        # check defined indexes
        setOfCoordinatesFileName = \
            prot._getPath("coordinates.sqlite")
        logger.info(os.path.abspath(setOfCoordinatesFileName))
        indexes = sorted([index[1] for index in
                          getIndex(setOfCoordinatesFileName)])
        for index, indexName in zip(indexes, indexesNames):
            self.assertEqual(index, 'index_' + indexName)

        # Test speed: based on loop in file protocol_extractparticles.py
        # for 600 mic and 100 part the values for the first
        # second and third  case where:
        # Loop with index: 5 sec
        # Loop no index: 8:01 min
        # Loop optimized code: 4 sec
        # for 6000 mic and 200 part the values for the first
        # Loop with index: 1:47 min
        # optimized Loop with index: 1:20 min
        # Loop no index: after several  hours I stopped the process

        SPEEDTEST = True
        if SPEEDTEST:  # code from protocol_particles. line 415
            testTimer = pwutils.Timer()
            testTimer.tic()
            for mic in micSet:
                micId = mic.getObjId()
                coordList = []
                for coord in coordSet.iterItems(where='_micId=%s' % micId):
                    coordList.append(coord.clone())
            testTimer.toc("Loop with INDEX took:")

            lastMicId = None
            testTimer.tic()
            for coord in coordSet.iterItems(orderBy='_micId',
                                            direction='ASC'):
                micId = coord.getMicId()
                if micId != lastMicId:
                    lastMicId = micId
                    coordList = []
                coordList.append(coord.clone())
            testTimer.toc("Loop with INDEX and proper code, took:")

            # delete INDEX, this will not work
            # if database is not sqlite
            conn = sqlite3.connect(setOfCoordinatesFileName)
            cur = conn.cursor()
            for index in indexesNames:
                cur.execute("DROP INDEX index_%s" % index)
            cur.close()
            conn.close()

            testTimer.tic()
            for mic in micSet:
                micId = mic.getObjId()
                coordList = []
                for coord in coordSet.iterItems(where='_micId=%s' % micId):
                    coordList.append(coord.clone())
            testTimer.toc("Loop with NO INDEX took:")

            lastMicId = None
            testTimer.tic()
            for coord in coordSet.iterItems(orderBy='_micId',
                                            direction='ASC'):
                micId = coord.getMicId()
                if micId != lastMicId:
                    lastMicId = micId
                    coordList = []
                coordList.append(coord.clone())
            testTimer.toc("Loop with NO INDEX but proper code, took:")


class TestSetOfClasses2D(BaseTest):
    _labels = [SMALL, WEEKLY]

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('model')
        cls.selectionFn = cls.dataset.getFile('classesSelection')

    def test_basics(self):
        """ Load an existing SetOfClasses and test basic properties 
        such us: _mapperPath, iteration and others.
        """
        classes2DSet = emobj.SetOfClasses2D(filename=self.selectionFn)
        # Check the _mapperPath was properly set
        self.assertEqual(classes2DSet._mapperPath.get(), '%s, ' % self.selectionFn)

        cls1 = classes2DSet.getFirstItem()
        self.assertEqual(cls1._mapperPath.get(), '%s,Class001' % self.selectionFn)

        img1 = cls1.getFirstItem()
        self.assertEqual(img1.getObjId(), 1)  # First image of first class is 1 in this test

        images = [[1, 3, 4, 5, 7, 9, 11, 14, 16, 17, 21, 22, 24, 26, 28, 29, 30,
                   35, 37, 38, 39, 40, 42, 45, 54, 57, 62, 65, 67, 69, 70, 71, 74],
                  [2, 8, 10, 12, 13, 15, 19, 20, 23, 25, 27, 33, 34, 41, 43, 44,
                   46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 58, 59, 60, 61, 63,
                   64, 68, 72, 73, 75, 76],
                  [18, 31, 32],
                  [6, 36, 66]]

        classes2DSet.close()
        # Check iteration is working properly even after 
        # a close operation, it should open automatically
        for i, cls in enumerate(classes2DSet):
            l = images[i]
            for j, img in enumerate(cls):
                self.assertEquals(img.getObjId(), l[j])

        for i, rep in enumerate(classes2DSet.iterRepresentatives()):
            self.assertIsNotNone(rep.getLocation())

        # Check the SetOfClasses.iterClassItems method
        allImages = [img for imgList in images for img in imgList]
        idsImages = [img.getObjId()
                     for img in classes2DSet.iterClassItems(iterDisabled=True)]
        self.assertEqual(allImages, idsImages)

    def test_subsetsFromSelection(self):
        """ From a sqlite file of a SetOfClasses2D, with some
        classes and element disabled, we want to create a 
        subset of images and classes.
        """
        classes2DSet = emobj.SetOfClasses2D(filename=self.selectionFn)

        imgSet = emobj.SetOfParticles(filename=':memory:')
        # We are going to iterate over the enabled items and create
        # a new set of images
        imgSet.appendFromClasses(classes2DSet)
        # Since we have disabled two classes (6 images) and 1 images
        # from class 1 and 2 (2 images), the final size of the set
        # should be 68
        sizes = [32, 36]
        self.assertEqual(imgSet.getSize(), sum(sizes))
        imgSet.clear()  # Close db connection and clean data

        # Now create a subset of classes and check the number
        # of images per class
        clsSet = emobj.SetOfClasses2D(filename=':memory:')
        clsSet.appendFromClasses(classes2DSet)
        for i, cls in enumerate(clsSet):
            self.assertEqual(cls.getSize(), sizes[i])
        clsSet.clear()  # Close db connection and clean data


class TestTransform(BaseTest):

    def test_scale(self):
        """ Check Scale storage in transformation class
        """
        t = emobj.Transform()
        m = t.getMatrix()
        m[0, 3] = 2
        m[1, 3] = 4
        m[2, 3] = 6
        m[3, 3] = 5
        t.scale(0.5)

        self.assertAlmostEqual(m[0, 3], 1)
        self.assertAlmostEqual(m[1, 3], 2)
        self.assertAlmostEqual(m[2, 3], 3)
        self.assertAlmostEqual(m[3, 3], 1)

    def test_scaleShifts(self):
        """ Check Scale 2D shifts in transformation class
        """
        t = emobj.Transform()
        m = t.getMatrix()
        m[0, 3] = 2
        m[1, 3] = 4
        m[2, 3] = 6
        m[3, 3] = 5
        t.scaleShifts(0.5)

        self.assertAlmostEqual(m[0, 3], 1)
        self.assertAlmostEqual(m[1, 3], 2)
        self.assertAlmostEqual(m[2, 3], 3)
        self.assertAlmostEqual(m[3, 3], 5)

    def test_clone(self):
        """ Check that cloning the Transform will 
        also copy the values of underlying numpy matrix.
        """
        t = emobj.Transform()
        m = t.getMatrix()
        m[0, 3] = 2
        m[1, 3] = 4

        t2 = t.clone()
        m2 = t2.getMatrix()
        self.assertTrue(np.allclose(m, m2, rtol=1e-2))

        p = emobj.Particle()
        p.setTransform(t)

        p2 = p.clone()
        m3 = p2.getTransform().getMatrix()
        self.assertTrue(np.allclose(m, m3, rtol=1e-2))


class TestCopyItems(BaseTest):
    _labels = [SMALL, WEEKLY]

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('model')

    def test_listAttribute(self):
        """ 
        Test the use of copyItems function where the items have
        a property that is a list.
        """

        inFn = self.dataset.getFile('particles/particles_nma.sqlite')
        outFn = self.getOutputPath('particles.sqlite')

        inputSet = emobj.SetOfParticles(filename=inFn)
        inputSet.setSamplingRate(1.0)

        outputSet = emobj.SetOfParticles(filename=outFn)
        outputSet.copyInfo(inputSet)

        outputSet.copyItems(inputSet,
                            updateItemCallback=self._updateItem)

        outputSet.write()
        outputSet.close()

        # Read again the written set of particles
        checkSet = emobj.SetOfParticles(filename=outFn)

        # check that the first item (as the rest)
        # have the added _list attribute with the correct values
        particle = checkSet.getFirstItem()
        self.assertTrue(particle.hasAttribute('_list'))
        for v1, v2 in zip([1.0, 2.0], particle._list):
            self.assertAlmostEqual(v1, float(v2))

        # check that copied items have exactly
        # the same attributes than the input set
        # except for the _list
        for i1, i2 in zip(inputSet, checkSet):
            self.assertTrue(i2.equalAttributes(i1, ignore=['_list']))

    def _updateItem(self, item, row):
        item._list = emobj.CsvList()
        item._list.set([1.0, 2.0])


class TestCoordinatesTiltPair(BaseTest):
    # TODO: A proper test for CoordinatesTiltPair is missing
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def test_mapper(self):
        """ test that indexes are created when a
        setOfCoordinates is created """
        MICNUMBER = 10
        indexesNames = ['_untilted._micId',
                        '_tilted._micId']
        prot = createDummyProtocol("dummy_protocol")

        # create set of untilted and tilted micrographs
        uMicSet = emobj.SetOfMicrographs(filename=":memory:")
        tMicSet = emobj.SetOfMicrographs(filename=":memory:")
        mic = emobj.Micrograph()
        for i in range(MICNUMBER):
            mic.setLocation(i, "umic_%06d.mrc" % i)
            uMicSet.append(mic)
            mic.cleanObjId()
            mic.setLocation(i, "tmic_%06d.mrc" % i)
            tMicSet.append(mic)
            mic.cleanObjId()
        # TODO I do not see any example on how to create
        # a set of CoordinatesTiltPair. I can image
        # that a need to create a set of micTiltPairs
        # and two sets of coordinates but the person who
        # added that data type Should provide a clear test
        # when this is done then I will finish the test_mapper


class TestCTFModel(TestCase):
    def test_stringContext(self):
        ctf = emobj.CTFModel()
        # All values to None should be printable without exception
        logger.info("Empty Ctf str converts to %s" % ctf)

        ctf.setDefocusU(2.4545)
        logger.info("Ctf with defocusU converts to %s" % ctf)
        self.assertTrue("2.45," in str(ctf))


class TestSequenceHandler(pwtests.BaseTest):
    NAME = 'USER_SEQ'
    DESCRIPTION = 'User description'
    UNIPROTID1 = 'P12345'
    UNIPROTID2 = 'P03252'

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)

    def testImportStructureAminoacidSequence1(self):
        """
        Import the sequence of chain B of atomic structure 3lqd.cif
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    emprot.ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID1,
                }
        prot1 = self.newProtocol(emprot.ProtImportSequence, **args)
        self.launchProtocol(prot1)
        sequence1 = prot1.outputSequence

        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    emprot.ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID2,
                }
        prot2 = self.newProtocol(emprot.ProtImportSequence, **args)
        self.launchProtocol(prot2)
        sequence2 = prot2.outputSequence

        #Saving Sequences into files
        outFile = prot1._getPath('sequence.fasta')
        sequence1.exportToFile(outFile)
        # sequence2.exportToFile(prot2._getPath('sequence.genbank'))

        #Reading exported sequence
        readSequence = SequenceHandler()
        sequencesDics = readSequence.readSequencesFromFile(outFile)
        self.assertEqual(sequence1.getSequence(), sequencesDics[0]['sequence'])

        #Appending sequence to existing sequence file
        sequence2.appendToFile(outFile)

        # Reading appended sequences
        readSequence = SequenceHandler()
        sequencesDics = readSequence.readSequencesFromFile(outFile)
        self.assertEqual(sequence1.getSequence(), sequencesDics[0]['sequence'])
        self.assertEqual(sequence2.getSequence(), sequencesDics[1]['sequence'])


