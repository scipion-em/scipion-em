# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from os.path import exists, basename, abspath, relpath, join
from os import stat
from numpy import array
from numpy.linalg import norm

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pwem.objects as emobj
import pwem.convert as emconv
from pwem import emlib

from .base import ProtImportFiles
from .images import ProtImportImages
from ...convert import Ccp4Header


class ProtImportVolumes(ProtImportImages):
    """Protocol to import a set of volumes to the project"""
    _outputClassName = 'SetOfVolumes'
    _label = 'import volumes'
    IMPORT_FROM_EMDB = 1

    def __init__(self, **args):
        ProtImportImages.__init__(self, **args)

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        return ['local file', 'EMDBid']

    def _defineAcquisitionParams(self, form):
        """ Define acquisition parameters, it can be overriden
        by subclasses to change what parameters to include.
        """
        param = form.getParam('importFrom')
        # Customize the help of this parameter with specific information
        # of the import particles
        param.help.set('Import 3D map data from EMDB server or local file')
        form.addParam('emdbId', params.IntParam,
                      condition='importFrom != IMPORT_FROM_FILES',
                      label="EMDB map ID (integer)", allowsNull=True,
                      help='Type an EMDB ID (Integer). ')
        form.addParam('setHalfMaps', params.BooleanParam,
                      condition='importFrom == IMPORT_FROM_FILES',
                      label='Set half maps',
                      help='Option YES:\nAssign two half maps to the imported map.',
                      default=False)
        form.addParam('half1map', params.PathParam,
                      label='Path half map1', help='Select first half map',
                      condition='setHalfMaps')
        form.addParam('half2map', params.PathParam,
                      label='Path half map2', help='Select second half map',
                      condition='setHalfMaps')
        form.addParam('samplingRate', params.FloatParam,
                      condition='importFrom == IMPORT_FROM_FILES',
                      label=pwutils.Message.LABEL_SAMP_RATE)
        form.addParam('setOrigCoord', params.BooleanParam,
                      condition='importFrom == IMPORT_FROM_FILES',
                      label="Set origin of coordinates",
                      help="Option YES:\nA new volume will be created with "
                           "the "
                           "given ORIGIN of coordinates. This ORIGIN will be "
                           "set in the map file header.\nThe ORIGIN of "
                           "coordinates will be placed at the center of the "
                           "whole volume if you select n(x)/2, n(y)/2, "
                           "n(z)/2 as "
                           "x, y, z coordinates (n(x), n(y), n(z) are the "
                           "dimensions of the whole volume). However, "
                           "selecting "
                           "0, 0, 0 as x, y, z coordinates, the volume will be "
                           "placed at the upper right-hand corner.\n\n"
                           "Option NO:\nThe ORIGIN of coordinates will be "
                           "placed at the center of the whole volume ("
                           "coordinates n(x)/2, n(y)/2, n(z)/2 by default). "
                           "This "
                           "ORIGIN will NOT be set in the map file header.\n\n"
                           "WARNING: In case you want to process "
                           "the volume with programs requiring a specific "
                           "symmetry regarding the origin of coordinates, "
                           "for example the protocol extract unit "
                           "cell, check carefully that the coordinates of the "
                           "origin preserve the symmetry of the whole volume. "
                           "This is particularly relevant for loading "
                           "fragments/subunits of the whole volume.\n",
                      default=False)
        line = form.addLine('Offset',
                            help="A wizard will suggest you possible "
                                 "coordinates for the ORIGIN. In MRC volume "
                                 "files, the ORIGIN coordinates will be "
                                 "obtained from the file header.\n "
                                 "In case you prefer set your own ORIGIN "
                                 "coordinates, write them here. You have to "
                                 "provide the map center coordinates in "
                                 "Angstroms (pixels x sampling).\n",
                            condition='setOrigCoord')
        # line.addParam would produce a nicer looking form
        # but them the wizard icon is drawn outside the visible
        # window. Until this bug is fixed form is a better option
        form.addParam('x', params.FloatParam, condition='setOrigCoord',
                      label="x", help="offset along x axis (Angstroms)")
        form.addParam('y', params.FloatParam, condition='setOrigCoord',
                      label="y", help="offset along y axis (Angstroms)")
        form.addParam('z', params.FloatParam, condition='setOrigCoord',
                      label="z", help="offset along z axis (Angstroms)")

    def _insertAllSteps(self):
        self._insertFunctionStep('importVolumesStep',
                                 self.getPattern(),
                                 self.samplingRate.get(),
                                 self.setOrigCoord.get())

    # --------------------------- STEPS functions -----------------------------

    def importVolumesStep(self, pattern, samplingRate, setOrigCoord=False):
        """ Copy images matching the filename pattern
        Register other parameters.
        """
        volSet = self._createSetOfVolumes()
        vol = emobj.Volume()

        if self.importFrom == self.IMPORT_FROM_FILES:
            self.info("Using pattern: '%s'" % pattern)

            # Create a Volume template object
            vol.setSamplingRate(samplingRate)

            imgh = emlib.image.ImageHandler()

            volSet.setSamplingRate(samplingRate)

            for fileName, fileId in self.iterFiles():
                x, y, z, n = imgh.getDimensions(fileName)
                if fileName.endswith('.mrc') or fileName.endswith('.map'):
                    fileName += ':mrc'
                    if z == 1 and n != 1:
                        zDim = n
                        n = 1
                    else:
                        zDim = z
                else:
                    zDim = z
                origin = emobj.Transform()
                if setOrigCoord:
                    origin.setShiftsTuple(self._getOrigCoord())
                else:
                    origin.setShifts(x / -2. * samplingRate,
                                     y / -2. * samplingRate,
                                     zDim / -2. * samplingRate)

                vol.setOrigin(origin)  # read origin from form

                if self.copyFiles or setOrigCoord:
                    newFileName = abspath(self._getVolumeFileName(fileName, "mrc"))
                    emconv.Ccp4Header.fixFile(fileName, newFileName, origin.getShifts(),
                                              samplingRate, emconv.Ccp4Header.ORIGIN)
                    if self.setHalfMaps.get():
                        newFileName1 = abspath(self._getVolumeFileName(self.half1map.get(), "mrc"))
                        emconv.Ccp4Header.fixFile(self.half1map.get(), newFileName1, origin.getShifts(),
                                                  samplingRate, emconv.Ccp4Header.ORIGIN)
                        newFileName2 = abspath(self._getVolumeFileName(self.half2map.get(), "mrc"))
                        emconv.Ccp4Header.fixFile(self.half2map.get(), newFileName2, origin.getShifts(),
                                                  samplingRate, emconv.Ccp4Header.ORIGIN)

                        vol.setHalfMaps([relpath(newFileName1), relpath(newFileName2)])
                else:
                    newFileName = abspath(self._getVolumeFileName(fileName))

                    if fileName.endswith(':mrc'):
                        fileName = fileName[:-4]

                    pwutils.createAbsLink(fileName, newFileName)
                    if self.setHalfMaps.get():
                        pwutils.createAbsLink(self.half1map.get(),
                                              abspath(self._getVolumeFileName(self.half1map.get())))
                        pwutils.createAbsLink(self.half2map.get(),
                                              abspath(self._getVolumeFileName(self.half2map.get())))

                        vol.setHalfMaps([relpath(self._getVolumeFileName(self.half1map.get())),
                                         relpath(self._getVolumeFileName(self.half2map.get()))
                                         ])

                # Make newFileName relative
                # https://github.com/I2PC/scipion/issues/1935
                newFileName = relpath(newFileName)
                if n == 1:
                    vol.cleanObjId()
                    vol.setFileName(newFileName)
                    volSet.append(vol)
                else:
                    for index in range(1, n + 1):
                        vol.cleanObjId()
                        vol.setLocation(index, newFileName)
                        volSet.append(vol)
        else:  # import from EMDB
            self.info("Downloading map with ID = %s" % self.emdbId)
            try:
                localFileName, sampling, origin = \
                    fetch_emdb_map(self.emdbId,
                                   self._getExtraPath(),
                                   self._getTmpPath())
            except Exception as e:
                print(e)
                return
            # open volume and fill sampling and origin
            vol.setSamplingRate(sampling)
            vol.setFileName(localFileName)
            from pwem.objects.data import Transform
            originMat = Transform()
            originMat.setShifts(origin[0], origin[1], origin[2])
            vol.setOrigin(originMat)

        if volSet.getSize() > 1:
            self._defineOutputs(outputVolumes=volSet)
        else:
            self._defineOutputs(outputVolume=vol)

    # --------------------------- INFO functions ------------------------------

    def _getVolMessage(self):
        if self.hasAttribute('outputVolume'):
            return "Volume %s" % self.getObjectTag('outputVolume')
        else:
            return "Volumes %s" % self.getObjectTag('outputVolumes')

    def _summary(self):
        summary = []
        if self.hasAttribute('outputVolume') or \
                self.hasAttribute('outputVolumes'):
            summary.append("%s imported from:\n%s" % (self._getVolMessage(),
                                                      self.getPattern()))

            summary.append(u"Sampling rate: *%0.2f* (Å/px)" %
                           self.samplingRate.get())
        return summary

    def _methods(self):
        methods = []
        if self.hasAttribute('outputVolume') or \
                self.hasAttribute('outputVolumes'):
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getVolMessage(), self.samplingRate.get()), )
        return methods

    def _getVolumeFileName(self, fileName, extension=None):
        if extension is not None:
            baseFileName = "import_" + basename(fileName).split(".")[0] + ".%s" % extension
        else:
            baseFileName = "import_" + basename(fileName).split(":")[0]

        return self._getExtraPath(baseFileName)

    def _getOrigCoord(self):
        return -1. * self.x.get(), -1. * self.y.get(), -1. * self.z.get()

    def _validate(self):
        errors = super(ProtImportVolumes, self)._validate()
        if (not self.filesPattern.empty()) and self.setHalfMaps.get():
            errors.append("You can not use the options 'Pattern' "
                          "and 'Set half maps' simultaneously")

        return errors


class ProtImportPdb(ProtImportFiles):
    """ Protocol to import an atomic structure  to the project.
Format may be PDB or MMCIF"""
    _label = 'import atomic structure'
    IMPORT_FROM_ID = 0
    IMPORT_FROM_FILES = 1
    SYMPLIFIED_STRUCT = "symplified_atom_structure.cif"

    def __init__(self, **args):
        ProtImportFiles.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPdbData', params.EnumParam, choices=['id', 'file'],
                      label="Import atomic structure from",
                      default=self.IMPORT_FROM_ID,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Import mmCIF data from online server or local file')
        form.addParam('pdbId', params.StringParam,
                      condition='inputPdbData == IMPORT_FROM_ID',
                      label="Atomic structure ID ", allowsNull=True,
                      help='Type a mmCIF ID (four alphanumeric characters).')
        form.addParam('pdbFile', params.PathParam, label="File path",
                      condition='inputPdbData == IMPORT_FROM_FILES',
                      allowsNull=True,
                      help='Specify a path to desired atomic structure.')
        form.addParam('inputVolume', params.PointerParam, label="Input Volume",
                      pointerClass='Volume',
                      allowsNull=True,
                      help='Associate this volume to the mmCIF file.')

    def _insertAllSteps(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            self._insertFunctionStep('pdbDownloadStep')
        else:
            self._insertFunctionStep('createOutputStep', self.pdbFile.get())

    def pdbDownloadStep(self):
        """Download all pdb files in file_list and unzip them.
        """
        aSH = emconv.AtomicStructHandler()
        print("retrieving atomic structure with ID = %s" % self.pdbId.get())
        atomStructPath = aSH.readFromPDBDatabase(self.pdbId.get(),
                                                 type='mmCif',
                                                 dir=self._getExtraPath())
        self.createOutputStep(atomStructPath)

    #        downloadPdb(self.pdbId.get(), pdbPath, self._log)

    def createOutputStep(self, atomStructPath):
        """ Copy the PDB structure and register the output object.
        """
        if not exists(atomStructPath):
            raise Exception("Atomic structure not found at *%s*" % atomStructPath)

        baseName = basename(atomStructPath)
        localPath = abspath(self._getExtraPath(baseName))

        if str(atomStructPath) != str(localPath):  # from local file
            pwutils.copyFile(atomStructPath, localPath)

        localPath = relpath(localPath)

        pdb = emobj.AtomStruct()
        volume = self.inputVolume.get()

        # if a volume exists assign it to the pdb object
        # IMPORTANT: we DO need "if volume is not None"
        # because we need to persist the pdb object
        # before we can make the last source relation
        if volume is not None:
            pdb.setVolume(volume)

        pdb.setFileName(localPath)
        self._defineOutputs(outputPdb=pdb)

        if volume is not None:
            self._defineSourceRelation(volume, pdb)

    def _summary(self):
        if self.inputPdbData == self.IMPORT_FROM_ID:
            summary = ['Atomic structure imported from ID: *%s*' %
                       self.pdbId]
        else:
            summary = ['Atomic structure imported from file: *%s*' %
                       self.pdbFile]

        return summary

    def _validate(self):
        errors = []
        if (self.inputPdbData == self.IMPORT_FROM_FILES and not exists(
                self.pdbFile.get())):
            errors.append("Atomic structure not found at *%s*" %
                          self.pdbFile.get())
        # TODO: maybe also validate that if exists is a valid PDB file
        return errors


######################################


def fetch_emdb_map(id, directory, tmpDirectory):
    """ get map from emd
    :param id: 3D MAP ID in EMDB
    :return: local 3Dmap filename
    """
    import socket

    # get computer name and select server
    url_rest_api = "https://www.ebi.ac.uk/pdbe/api/emdb/entry/map/EMD-%d"
    hname = socket.gethostname()
    if hname.endswith('.edu') or hname.endswith('.gov'):
        site = 'ftp.wwpdb.org'
        url_pattern = 'ftp://%s/pub/emdb/structures/EMD-%s/map/%s'
    elif hname.endswith('.cn'):
        site = 'ftp.emdb-china.org'
        url_pattern = 'ftp://%s/structures/EMD-%s/map/%s'
    else:
        site = 'ftp.ebi.ac.uk'
        url_pattern = 'ftp://%s/pub/databases/emdb/structures/EMD-%s/map/%s'

    map_name = 'emd_%s.map' % id
    map_gz_name = map_name + '.gz'
    map_url = url_pattern % (site, id, map_gz_name)
    name = 'EMD-%d' % id
    minimum_map_size = 8192  # bytes
    url_rest_api = url_rest_api % id

    try:
        map_path, samplingAPI, originAPI = fetch_file(map_url,
                                                url_rest_api,
                                                name,
                                                minimum_map_size,
                                                directory,
                                                tmpDirectory,
                                                map_name
                                                )
    except Exception as e:
        raise Exception ("Cannot retrieve File from EMDB", e)

    originAPI = array(originAPI) * samplingAPI  # convert to Angstrom
    #check consistency between file header and rest API
    ccp4header = Ccp4Header(map_path, readHeader=True)
    samplingHeader = ccp4header.computeSampling()  # unit = A/px
    originHeader = array(ccp4header.getOrigin())   # unit = A

    if abs(samplingHeader - samplingAPI) >= 0.01:
        print("###########################\n"
              "WARNING: sampling rate stored in EMDB\n"
              "database and 3D map header file do not match\n"
              "API=%f, header=%f\n"
              "###########################\n" % (samplingAPI, samplingHeader))

    if norm(originHeader - originAPI) >= 0.1:
        print("###########################\n"
              "WARNING: origin  stored in EMDB\n"
              "database and 3D map header file do not match\n"
              "API=%f, header=%f\n"
              "###########################\n" % (originAPI, originHeader))
    return map_path, samplingAPI, originAPI


def fetch_file(url, url_rest_api, name,
               minimum_file_size=8192,
               save_dir='',
               tmp_dir='',
               save_name=''):
    """
    Download 3DMAPfile from EMDB

    :param name: EMB name id, format -> EMD-id
    :param url_rest_api: EMDB rest api, here we ay ask for 3D map properties
    :param tmp_dir: directory in which the compressed
                    file will be saved temporarily
    :param url:  3D map url
    :param minimum_file_size:
    :param save_dir: save file in this directory
    :param save_name:  local filename
    :return: local file name
    """
    import urllib.request
    import requests
    noCompressName = join(save_dir, save_name)
    compressName = join(tmp_dir, save_name + ".gz")

    try:
        urllib.request.urlretrieve(url, filename=compressName)
        # if retrieval fails retry another time
        if not exists(compressName):
            urllib.request.urlretrieve(url, filename=compressName)
            if not exists(compressName):
                raise Exception("Can not download file from EMDB")

        json_results = requests.get(url_rest_api).json()
        sampling_tag = "pixel_spacing"
        sampling_x_tag = "x"
        origin_tag = "origin"
        origin_x_tag = "column"
        origin_y_tag = "row"
        origin_z_tag = "section"

        # units A/px
        results = json_results[name][0]['map']
        sampling = results[sampling_tag][sampling_x_tag]['value']
        # units unknown may be pixels since this is integer
        x = results[origin_tag][origin_x_tag]
        y = results[origin_tag][origin_y_tag]
        z = results[origin_tag][origin_z_tag]
    except Exception as e:
        print("Error retriving data from EMDB", str(e))

    if stat(compressName).st_size < minimum_file_size:
        raise Exception("File Downloaded from EMDB is empty")

    gunzip(compressName, noCompressName)
    return noCompressName, sampling, (x, y, z)


def gunzip(gzpath, path):
    import gzip
    gzf = gzip.open(gzpath)
    f = open(path, 'wb')
    f.write(gzf.read())
    f.close()
    gzf.close()
