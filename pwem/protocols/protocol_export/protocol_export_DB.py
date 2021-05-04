# **************************************************************************
# *
# * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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


import os

import pyworkflow.protocol.params as params
from pwem.convert import Ccp4Header
from pyworkflow import VERSION_1_2

from pwem.convert.atom_struct import fromPDBToCIF, fromCIFTommCIF, \
    AtomicStructHandler
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pwem.objects import FSC
from pyworkflow.utils.path import copyFile


class ProtExportDataBases(EMProtocol):
    """ generates files for elements to submit structures to EMDB/PDB.
        Since mmcif/pdb is only partially supported by some software
        the protocol creates 4 versions of the atomic struct file with the hope that at least
        one of them will work.
    """
    _label = 'export emdb/pdb'
    _program = ""
    _lastUpdateVersion = VERSION_1_2
    VOLUMENAME = 'main_map.mrc'
    HALFVOLUMENAME = 'half_map_%d.mrc'
    COORDINATEFILENAME = 'coordinates.cif'
    ADDITIONALVOLUMEDIR = "addMaps"
    ADDITIONALVOLUMENAME = "map_%02d.mrc"
    MASKDIR = "masks"
    MASKNAME = "mask_%02d.mrc"
    SYMPLIFIED_STRUCT = "symplified_atom_structure.cif"

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('exportVolume', params.PointerParam,
                      label="Main EM map to export",
                      allowsNull=True,
                      pointerClass='Volume',
                      help='This EM map is mandatory for EMDB and it '
                           'will be exported using mrc format. '
                           'If this map is associated to their respective '
                           'half maps, they will be exported as well.')
        form.addParam('additionalVolumesToExport', params.BooleanParam,
                      default=False, label='Additional maps to export?',
                      help='Select YES if you want to add some more '
                           'EM maps to export.')
        form.addParam('exportAdditionalVolumes', params.MultiPointerParam,
                      label="Additional EM maps to export",
                      allowsNull=True,
                      condition='additionalVolumesToExport == True',
                      pointerClass='Volume',
                      help='These additional EM maps will be also exported '
                           'using mrc format.')
        form.addParam('exportFSC', params.PointerParam, label="FSC file to export",
                      allowsNull=True,
                      pointerClass='FSC, SetOfFSCs',
                      help='This FSCs will be exported using XML format')
        form.addParam('masksToExport', params.BooleanParam,
                      default=False, label='Masks to export?',
                      help='Select YES if you want to add some  '
                           'masks to export.')
        form.addParam('exportMasks', params.MultiPointerParam, label="Masks to export",
                      allowsNull=True, condition='masksToExport == True',
                      pointerClass='Mask',
                      help='These mask will be exported using mrc format')
        form.addParam('exportAtomStruct', params.PointerParam,
                      label="Atomic structure to export", allowsNull=True,
                      pointerClass='AtomStruct',
                      help='This atomic structure will be exported using mmCIF format')
        form.addParam('exportPicture', params.PathParam,
                      label="Image to export", allowsNull=True,
                      pointerClass='Image',
                      help='This image is mandatory for EMDB')
        form.addParam('filesPath', params.PathParam, important=True,
                      label="Export to directory",
                      help="Directory where the files will be generated.")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.dirName = self.filesPath.get()
        self._insertFunctionStep('createDirectoryStep', self.dirName)
        if self.exportVolume.get() is not None:
            self._insertFunctionStep('exportVolumeStep')
        if self.additionalVolumesToExport:
            self._insertFunctionStep('exportAdditionalVolumeStep')
        if self.exportFSC.get() is not None:
            self._insertFunctionStep('exportFSCStep')
        if self.masksToExport:
            self._insertFunctionStep('exportMasksStep')
        if self.exportAtomStruct.get() is not None:
            self._insertFunctionStep('exportAtomStructStep')
        if self.exportPicture.get() is not None:
            self._insertFunctionStep('exportImageStep')

    # --------------------------- STEPS functions -----------------------------

    def createDirectoryStep(self, dirPath):
        try:
            os.makedirs(dirPath)
        except OSError:
            if not os.path.isdir(dirPath):
                print("Can not create directory %s" % dirPath)
                raise

    def exportVolumeStep(self):
        inVolFileName = self.exportVolume.get().getFileName()
        inVol = self.exportVolume.get()
        outVolFileName = os.path.join(self.dirName, self.VOLUMENAME)
        shifts = inVol.getOrigin(force=True).getShifts()
        sampling = inVol.getSamplingRate()

        ccp4header = Ccp4Header(inVolFileName)
        ccp4header.fixFile(inVolFileName, outVolFileName, shifts,
                           sampling=sampling)

        # Do we have half volumes?
        if self.exportVolume.get().hasHalfMaps():
            ih = ImageHandler()
            for counter, half_map in enumerate(
                    self.exportVolume.get().getHalfMaps().split(','), 1):
                outVolFileName = os.path.join(self.dirName,
                                              self.HALFVOLUMENAME % counter)
                ccp4header.fixFile(half_map, outVolFileName, shifts,
                                   sampling=sampling)

    def exportAdditionalVolumeStep(self):
        outputDir = os.path.join(self.dirName, self.ADDITIONALVOLUMEDIR)
        self.createDirectoryStep(outputDir)
        ih = ImageHandler()
        for counter, map in enumerate(self.exportAdditionalVolumes, 1):
            map = map.get()
            inVolFileName = map.getFileName()
            outVolFileName = os.path.join(outputDir,
                                          self.ADDITIONALVOLUMENAME % counter)
            shifts = map.getOrigin(force=True).getShifts()
            sampling = map.getSamplingRate()

            ccp4header = Ccp4Header(inVolFileName)
            ccp4header.fixFile(inVolFileName, outVolFileName, shifts,
                               sampling=sampling)

    def exportFSCStep(self):
        exportFSC = self.exportFSC.get()
        if isinstance(self.exportFSC.get(), FSC):
            fscSet = self._createSetOfFSCs()
            fscSet.append(exportFSC)
        else:
            fscSet = exportFSC

        dirName = self.filesPath.get()
        for i, exportFSC in enumerate(fscSet, 1):
            x, y = exportFSC.getData()
            fnFSC = os.path.join(dirName, "fsc_%02d.xml" % i)
            fo = open(fnFSC, "w")
            fo.write('<fsc title="FSC(%s)" xaxis="Resolution (A-1)" '
                     'yaxis="Correlation Coefficient">\n' %
                     os.path.join(dirName, self.VOLUMENAME))
            for k in range(len(x)):
                fo.write("<coordinate>\n")
                fo.write("<x>%f</x>\n" % x[k])
                fo.write("<y>%f</y>\n" % y[k])
                fo.write("</coordinate>\n")
            fo.write("</fsc>\n")
            fo.close()

    def exportMasksStep(self):
        outputDir = os.path.join(self.dirName, self.MASKDIR)
        self.createDirectoryStep(outputDir)

        for counter, mask in enumerate(self.exportMasks, 1):
            mask = mask.get()
            inVolFileName = mask.getFileName()
            outVolFileName = os.path.join(outputDir,
                                          self.MASKNAME % counter)
            shifts = mask.getOrigin(force=True).getShifts()
            sampling = mask.getSamplingRate()

            ccp4header = Ccp4Header(inVolFileName)
            ccp4header.fixFile(inVolFileName, outVolFileName, shifts,
                               sampling=sampling)

    def exportAtomStructStep(self):
        exportAtomStruct = self.exportAtomStruct.get()
        originStructPath = exportAtomStruct.getFileName()
        dirName = self.filesPath.get()
        destinyStructPath = os.path.join(dirName, self.COORDINATEFILENAME)
        destinySympleStructPath = os.path.join(dirName, self.SYMPLIFIED_STRUCT)

        # save input atom struct with no change
        baseName = os.path.basename(originStructPath)
        localPath = os.path.abspath(os.path.join(dirName, baseName))
        copyFile(originStructPath, localPath)

        # call biopython to simplify atom struct and save it
        aSH = AtomicStructHandler()
        aSH.read(originStructPath)
        aSH.write(destinySympleStructPath)

        # if pdb convert to mmcif calling maxit twice
        if originStructPath.endswith(".pdb"):
            # convert pdb to cif using maxit program
            log = self._log
            fromPDBToCIF(originStructPath,
                         destinyStructPath, log)
            try:
                # convert cif to mmCIF by using maxit program
                fromCIFTommCIF(destinyStructPath,
                               destinyStructPath, log)
            except Exception as e:
                pass
        # if cif convert to mmcif using maxit
        elif originStructPath.endswith(".cif"):
            # convert cif to mmCIF by using maxit program
            log = self._log
            try:
                fromCIFTommCIF(originStructPath,
                               destinyStructPath, log)
            except Exception as e:
                pass

    def exportImageStep(self):
        imageBaseFileName = os.path.basename(self.exportPicture.get())
        outputFile = os.path.join(self.dirName, imageBaseFileName)
        copyFile(self.exportPicture.get(), outputFile)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        message = []
        fnPath = self.filesPath.get()
        if fnPath == "" or fnPath is None:
            message.append("You must set a path to export.")
        return message

    def _summary(self):
        message = "Data Available at : *%s*" % self.filesPath.get()
        return [message]

    def _methods(self):
        return []

    # --------------------------- UTILS functions ---------------------------------

    def getFnPath(self, label='volume'):
        return os.path.join(self.filesPath.get(),
                            self._getFileName(label))
