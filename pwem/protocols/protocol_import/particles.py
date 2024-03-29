# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from os.path import abspath

import pyworkflow.utils as pwutils
import pyworkflow.protocol.constants as pwcts
import pyworkflow.protocol.params as params

from pwem import Domain
import pwem.constants as emcts

from .images import ProtImportImages


class ProtImportParticles(ProtImportImages):
    """Protocol to import a set of particles to the project"""
    _label = 'import particles'
    _outputClassName = 'SetOfParticles'

    IMPORT_FROM_EMX = 1
    IMPORT_FROM_XMIPP3 = 2
    IMPORT_FROM_RELION = 3
    IMPORT_FROM_SCIPION = 4
    IMPORT_FROM_FREALIGN = 5
    IMPORT_FROM_EMAN = 6
    IMPORT_FROM_CRYOSPARC = 7

    importFormats = ['emx', 'xmipp3', 'relion', 'scipion', 'frealign', 'eman', 'cryosparc']
    importExts = ['emx', 'xmd', 'star', 'sqlite', 'par', 'lst', 'cs']
    alignTypeList = [emcts.ALIGN_2D, emcts.ALIGN_3D, emcts.ALIGN_PROJ,
                     emcts.ALIGN_NONE]

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.)
        """
        choices = ProtImportImages._getImportChoices(self)
        # Do not change the order of this list since
        # it is related to the constants defined
        return choices + self.importFormats
    
    def _defineImportParams(self, form):
        """
        Import files from: emx, xmipp3, relion, scipion, cryosparc formats.
        """
        param = form.getParam('importFrom')
        # Customize the help of this parameter with specific information
        # of the import particles
        param.help.set("You can import particles directly from the binary "
                       "files, or import from other packages formats. \n"
                       "Currently, we can import from: %s \n"
                       "Following are the expected import files for each one:\n"
                       "*emx*: particles.emx\n"
                       "*xmipp3*: images.xmd\n"
                       "*relion*: itXX_data.star\n"
                       "*scipion*: particles.sqlite\n"
                       "*eman*: particleSet.lst\n"
                       "*cryosparc*: particles.cs\n"
                       "" % ', '.join(self.importFormats))

        form.addParam('emxFile', params.FileParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_EMX,
                      label='Input EMX file',
                      help="Select the EMX file containing particles "
                           "information.\n See more about \n"
                           "[[http://i2pc.cnb.csic.es/emx][EMX format]]")

        form.addParam('alignType', params.EnumParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_EMX,
                      default=0,
                      choices=self.alignTypeList,
                      label='Alignment Type',
                      help="Is this a 2D alignment, a 3D alignment or a set of projections")
#
        form.addParam('mdFile', params.FileParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_XMIPP3,
                      label='Particles metadata file',
                      help="Select the particles Xmipp metadata file.\n"
                           "It is usually a images.xmd file result\n"
                           "from Xmipp protocols execution.")
        
        form.addParam('starFile', params.FileParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_RELION,
                      label='Star file',
                      help="Select a *_data.star file from a\n"
                           "previous Relion execution."
                           "To detect if the input particles contains alignment "
                           "information, it is required to have the "
                           "optimiser.star file corresponding to the data.star")
        
        form.addParam('ignoreIdColumn', params.BooleanParam, default=False,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_RELION,
                      label='Ignore ID column?',
                      help="Set this option to True to regenerate \n"
                           "the id's of the particles. By default \n"
                           "it is read from metadata file.        \n"
                           "This option can be useful when merging\n"
                           "different metadatas and id's are not  \n"
                           "longer unique.")
        
        form.addParam('sqliteFile', params.FileParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_SCIPION,
                      label='Particles sqlite file',
                      help="Select the particles sqlite file.\n")

        form.addParam('frealignLabel', params.LabelParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_FREALIGN,
                      label='For Frealign you need to import both stack and .par files.')  
        form.addParam('stackFile', params.FileParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_FREALIGN,
                      label='Stack file',
                      help="Select an stack file with the particles.")          
        form.addParam('parFile', params.FileParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_FREALIGN,
                      label='Param file',
                      help="Select a Frealign .par file with the refinement information.")
        form.addParam('lstFile', params.FileParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_EMAN,
                      label='Lst file',
                      help='Select a *.lst set file from EMAN2 project.')

        form.addParam('csFile', params.FileParam,
                      condition='(importFrom == %d)' % self.IMPORT_FROM_CRYOSPARC,
                      label='cs file',
                      help="Select a .cs file.\n"
                           "It is usually a .cs file result from cryoSPARC job "
                           "execution.")
        
    def _defineAcquisitionParams(self, form):
        group = ProtImportImages._defineAcquisitionParams(self, form)
        group.addParam('samplingRate', params.FloatParam,
                       label=pwutils.Message.LABEL_SAMP_RATE)

    def _insertAllSteps(self):
        importFrom = self.importFrom.get()
        ci = self.getImportClass()
        
        if ci is None:
            ProtImportImages._insertAllSteps(self)
        else:
            self._insertFunctionStep('importParticlesStep', importFrom,
                                     self.importFilePath)
            
    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        if self.importFrom == self.IMPORT_FROM_EMX:
            EmxImport = Domain.importFromPlugin('emxlib.convert', 'EmxImport',
                                                errorMsg='Emx is needed to import .emx files',
                                                doRaise=True)
            self.importFilePath = abspath(self.emxFile.get('').strip())
            return EmxImport(self, self.importFilePath,
                             self.alignTypeList[self.alignType.get()])

        elif self.importFrom == self.IMPORT_FROM_XMIPP3:
            XmippImport = Domain.importFromPlugin('xmipp3.convert', 'XmippImport',
                                                  'Xmipp is needed to import .xmd files',
                                                  doRaise=True)
            self.importFilePath = self.mdFile.get('').strip()
            return XmippImport(self, self.mdFile.get())

        elif self.importFrom == self.IMPORT_FROM_RELION:
            RelionImport = Domain.importFromPlugin('relion.convert', 'RelionImport',
                                                   errorMsg='Relion is needed to import .star files',
                                                   doRaise=True)
            self.importFilePath = self.starFile.get('').strip()
            return RelionImport(self, self.starFile.get())
        elif self.importFrom == self.IMPORT_FROM_SCIPION:
            from .dataimport import ScipionImport
            self.importFilePath = self.sqliteFile.get('').strip()
            return ScipionImport(self, self.importFilePath)    
        elif self.importFrom == self.IMPORT_FROM_FREALIGN:
            self.importFilePath = self.parFile.get('').strip()
            GrigorieffLabImportParticles = Domain.importFromPlugin(
                'cistem.convert', 'GrigorieffLabImportParticles',
                errorMsg='Cistem is needed to import .stk files',
                doRaise=True)
            return GrigorieffLabImportParticles(self, self.parFile.get(),
                                                self.stackFile.get())
        elif self.importFrom == self.IMPORT_FROM_EMAN:
            self.importFilePath = self.lstFile.get('').strip()
            EmanImport = Domain.importFromPlugin('eman2.convert', 'EmanImport',
                                                 doRaise=True)
            return EmanImport(self, self.lstFile.get())

        elif self.importFrom == self.IMPORT_FROM_CRYOSPARC:
            cryoSPARCImport = Domain.importFromPlugin('cryosparc2.convert', 'cryoSPARCImport',
                                                      'cryoSPARC is needed to import .cs files',
                                                      doRaise=True)
            self.importFilePath = self.csFile.get('').strip()
            return cryoSPARCImport(self, self.csFile.get())
        else:
            self.importFilePath = ''
            return None 
                    
    def setSamplingRate(self, imgSet):
        imgSet.setSamplingRate(self.samplingRate.get())
    
    def importParticlesStep(self, importFrom, *args):
        ci = self.getImportClass()
        ci.importParticles()
        
        summary = "Import from *%s* file:\n" % self.getEnumText('importFrom')
        summary += self.importFilePath + '\n'
        
        if self.hasAttribute('outputParticles'):
            particles = self.outputParticles
            summary += ' Particles: *%d* ' % particles.getSize()
            summary += ('(ctf=%s, alignment=%s, phaseFlip=%s)\n'
                        % (particles.hasCTF(), particles.getAlignment(),
                           particles.isPhaseFlipped()))

        # EMX files can contain only Coordinates information
        if self.hasAttribute('outputCoordinates'):
            summary += '   Coordinates: *%d* \n' % (self.outputCoordinates.getSize())

        # EMX files can contain only Coordinates information
        if self.hasAttribute('outputMicrographs'):
            summary += '   Micrographs: *%d* \n' % (self.outputMicrographs.getSize())
        
        if self.copyFiles:
            summary += '\n_WARNING_: Binary files copied into project (extra disk space)'
            
        self.summaryVar.set(summary)

    def _validateFileExtension(self):
        """ Simple check about the expected file extension. """
        # Since 'from files' is index 0, we need to subtract 1 to
        # get the proper index
        i = self.importFrom.get() - 1
        ext = self.importExts[i]
        if not self.importFilePath.endswith(ext):
            return ["Expected *%s* file extension for importing from *%s*" %
                    (ext, self.importFormats[i])]
        else:
            return []

    def _validate(self):
        ci = self.getImportClass()
        if ci is None:
            return ProtImportImages._validate(self)
        else:
            errors = self._validateFileExtension()
            if errors:
                return errors
            else:
                return ci.validateParticles()
    
    def _summary(self):
        if self.importFrom == self.IMPORT_FROM_FILES:
            return ProtImportImages._summary(self)
        else:
            return [self.summaryVar.get('')]


class ProtImportAverages(ProtImportParticles):
    """Protocol to import a set of averages to the project"""
    _label = 'import averages'
    _outputClassName = 'SetOfAverages'    

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formas such as: xmipp3, eman2, relion...etc.
        """
        choices = ProtImportImages._getImportChoices(self)
        return choices
            
    def _defineAcquisitionParams(self, form):
        form.addParam('samplingRate', params.FloatParam, default=1.,
                      label=pwutils.Message.LABEL_SAMP_RATE)
        group = ProtImportImages._defineAcquisitionParams(self, form)
        group.expertLevel.set(pwcts.LEVEL_ADVANCED)
