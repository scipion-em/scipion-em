# **************************************************************************
# *
# * Authors:     Pablo Conesa (pconesa@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
""" This module is for Actual wizards to be able to be discovered from the
Domain. wizard.py is left for wizard models and base classes."""
import decimal
import os
import requests

import pyworkflow.object as pwobj
from pyworkflow.gui import dialog
import pyworkflow.wizard as pwizard
from pyworkflow.gui.tree import ListTreeProviderString

import pwem.convert as emconv

from pwem.wizards.wizard import EmWizard, FormulaDialog
import pwem.protocols as emprot
import pwem.objects as emobj


class ImportAcquisitionWizard(EmWizard):
    _targets = [(emprot.ProtImportImages, ['acquisitionWizard'])]

    def show(self, form, *params):
        try:
            acquisitionInfo = form.protocol.loadAcquisitionInfo()
            if isinstance(acquisitionInfo, dict):
                # If acquisitionInfo is None means something is wrong.
                # Now, let's try to show a meaningful error message.
                self._setAcquisition(form, acquisitionInfo)
            else:
                # If not dict, it should be an error message
                dialog.showError("Input error", acquisitionInfo, form.root)

        except FileNotFoundError as e:
            dialog.showInfo("File not found", "Metadata file with acquisition not found.\n\n %s" % e, form.root)

    @classmethod
    def _setAcquisition(cls, form, acquisitionInfo):
        """ Ask whether to set the AcquisitionInfo to the protocol parameters.
        Params:
            acquisitionInfo: Should be a dictionary with acquisition values.
                If None, show an error.
        """
        msg = ''
        for k, v in acquisitionInfo.items():
            msg += '%s = %s\n' % (k, v)
        msg += '\n*Do you want to use detected acquisition values?*'
        response = dialog.askYesNo("Import acquisition",
                                   msg, form.root)
        if response:
            prot = form.protocol
            comment = ''

            for k, v in acquisitionInfo.items():
                if prot.hasAttribute(k):
                    form.setVar(k, v)
                else:
                    comment += "%s = %s\n" % (k, v)
            if comment:
                prot.setObjComment(comment)


class ImportCoordinatesBoxSizeWizard(pwizard.Wizard):
    _targets = [(emprot.ProtImportCoordinates, ['boxSize']),
                (emprot.ProtImportCoordinatesPairs, ['boxSize'])]

    @classmethod
    def _getBoxSize(cls, protocol):
        return protocol.getDefaultBoxSize()

    @classmethod
    def show(cls, form, *params):
        form.setVar('boxSize', cls._getBoxSize(form.protocol))


class ImportOriginVolumeWizard(pwizard.Wizard):
    _targets = [(emprot.ProtImportVolumes, ['x', 'y', 'z'])]

    def show(self, form, *params):
        protocol = form.protocol
        filesPath = protocol.filesPath.get()
        filesPattern = protocol.filesPattern.get()
        if filesPattern:
            fullPattern = os.path.join(filesPath, filesPattern)
        else:
            fullPattern = filesPath

        sampling = protocol.samplingRate.get()
        for fileName, fileId in protocol.iterFiles():
            inputVol = emobj.Volume()
            inputVol.setFileName(fileName)
            if ((str(fullPattern)).endswith('mrc') or
                    (str(fullPattern)).endswith('map')):
                ccp4header = emconv.Ccp4Header(fileName, readHeader=True)
                x, y, z = ccp4header.getOrigin(changeSign=True)  # In Angstroms
            else:
                x, y, z = self._halfOriginCoordinates(inputVol, sampling)

            form.setVar('x', x)
            form.setVar('y', y)
            form.setVar('z', z)

    @classmethod
    def _halfOriginCoordinates(cls, volume, sampling):
        xdim, ydim, zdim = volume.getDim()
        if zdim > 1:
            zdim = zdim / 2.
        x = xdim / 2. * sampling
        y = ydim / 2. * sampling
        z = zdim * sampling
        return x, y, z


class ChangeOriginSamplingWizard(pwizard.Wizard):
    _targets = [(emprot.ProtOrigSampling, ['x', 'y', 'z', 'samplingRate'])]

    def show(self, form, *params):
        protocol = form.protocol
        vol = protocol.inVolume.get()
        fullPattern = vol.getFileName()
        sampling = vol.getSamplingRate()
        if ((str(fullPattern)).endswith('mrc') or
                (str(fullPattern)).endswith('map')):
            ccp4header = emconv.Ccp4Header(fullPattern, readHeader=True)
            x, y, z = ccp4header.getOrigin(changeSign=True)  # In Angstroms
        else:
            x, y, z = \
                ImportOriginVolumeWizard._halfOriginCoordinates(vol, sampling)

        form.setVar('x', round(x, 3))
        form.setVar('y', round(y, 3))
        form.setVar('z', round(z, 3))
        form.setVar('samplingRate', round(sampling, 3))


class GetStructureChainsWizard(pwizard.Wizard):
    """Load an atomic structure, parse chain related information as
       name, number of residues, list of aminoacids (or other residues)"""
    _targets = [(emprot.ProtImportSequence, ['inputStructureChain'])
                # NOTE: be careful if you change this class since
                # chimera-wizard inherits from it.
                # (ChimeraModelFromTemplate, ['inputStructureChain'])
                # (atomstructutils, ['inputStructureChain'])
                ]

    @classmethod
    def getModelsChainsStep(cls, protocol):
        """ Returns (1) list with the information
           {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
           (2) list with residues, position and chain (modelsFirstResidue)"""
        structureHandler = emconv.AtomicStructHandler()
        fileName = ""
        if hasattr(protocol, 'pdbId'):
            if protocol.pdbId.get() is not None:
                pdbID = protocol.pdbId.get()
                url = "https://www.rcsb.org/structure/"
                URL = url + ("%s" % pdbID)
                try:
                    response = requests.get(URL)
                except:
                    raise Exception("Cannot connect to PDB server")
                if (response.status_code >= 400) and (response.status_code < 500):
                    raise Exception("%s is a wrong PDB ID" % pdbID)
                fileName = structureHandler.readFromPDBDatabase(
                    os.path.basename(pdbID), dir="/tmp/")
            else:
                fileName = protocol.pdbFile.get()
        else:
            if protocol.pdbFileToBeRefined.get() is not None:
                fileName = os.path.abspath(protocol.pdbFileToBeRefined.get(
                ).getFileName())

        structureHandler.read(fileName)
        structureHandler.getStructure()
        # listOfChains, listOfResidues = structureHandler.getModelsChains()
        return structureHandler.getModelsChains()

    def editionListOfChains(self, listOfChains):
        self.chainList = []
        for model, chainDic in listOfChains.items():
            for chainID, lenResidues in chainDic.items():
                self.chainList.append(
                    '{"model": %d, "chain": "%s", "residues": %d}' %
                    (model, str(chainID), lenResidues))

    def show(self, form, *params):
        protocol = form.protocol
        try:
            listOfChains, listOfResidues = self.getModelsChainsStep(protocol)
        except Exception as e:
            print("ERROR: ", e)
            return

        self.editionListOfChains(listOfChains)
        finalChainList = []
        for i in self.chainList:
            finalChainList.append(pwobj.String(i))
        provider = ListTreeProviderString(finalChainList)
        dlg = dialog.ListDialog(form.root, "Model chains", provider,
                                "Select one of the chains (model, chain, "
                                "number of chain residues)")
        form.setVar('inputStructureChain', dlg.values[0].get())

class PythonFormulaeWizard(pwizard.Wizard):
    """Assist in the creation of python formula to be evaluated. In Steps"""
    _targets = [(emprot.ProtSetFilter, ['formula']),
                (emprot.ProtSetEditor, ['formula'])]

    def show(self, form, *params):

        d = FormulaDialog(form.root, form.protocol.inputSet.get(), formula=form.protocol.formula.get())

        # If accepted
        if d.resultYes():
            form.setVar('formula',d.getFormula())
            
class PythonTopRankWizard(pwizard.Wizard):
    """Assist in the creation of python formula to be evaluated. In Steps"""
    _targets = [(emprot.ProtSetFilter, ['rankingField'])]

    def getInputAttributes(self, form):
        attrNames=[]
        item = form.protocol.inputSet.get().getFirstItem()
        for key, attr in item.getAttributesToStore():
            attrNames.append(key)
        return attrNames

    def show(self, form, *params):
        attrsList = self.getInputAttributes(form)
        finalAttrsList = []
        for i in attrsList:
            finalAttrsList.append(pwobj.String(i))
        provider = ListTreeProviderString(finalAttrsList)

        dlg = dialog.ListDialog(form.root, "Filter set", provider,
                                "Select one of the attributes")
        form.setVar('rankingField', dlg.values[0].get())
