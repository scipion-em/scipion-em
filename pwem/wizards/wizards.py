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
import os, json
import requests

import pyworkflow.object as pwobj
from pyworkflow.gui import dialog
import pyworkflow.wizard as pwizard
from pyworkflow.gui.tree import ListTreeProviderString

import pwem.convert as emconv

from pwem.wizards.wizard import EmWizard, FormulaDialog, ColorScaleWizardBase, VariableWizard
import pwem.protocols as emprot
import pwem.viewers as emview
import pwem.objects as emobj
from pwem.constants import RESIDUES3TO1, RESIDUES1TO3


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
  # WARNING: this wizard is deprecated. Instead, SelectChainWizard must be used, where inputs, targets and outputs
  # can be specified.
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


class SelectChainWizard(VariableWizard):
    '''Opens the input AtomStruct and allows you to select one of the present chains'''
    _targets, _inputs, _outputs = [], {}, {}

    @classmethod
    def getModelsChainsStep(cls, protocol, inputObj):
        """ Returns (1) list with the information
           {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
           (2) list with residues, position and chain (modelsFirstResidue)"""
        structureHandler = emconv.AtomicStructHandler()
        if type(inputObj) == str:
            if os.path.exists(inputObj):
                fileName = inputObj
            else:
                pdbID = inputObj
                url = "https://www.rcsb.org/structure/"
                URL = url + ("%s" % pdbID)
                try:
                    response = requests.get(URL)
                except:
                    raise Exception("Cannot connect to PDB server")
                if (response.status_code >= 400) and (response.status_code < 500):
                    raise Exception("%s is a wrong PDB ID" % pdbID)
                fileName = structureHandler.readFromPDBDatabase(os.path.basename(pdbID), dir="/tmp/")

        elif str(type(inputObj).__name__) == 'SchrodingerAtomStruct':
            fileName = os.path.abspath(inputObj.convert2PDB())
        else:
            fileName = os.path.abspath(inputObj.getFileName())

        structureHandler.read(fileName)
        structureHandler.getStructure()
        return structureHandler.getModelsChains()

    def editionListOfChains(self, listOfChains):
        chainList = []
        for model, chainDic in listOfChains.items():
            for chainID, lenResidues in chainDic.items():
                chainList.append(
                    '{"model": %d, "chain": "%s", "residues": %d}' %
                    (model, str(chainID), lenResidues))
        return chainList

    def show(self, form, *params):
        inputParams, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        try:
            inputObj = getattr(protocol, inputParams[0]).get()
            listOfChains, listOfResidues = self.getModelsChainsStep(protocol, inputObj)
        except Exception as e:
            print("ERROR: ", e)
            return

        chainList = self.editionListOfChains(listOfChains)
        finalChainList = []
        for i in chainList:
            finalChainList.append(pwobj.String(i))
        provider = ListTreeProviderString(finalChainList)
        dlg = dialog.ListDialog(form.root, "Model chains", provider,
                                "Select one of the chains (model, chain, "
                                "number of chain residues)")
        form.setVar(outputParam[0], dlg.values[0].get())

SelectChainWizard().addTarget(protocol=emprot.ProtImportSequence,
                              targets=['inputStructureChain'],
                              inputs=[['pdbId', 'pdbFile']],
                              outputs=['inputStructureChain'])


class SelectResidueWizard(SelectChainWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def editionListOfResidues(self, modelsFirstResidue, model, chain):
      residueList = []
      for modelID, chainDic in modelsFirstResidue.items():
        if int(model) == modelID:
          for chainID, seq_number in chainDic.items():
            if chain == chainID:
              for i in seq_number:
                residueList.append('{"index": %d, "residue": "%s"}' % (i[0], str(i[1])))
      return residueList

    def getResidues(self, form, inputObj, chainStr):
      protocol = form.protocol

      if type(inputObj) == str:
          # Select the residues if the input structure parameter is a str (PDB id or PDB file)
          structureHandler = emconv.AtomicStructHandler()
          if os.path.exists(inputObj):
              fileName = inputObj
          else:
              pdbID = inputObj
              url = "https://www.rcsb.org/structure/"
              URL = url + ("%s" % pdbID)
              try:
                response = requests.get(URL)
              except:
                raise Exception("Cannot connect to PDB server")
              if (response.status_code >= 400) and (response.status_code < 500):
                raise Exception("%s is a wrong PDB ID" % pdbID)

              fileName = structureHandler.readFromPDBDatabase(os.path.basename(pdbID), dir="/tmp/")

          structureHandler.read(fileName)
          structureHandler.getStructure()
          modelsLength, modelsFirstResidue = structureHandler.getModelsChains()

          struct = json.loads(chainStr)  # From wizard dictionary
          chain, model = struct["chain"].upper().strip(), int(struct["model"])

          residueList = self.editionListOfResidues(modelsFirstResidue, model, chain)
          finalResiduesList = []
          for i in residueList:
              finalResiduesList.append(emobj.String(i))

      elif issubclass(type(inputObj), emobj.AtomStruct):
          try:
            modelsLength, modelsFirstResidue = self.getModelsChainsStep(protocol, inputObj)
          except Exception as e:
            print("ERROR: ", e)
            return
          struct = json.loads(chainStr)  # From wizard dictionary
          chain, model = struct["chain"].upper().strip(), int(struct["model"])

          residueList = self.editionListOfResidues(modelsFirstResidue, model, chain)
          finalResiduesList = []
          for i in residueList:
            finalResiduesList.append(emobj.String(i))

      elif issubclass(type(inputObj), emobj.Sequence) or str(type(inputObj).__name__) == 'SequenceVariants':
          finalResiduesList = []
          for i, res in enumerate(inputObj.getSequence()):
            if res in RESIDUES1TO3:
                res3 = RESIDUES1TO3[res]
            else:
                res3 = res
            stri = '{"index": %s, "residue": "%s"}' % (i + 1, res3)
            finalResiduesList.append(emobj.String(stri))

      return finalResiduesList

    def getSequence(self, finalResiduesList, idxs):
        roiStr, inSeq = '', False
        for residue in finalResiduesList:
            resDic = json.loads(residue.get())
            if resDic['residue'] in RESIDUES3TO1:
                if resDic['index'] == idxs[0]:
                  inSeq = True
                if resDic['index'] == idxs[-1]:
                  inSeq = False
                  roiStr += RESIDUES3TO1[resDic['residue']]
                  break

                if inSeq:
                    roiStr += RESIDUES3TO1[resDic['residue']]
        return roiStr

    def show(self, form, *params):
      inputParams, outputParam = self.getInputOutput(form)
      protocol = form.protocol
      inputObj = getattr(protocol, inputParams[0]).get()
      if len(inputParams) < 2:
          # For sequence objects, with no chain
          chainStr = None
      else:
          chainStr = getattr(protocol, inputParams[1]).get()
      finalResiduesList = self.getResidues(form, inputObj, chainStr)

      provider = ListTreeProviderString(finalResiduesList)
      dlg = dialog.ListDialog(form.root, "Chain residues", provider,
                              "Select one residue (residue number, "
                              "residue name)")

      idxs = [json.loads(dlg.values[0].get())['index'], json.loads(dlg.values[-1].get())['index']]
      roiStr = self.getSequence(finalResiduesList, idxs)


      intervalStr = '{"index": "%s-%s", "residues": "%s"}' % (idxs[0], idxs[1], roiStr)
      form.setVar(outputParam[0], intervalStr)


class PythonFormulaeWizard(pwizard.Wizard):
    """Assist in the creation of python formula to be evaluated. In Steps"""
    _targets = [(emprot.ProtSetFilter, ['formula']),
                (emprot.ProtSetEditor, ['formula'])]

    def show(self, form, *params):

        d = FormulaDialog(form.root, form.protocol.inputSet.get(), formula=form.protocol.formula.get())

        # If accepted
        if d.resultYes():
            form.setVar('formula',d.getFormula())


class PythonFormulaWizardX1(pwizard.Wizard):
    """Assist in the creation of python formula to be evaluated. In Steps"""
    _targets = [(emprot.ProtMathematicalOperator, ['attribute1'])]

    def show(self, form, *params):

        d = FormulaDialog(form.root, form.protocol.inputSet1.get(), formula=form.protocol.attribute1.get())

        # If accepted
        if d.resultYes():
            form.setVar('attribute1', d.getFormula())


class PythonFormulaWizardX2(pwizard.Wizard):
    """Assist in the creation of python formula to be evaluated. In Steps"""
    _targets = [(emprot.ProtMathematicalOperator, ['attribute2'])]

    def show(self, form, *params):

        d2 = FormulaDialog(form.root, form.protocol.inputSet2.get(), formula=form.protocol.attribute2.get())

        # If accepted
        if d2.resultYes():
            form.setVar('attribute2', d2.getFormula())


class PythonTopRankWizard(pwizard.Wizard):
    """Assist in the creation of python formula to be evaluated. In Steps"""
    _targets = [(emprot.ProtSetFilter, ['rankingField'])]

            

class SelectAttributeWizard(VariableWizard):
    """Wizard to select attributes stored in a scipion object or set """
    _targets, _inputs, _outputs = [], {}, {}

    def getFirstItem(self, form, inputParam):
        inputPointer = getattr(form.protocol, inputParam)
        if issubclass(inputPointer.__class__, pwobj.PointerList):
            inputPointer = inputPointer[0]

        inputSet = inputPointer.get()
        if issubclass(inputSet.__class__, pwobj.Set):
            item = inputSet.getFirstItem()
        elif issubclass(inputSet.__class__, pwobj.Object):
            item = inputSet
        return item

    def getInputAttributes(self, form, inputParam):
      attrNames = []
      item = self.getFirstItem(form, inputParam[0])
      for key, attr in item.getAttributesToStore():
        attrNames.append(key)
      return attrNames

    def show(self, form, *params):
      inputParam, outputParam = self.getInputOutput(form)
      attrsList = self.getInputAttributes(form, inputParam)
      finalAttrsList = []
      for i in attrsList:
        finalAttrsList.append(pwobj.String(i))
      provider = ListTreeProviderString(finalAttrsList)
      dlg = dialog.ListDialog(form.root, "Filter set", provider,
                            "Select one of the attributes")
      form.setVar(outputParam[0], dlg.values[0].get())


class ColorScaleWizardRMSD(ColorScaleWizardBase):
    _targets = ColorScaleWizardBase.defineTargets(emview.ChimeraAttributeViewer)


#Defining target for the SelectAttributeWizard
SelectAttributeWizard().addTarget(protocol=emprot.ProtSetFilter,
                                  targets=['rankingField'],
                                  inputs=['inputSet'],
                                  outputs=['rankingField'])
