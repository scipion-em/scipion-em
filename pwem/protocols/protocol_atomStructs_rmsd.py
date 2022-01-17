# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os, shutil
import numpy as np
import itertools

from pyworkflow.protocol.params import FloatParam, PointerParam, EnumParam, BooleanParam, \
    STEPS_PARALLEL, StringParam
from pyworkflow.object import Float, String
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb, toCIF, AtomicStructHandler, addScipionAttribute
from pwem.objects import AtomStruct, SetOfAtomStructs

class ProtRMSDAtomStructs(EMProtocol):
    """
    Protocol to calculate the RMSD between all pairs of atom structures in a set of them.
    It calculates the overall RMSD for all of them and for each of their residues to validate
    their associated volume
    """

    _label = 'RMSD validate map'
    _ATTRNAME = 'perResidueRMSD'
    _OUTNAME = 'outputAtomStructs'
    _possibleOutputs = {_OUTNAME: SetOfAtomStructs}

    averageRMSDs = []
    combinedResRMSDs = []

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructureSet', PointerParam,
                      pointerClass="SetOfAtomStructs", allowsNull=False,
                      label='Input atomic structures.',
                      help="Set the atomic structure to be processed.\n"
                           "Supported formats are PDB or mmCIF; this last one"
                           " is especially useful for very large structures.")
        form.addParam('chains', StringParam, default='',
                      label='Chains to perform the RMSD on: ',
                      help='Comma-separated chains to perform the RMSD on.\nIf empty, all chains will be used')

        group = form.addGroup('Atoms')
        group.addParam('considerAtoms', EnumParam, default=0,
                       choices=['CA', 'BackBone', 'All'], label='Atoms to take into account: ',
                       help='Atoms to use in the RMSD calculation')
        group.addParam('weightbb', FloatParam, default=1.0,
                       label='Backbone weight: ', condition='considerAtoms==2',
                       help='Set the resolution of the input volume.')
        group.addParam('hydrogen', BooleanParam, default=False,
                       label='Consider hydrogens: ',
                       help='Consider hydrogens into calculation')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        convId = self._insertFunctionStep('convertInputStep', prerequisites=[])

        rmsdIds = []
        for comboFn in itertools.combinations(self.getInputFileNames(), 2):
            rmsdIds.append(self._insertFunctionStep('calculateRMSDStep', comboFn, prerequisites=[convId]))

        self._insertFunctionStep('createOutputStep', prerequisites=rmsdIds)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self):
        for inFn in self.getInputFileNames():
            if not inFn.endswith('.pdb'):
                name, ext = os.path.splitext(inFn)
                toPdb(os.path.abspath(inFn), self._getExtraPath(os.path.basename(inFn).replace(ext, '.pdb')))
            else:
                shutil.copy(inFn, self._getExtraPath(os.path.basename(inFn)))


    def calculateRMSDStep(self, combo):
        pdbFiles = self.getConvertedFile(combo[0]), self.getConvertedFile(combo[1])

        if self.chains.get().strip() == '':
            chainlist = None
        else:
            chainlist = self.chains.get().split(',')
        P, Pres, self.Preslist = get_coordinates(pdbFiles[0], self.hydrogen.get(), chainlist, self.getRMSDAtoms(),
                                            self.weightbb.get())
        Q, Qres, Qreslist = get_coordinates(pdbFiles[1], self.hydrogen.get(), chainlist, self.getRMSDAtoms(),
                                            self.weightbb.get())
        # checking that files have same number of residues
        if (len(Pres) != len(Qres)):
            exit("Error: files [%s, %s] do not have same number of residues, %i vs. %i" % (
            pdbFiles[0], pdbFiles[1], len(Pres), len(Qres)))
        else:
            self.numRes = len(Qres)

        # Calculate per-residue RMSD
        per_res_rmsd = []
        for i in range(len(Pres)):
            res1 = np.asarray(Pres[i])
            res2 = np.asarray(Qres[i])
            per_res_rmsd.append(rmsd(res1, res2))

        # Calculate overall RMSD
        rmsdval = rmsd(P, Q)
        self.averageRMSDs.append(rmsdval)
        self.combinedResRMSDs.append(per_res_rmsd)
        print("RMSD %s to %s: %.4f" % (os.path.basename(pdbFiles[0]), os.path.basename(pdbFiles[1]), rmsdval))
    
    def createOutputStep(self):
        outStructFileBase = self._getPath('{}.cif')
        overFinalRMSD = np.mean(self.averageRMSDs)
        print("Overall RMSD:", overFinalRMSD)
        rmsdAttrDic = self.getRMSDAttributeDic()

        ASH = AtomicStructHandler()
        outSet = SetOfAtomStructs.create(self._getPath())
        for AS in self.inputStructureSet.get():
            outStructFileName = outStructFileBase.format(os.path.splitext(os.path.basename(AS.getFileName()))[0])
            inpAS = toCIF(AS.getFileName(), self._getTmpPath('inputStruct.cif'))
            cifDic = ASH.readLowLevel(inpAS)
            cifDic = addScipionAttribute(cifDic, rmsdAttrDic, self._ATTRNAME)
            ASH._writeLowLevel(outStructFileName, cifDic)

            outAS = AS.clone()
            outAS.setFileName(outStructFileName)
            outSet.append(outAS.clone())

        outSet.overallRMSD = Float(overFinalRMSD)
        self._defineOutputs(**{self._OUTNAME:outSet})


    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        try:
            summary.append('Overall RMSD: {:.4f}\n'.format(float(self.outputAtomStructs.overallRMSD)))
        except:
            summary = ["Overall RMSD not yet computed"]
        return summary

    # --------------------------- UTILS functions ----------------------------
    def getInputFileNames(self):
        fns = []
        for AS in self.inputStructureSet.get():
            fns.append(AS.getFileName())
        return fns
    
    def getInputVolume(self):
        return self.inputStructureSet.get().getFirstItem().getVolume()

    def getConvertedFile(self, inFile):
        name, ext = os.path.splitext(os.path.basename(inFile))
        for file in os.listdir(self._getExtraPath()):
            if name in file:
                return self._getExtraPath(file)

    def getRMSDAtoms(self):
        if self.considerAtoms.get() == 0:
            return 'CA'
        elif self.considerAtoms.get() == 1:
            return 'full'
        else:
            return None

    def getRMSDAttributeDic(self):
        '''Return a dictionary with {spec: value}
        "spec" should be a chimera specifier'''
        combinedResAvgs = []
        # get average per-residue RMSD
        for i in range(self.numRes):
            rmsdval = 0
            for struct in self.combinedResRMSDs:
                rmsdval += struct[i]
            # this gives us mean RMSD for this residue
            rmsdval /= len(self.combinedResRMSDs)
            combinedResAvgs.append(rmsdval)

        attrDic = {}
        for i in range(len(combinedResAvgs)):
            resId = '{}:{}'.format(self.Preslist[i][4], self.Preslist[i][5:].strip())
            attrDic[resId] = str(round(combinedResAvgs[i], 4))
        return attrDic


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i] - w[i]) ** 2.0 for i in range(D)])
    return np.sqrt(rmsd / N)


def get_coordinates(filename, hydrogen, chain=None, backbone=None, weight=1.0):
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.
    """
    atomlist = ["C", "N", "O", "S", "P"]
    if hydrogen:
        atomlist.append("H")

    x_column = None
    V = []

    # resarray contains the coordinates
    resarray = []

    # reslist contains the residue name and number
    reslist = []

    prevres = ""
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            #			if line.startswith("TER") or line.startswith("END"):
            if line.startswith("END"):
                break
            if line.startswith("ATOM"):
                # establish backbone atoms
                backatoms = ['N', 'CA', 'C', 'O']

                # if specific chain is specified, only process that chain
                if chain is not None:
                    if line[21] not in chain: continue

                # see if it's a nucleic acid
                if backbone == 'CA': backatoms = ['CA']
                if line[17:20].strip() in ['DA', 'DC', 'DG', 'DT', 'DI', 'A', 'C', 'G', 'U', 'I']:
                    # add P for RNA/DNA
                    backatoms = ["P", "OP1", "OP2", "O3'", "O5'", "C3'", "C4'", "C5'"]
                    if backbone == 'CA': backatoms = ["C3'", "C4'"]
                tokens = line.split()
                # Try to get the residue
                try:
                    res = line[17:27]
                    # if moved to next residue:
                    if res != prevres:
                        reslist.append(res)
                        if prevres != "":
                            resarray.append(presarray)
                        presarray = []
                        prevres = res
                except:
                    exit("Error parsing residue type for the following line: \n%s" % line)
                # Try to get the atomtype
                try:
                    atom = tokens[-1][0]
                    if atom not in atomlist: continue

                    # check if a backbone atom
                    multiplier = weight
                    if line[12:16].strip() not in backatoms:
                        if backbone is not None: continue
                        multiplier = 1.0 / weight

                except:
                    exit("Error parsing atomtype for the following line: \n%s" % line)

                try:
                    x = line[30:38]
                    y = line[38:46]
                    z = line[46:54]
                    V.append(np.asarray([x, y, z], dtype=float) * multiplier)
                    presarray.append(np.asarray([x, y, z], dtype=float) * multiplier)
                except:
                    exit("Error parsing input for the following line: \n%s" % line)
        resarray.append(presarray)

    V = np.asarray(V)
    return V, resarray, reslist