# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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
#
# cif to pdb conversion is based on  cif2pdb.py by Spencer Bliven
# <spencer.bliven@gmail.com>
#
# see http://biopython.org/DIST/docs/tutorial/Tutorial.html for a description
# on the object "structure" and other Bio.xxxx modules

from Bio.PDB.Dice import ChainSelector
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO, MMCIFIO
from Bio.PDB.mmcifio import mmcif_order
from Bio.PDB import Superimposer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import Entity
from Bio.PDB import PDBList
from collections import OrderedDict
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Polypeptide import three_to_one
from Bio.Seq import Seq
from pwem import Plugin, MAXIT_HOME
from pwem.convert.transformations import translation_from_matrix
import mmap
import re
import os
import numpy
import pyworkflow.utils as pwutils
import shutil


class OutOfChainsError(Exception):
    pass


class scipionMMCIFIO(MMCIFIO):

    # this class was needed before including maxit
    # as converter.
    # TODO: Delete it after a more extensive testing period
    # DATE mar feb 18 17:48:57 CET 2020

    def _save_dict_delete(self, out_file):
        # Form dictionary where key is first part of mmCIF key and value is list
        # of corresponding second parts
        key_lists = {}
        for key in self.dic:
            if key == "data_":
                data_val = self.dic[key]
            else:
                s = re.split(r"\.", key)
                if len(s) == 2:
                    if s[0] in key_lists:
                        key_lists[s[0]].append(s[1])
                    else:
                        key_lists[s[0]] = [s[1]]
                else:
                    raise ValueError("Invalid key in mmCIF dictionary: " + key)

        # Re-order lists if an order has been specified
        # Not all elements from the specified order are necessarily present
        for key, key_list in key_lists.items():
            # P3
            # for key, key_list in list(key_lists.items()):
            if key in mmcif_order:
                inds = []
                for i in key_list:
                    try:
                        inds.append(mmcif_order[key].index(i))
                    # Unrecognised key - add at end
                    except ValueError:
                        inds.append(len(mmcif_order[key]))
                key_lists[key] = [k for _, k in sorted(zip(inds, key_list))]

        # Write out top data_ line
        if data_val:
            out_file.write("data_" + data_val + "\n#\n")
            # ESCRIBIR POLYSEQTABLE
            out_file.write("""loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero\n#\n""")
            total_seqs = []
            counter_chain = 1
            for model in self.structure:
                for chain in model:
                    min_val = int(chain.get_unpacked_list()[0].id[1])
                    # max_val = min_val + len(chain.get_unpacked_list())
                    #                    not_are = range(1, min_val-1)
                    #                    for i in range(min_val)[:-1]:
                    #                        not_are.append(i + 1)
                    counter = 1
                    if min_val > 1:
                        for counter in range(1, min_val):
                            if len(chain.get_unpacked_list()[0].resname.strip()) == 3:
                                total_seqs.append((counter_chain, str(counter), "XAA", "n"))
                            elif len(chain.get_unpacked_list()[0].resname.strip()) == 2:
                                total_seqs.append((counter_chain, str(counter), "DX", "n"))
                            elif len(chain.get_unpacked_list()[0].resname.strip()) == 1:
                                total_seqs.append((counter_chain, str(counter), "X", "n"))
                            counter += 1
                    for residue in chain:
                        if len(chain.get_unpacked_list()[0].resname.strip()) == 3 and \
                                is_aa(residue.get_resname(), standard=True):  # aminoacids
                            total_seqs.append((counter_chain, residue.id[1], residue.get_resname(), "n"))
                        elif len(chain.get_unpacked_list()[0].resname.strip()) == 2 and \
                                residue.get_resname()[1] in ['A', 'C', 'G', 'T']:
                            total_seqs.append((counter_chain, counter, residue.get_resname(), "n"))
                        elif len(chain.get_unpacked_list()[0].resname.strip()) == 1 and \
                                residue.get_resname() in ['A', 'C', 'G', 'U']:
                            total_seqs.append((counter_chain, counter, residue.get_resname(), "n"))
                        counter += 1
                    counter_chain += 1
            for item in total_seqs:
                item = list(item)
                out_file.write("""%s %s %s  %s\n""" %
                               (str(item[0]), item[1], item[2], item[3]))

        for key, key_list in key_lists.items():
            #  for key, key_list in list(key_lists.items()):
            # Pick a sample mmCIF value, which can be a list or a single value
            sample_val = self.dic[key + "." + key_list[0]]
            n_vals = len(sample_val)
            # Check the mmCIF dictionary has consistent list sizes
            for i in key_list:
                val = self.dic[key + "." + i]
                if (isinstance(sample_val, list) and (isinstance(val, str) or len(val) != n_vals)) or (
                        isinstance(sample_val, str) and isinstance(val, list)):
                    raise ValueError("Inconsistent list sizes in mmCIF dictionary: " + key + "." + i)
            # If the value is a single value, write as key-value pairs
            if isinstance(sample_val, str):
                m = 0
                # Find the maximum key length
                for i in key_list:
                    if len(i) > m:
                        m = len(i)
                for i in key_list:
                    out_file.write(
                        "{k: <{width}}".format(k=key + "." + i, width=len(key) + m + 4) + self._format_mmcif_col(
                            self.dic[key + "." + i], len(self.dic[key + "." + i])) + "\n")
            # If the value is a list, write as keys then a value table
            elif isinstance(sample_val, list):
                out_file.write("loop_\n")
                col_widths = {}
                # Write keys and find max widths for each set of values
                for i in key_list:
                    out_file.write(key + "." + i + "\n")
                    col_widths[i] = 0
                    for val in self.dic[key + "." + i]:
                        len_val = len(val)
                        # If the value requires quoting it will add 2 characters
                        if self._requires_quote(val) and not self._requires_newline(val):
                            len_val += 2
                        if len_val > col_widths[i]:
                            col_widths[i] = len_val
                # Technically the max of the sum of the column widths is 2048

                # Write the values as rows
                for i in range(n_vals):
                    for col in key_list:
                        out_file.write(self._format_mmcif_col(self.dic[key + "." + col][i], col_widths[col] + 1))
                    out_file.write("\n")
            else:
                raise ValueError("Invalid type in mmCIF dictionary: " + str(type(sample_val)))
            out_file.write("#\n")


class AtomicStructHandler:
    """ Class that contain utilities to handle pdb/cif files"""
    PDB = 0
    CIF = 1

    def __init__(self, fileName=None, permissive=1):
        # The PERMISSIVE flag indicates that a number of common problems
        # associated with PDB files will be ignored
        self.permissive = permissive
        self.pdbParser = None
        self.cifParser = None
        self.ioPDB = None
        self.ioCIF = None
        self.structure = None
        self._readDone = False
        if fileName is not None:
            self.read(fileName)

    def readFromPDBDatabase(self, pdbID, dir=None, type='mmCif'):
        """
        Retrieve structure from PDB
        :param pdbID:
        :param dir: save structure in this directory
        :param type:  mmCif or pdb
        :return: filename with pdb file
        """
        if dir is None:
            dir = os.getcwd()
        pdbl = PDBList(pdb=dir)
        fileName = pdbl.retrieve_pdb_file(pdbID, pdir=dir, file_format=type)
        return os.path.abspath(fileName)

    def checkLabelInFile(self, fileName, label):
        with open(fileName, "r") as f:
            s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            if s.find(label) != -1:
                return True
            else:
                return False

    def getStructure(self):
        """return structure information, model, chain,
           residues, atoms..."""
        return self.structure

    def read(self, fileName):
        """ Read and parse file."""
        # biopython assigns an ID to any read structure
        structure_id = os.path.basename(fileName)
        structure_id = structure_id[:4] if len(structure_id) > 4 else "1xxx"

        if fileName.endswith(".pdb") or fileName.endswith(".ent"):
            if self.pdbParser is None:
                self.pdbParser = PDBParser(PERMISSIVE=self.permissive)
            parser = self.pdbParser
            self.type = self.PDB
        else:
            if self.cifParser is None:
                self.cifParser = MMCIFParser()
            parser = self.cifParser
            self.type = self.CIF

        self.structure = parser.get_structure(structure_id, fileName)
        self._readDone = True

    def checkRead(self):
        if self._readDone:
            return True
        else:
            print("you must read the pdb file first")
            exit(0)

    def getModelsChains(self):
        """
        given an atomic structure returns two dictionaries:
            (1) for all models and respective chains (chainID and length of residues)
            (2) for each chain list of residues
        """
        self.checkRead()
        listOfChains = OrderedDict()
        listOfResidues = OrderedDict()

        for model in self.structure:
            chainDicLength = OrderedDict()
            chainDicFirstResidue = OrderedDict()
            for chain in model:
                if len(chain.get_unpacked_list()[0].resname.strip()) == 1:  # RNA
                    seq = list()
                    seq_number = list()
                    for residue in chain:
                        if residue.get_resname() in ['A', 'C', 'G', 'U']:
                            seq.append(residue.get_resname())
                        else:
                            seq.append("X")
                        seq_number.append((residue.get_id()[1], residue.get_resname()))
                elif len(chain.get_unpacked_list()[0].resname.strip()) == 2:  # DNA
                    seq = list()
                    seq_number = list()
                    for residue in chain:
                        if residue.get_resname()[1] in ['A', 'C', 'G', 'T']:
                            seq.append(residue.get_resname()[1])
                        else:
                            seq.append("X")
                        seq_number.append((residue.get_id()[1], residue.get_resname()))
                elif len(chain.get_unpacked_list()[0].resname.strip()) == 3:  # Protein
                    seq = list()
                    seq_number = list()
                    counter = 0
                    for residue in chain:
                        if is_aa(residue.get_resname(), standard=True):  # aminoacids
                            seq.append(three_to_one(residue.get_resname()))
                            counter += 1
                        else:
                            seq.append("X")
                        seq_number.append((residue.get_id()[1], residue.get_resname()))
                    if counter == 0:  # HETAM
                        for residue in chain:
                            seq.append(residue.get_resname())
                while seq[-1] == "X":
                    del seq[-1]
                while seq[0] == "X":
                    del seq[0]
                chainDicLength[chain.id] = len(seq)
                chainDicFirstResidue[chain.id] = seq_number
            listOfChains[model.id] = chainDicLength
            listOfResidues[model.id] = chainDicFirstResidue

        return listOfChains, listOfResidues

    def getSequenceFromChain(self, modelID, chainID):
        self.checkRead()
        seq = list()
        for model in self.structure:
            if model.id == modelID:
                for chain in model:
                    if str(chain.id) == chainID:
                        if len(chain.get_unpacked_list()[0].resname) == 1:
                            print("Your sequence is a nucleotide sequence ("
                                  "RNA)\n")
                            # alphabet = IUPAC.IUPACAmbiguousRNA._upper()
                            for residue in chain:
                                # Check if the residue belongs to the
                                # standard RNA and add those residues to the
                                # seq
                                if residue.get_resname() in ['A', 'C',
                                                             'G', 'U']:
                                    seq.append(residue.get_resname())
                                else:
                                    seq.append("X")
                        elif len(chain.get_unpacked_list()[0].resname) == 2:
                            print("Your sequence is a nucleotide sequence ("
                                  "DNA)\n")
                            # alphabet = IUPAC.ExtendedIUPACDNA._upper()
                            for residue in chain:
                                # Check if the residue belongs to the
                                # standard DNA and add those residues to the
                                # seq
                                if residue.get_resname()[1] in ['A', 'C',
                                                                'G', 'T']:
                                    seq.append(residue.get_resname()[1])
                                else:
                                    seq.append("X")
                        elif len(chain.get_unpacked_list()[0].resname) == 3:
                            counter = 0
                            for residue in chain:
                                if is_aa(residue.get_resname(), standard=True):
                                    # alphabet = IUPAC.ExtendedIUPACProtein._upper()
                                    # The test checks if the amino acid
                                    # is one of the 20 standard amino acids
                                    # Some proteins have "UNK" or "XXX", or other symbols
                                    # for missing or unknown residues
                                    seq.append(three_to_one(residue.get_resname()))
                                    counter += 1
                                else:
                                    seq.append("X")
                            if counter != 0:  # aminoacids
                                print("Your sequence is an aminoacid sequence")
                            else:  # HETAM
                                print("Your sequence is a HETAM sequence")
                                for residue in chain:
                                    seq.append(residue.get_resname())
                        while seq[-1] == "X":
                            del seq[-1]
                        while seq[0] == "X":
                            del seq[0]
                        # return Seq(str(''.join(seq)), alphabet=alphabet)
                        return Seq(str(''.join(seq)))

    def getFullID(self, model_id='0', chain_id=None):
        """
        assign a label to a sequence obtained from a PDB file
        :parameter
        model_id
        chain_id
        :return: string with label
        """
        self.checkRead()  # cehck we have read the structure
        label = "%s" % self.structure.get_id()  # PDB_ID
        # check how many model are in sequence if > 1 add to id
        models = self.getModelsChains()
        if len(models) > 1:
            label = label + "__%s" % str(model_id)
        # if chain_id is None do not add it to te name
        if chain_id is not None:
            label = label + "_%s" % str(chain_id)
        return label

    def readLowLevel(self, fileName):
        """ Return a dictionary with all mmcif fields. you should parse them
            Example: get the list of the y coordinates of all atoms
              dict = readLowLevel("kk.pdb")
              y_list = dict['_atom_site.Cartn_y']
        """

        if fileName.endswith(".pdb"):
            print("Low level access to PDB is not implemented")
        else:
            dict = MMCIF2Dict(fileName)
        return dict

    def _write(self, fileName):
        """ Do not use this function use toPDB or toCIF, they take care of some
        compatibiity issues"""
        if fileName.endswith(".pdb") or fileName.endswith(".ent"):
            if self.ioPDB is None:
                self.ioPDB = PDBIO()
            io = self.ioPDB
        else:
            if self.ioCIF is None:
                self.ioCIF = scipionMMCIFIO()
                # self.ioCIF = MMCIFIO()
            io = self.ioCIF
        io.set_structure(self.structure)
        io.save(fileName)

    def _writeLowLevel(self, fileName, dict):
        """ write a dictionary as cif file
        """

        if fileName.endswith(".pdb"):
            print("Low level access to PDB is not implemented")
        else:
            if self.ioCIF is None:
                self.ioCIF = MMCIFIO()
            io = self.ioCIF
        io.set_dict(dict)
        io.save(fileName)

    def _intToChain(self, i, base=62):
        """
        int_to_chain(int,int) -> str
        Converts a positive integer to a chain ID. Chain IDs include uppercase
        characters, numbers, and optionally lowercase letters.
        i = a positive integer to convert
        base = the alphabet size to include. Typically 36 or 62.
        """
        if i < 0:
            raise ValueError("positive integers only")
        if base < 0 or 62 < base:
            raise ValueError("Invalid base")

        quot = int(i) // base
        rem = i % base
        if rem < 26:
            letter = chr(ord("A") + rem)
        elif rem < 36:
            letter = str(rem - 26)
        else:
            letter = chr(ord("a") + rem - 36)
        if quot == 0:
            return letter
        else:
            return self._intToChain(quot - 1, base) + letter

    def renameChains(self, structure):
        """Renames chains to be one-letter chains

        Existing one-letter chains will be kept.
        Multi-letter chains will be truncated
        or renamed to the next available letter of the alphabet.

        If more than 62 chains are present in the structure,
        raises an OutOfChainsError

        Returns a map between new and old chain IDs, as well as modifying
        the input structure

        """
        next_chain = 0  #
        # single-letters stay the same
        chainmap = {c.id: c.id
                    for c in structure.get_chains() if len(c.id) == 1}
        for o in structure.get_chains():
            if len(o.id) != 1:
                if o.id[0] not in chainmap:
                    chainmap[o.id[0]] = o.id
                    o.id = o.id[0]
                else:
                    c = self._intToChain(next_chain)
                    while c in chainmap:
                        next_chain += 1
                        c = self._intToChain(next_chain)
                        if next_chain >= 62:
                            raise OutOfChainsError()
                    chainmap[c] = o.id
                    o.id = c
        return chainmap

    def write(self, fileName):
        if fileName.endswith(".pdb") or fileName.endswith(".ent"):
            self.writeAsPdb(fileName)
        else:
            self.writeAsCif(fileName)

    def writeAsPdb(self, pdbFile):
        """ Save structure as PDB. Be aware that this is not a lossless conversion
        Returns False is conversion is not possible. True otherwise
        """
        # check input is not PDB
        if self.type == self.PDB:
            pass
        else:
            # rename long chains
            try:
                chainmap = self.renameChains(self.structure)
            except OutOfChainsError:
                print("Too many chains to represent in PDB format")
                return False

            for new, old in chainmap.items():
                # for new, old in list(chainmap.items()):
                if new != old:
                    print("Renaming chain {0} to {1}".format(old, new))

        self._write(pdbFile)

        return True

    def writeAsCif(self, cifFile):
        """ Save structure as CIF.
            Be aware that this is not a lossless conversion
        """
        self._write(cifFile)

    def getBoundingBox(self, expand=3):
        """Get bounding box (angstroms) for atom struct.
         Only alpha carbons are taken into account.
         parameter: expand.- make box larger addind this factor
         """
        # Use the first model
        ref_model = self.getStructure()[0]
        # init bounding box volues
        xmin = 100000; ymin = xmin; zmin = xmin
        xmax = -100000; ymax = xmax; zmax = xmax
        # iterate for all aminoacids
        for chain in ref_model:
            for residue in chain:
                (x,y,z) = residue["CA"].get_coord()
                if x < xmin: xmin = x
                if x > xmax: xmax = x

                if y < ymin: ymin = y
                if y > ymax: ymax = y

                if z < zmin: zmin = z
                if z > zmax: zmax = z
            # DEBUG
            # with open("/tmp/kk.bild", "w") as file:
            #     file.write(".transparency 0.8\n"
            #                ".color red\n"
            #                ".box %f %f %f %f %f %f\n" %
            #                (xmin - expand , ymin - expand, zmin - expand,
            #                xmax + expand, ymax + expand, zmax + expand))
        return [[xmin - expand, ymin - expand, zmin - expand],
                [xmax + expand, ymax + expand, zmax + expand]]


    def getTransformMatrix(self, atomStructFn, startId=-1, endId=-1):
        """find matrix that Superimposes two atom structures.
        this matrix moves atomStructFn to self """
        if endId == -1:
            endId = 10000000
        # load second atom structure
        aSH = AtomicStructHandler()
        aSH.read(atomStructFn)

        # Use the first model in the atom struct for alignment
        ref_model = self.getStructure()[0]
        sample_model = aSH.getStructure()[0]

        # Make a list of the atoms (in the structures) you wish to align.
        # In this case we use CA atoms whose index is in the
        # specified range (starId, endId)
        ref_atoms = []
        sample_atoms = []

        # Now get a list with CA atoms (alpha-carbon)
        # I assume that both atom structs have the same number of chains
        # and in the same order.
        for ref_chain, sample_chain in zip(ref_model, sample_model):
            # create a set with the id of all residues
            # for the current chain
            ref_set_id = set(res.get_id()[1] for res in ref_chain)
            sample_set_id = set(res.get_id()[1] for res in sample_chain)
            # keep the intersection as an ordered list (sets are not ordered)
            ref_set_id.intersection_update(sample_set_id)
            ref_list_id = sorted(ref_set_id)
            # delete AA smaller or large than the predefined values
            ref_list_id[:] = [x for x in ref_list_id
                              if x >= startId and x <= endId]

            # Iiterate through all residues in each chain a store CA atoms
            for id in ref_list_id:
                ref_atoms.append(ref_chain[id]['CA'])
                sample_atoms.append(sample_chain[id]['CA'])
        # Now we initiate the superimposer:
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_atoms, sample_atoms)
        # super_imposer.apply(sample_model.get_atoms())
        (rot, trans) = super_imposer.rotran
        # DEBUG, uncomment next two lines to see
        # transformation applied to atomStructFn
        # aSH.getStructure().transform(rot, trans)
        # aSH.write("/tmp/pp.cif")

        # convert 3x3 rotation matrix to homogeneous matrix
        rot = numpy.transpose(rot)  # scipion and biopython use
                                    # different conventions
        tmp = numpy.r_[rot, numpy.zeros((1, 3))]
        mat = numpy.c_[tmp, numpy.array(
            [[trans[0]], [trans[1]], [trans[2]], [1]])]
        return mat, super_imposer.rms

    def centerOfMass(self, geometric=False):
        """
        Returns gravity [default] or geometric center of mass of an Entity
        (anything with a get_atoms function in biopython.
        Geometric assumes all masses are equal (geometric=True)
        """
        entity = self.structure
        # Structure, Model, Chain, Residue
        if isinstance(entity, Entity.Entity):
            atom_list = entity.get_atoms()
        # List of Atoms
        elif hasattr(entity, '__iter__') and \
                [x for x in entity if x.level == 'A']:
            atom_list = entity
        else:  # Some other weirdo object
            raise ValueError("Center of Mass can only be calculated "
                             "from the following objects:\n"
                             "Structure, Model, Chain, Residue, "
                             "list of Atoms.")

        masses = []
        positions = [[], [], []]

        for atom in atom_list:
            masses.append(atom.mass)

            for i, coord in enumerate(atom.coord.tolist()):
                positions[i].append(coord)

        # If there is a single atom with undefined mass complain loudly.
        if 'ukn' in set(masses) and not geometric:
            raise ValueError("Some Atoms don't have an element assigned.\n"
                             "Try adding them manually or calculate the "
                             "geometrical center of mass instead.")

        if geometric:
            return [sum(coord_list) / len(masses) for coord_list in positions]
        else:
            w_pos = [[], [], []]
            for atom_index, atom_mass in enumerate(masses):
                w_pos[0].append(positions[0][atom_index] * atom_mass)
                w_pos[1].append(positions[1][atom_index] * atom_mass)
                w_pos[2].append(positions[2][atom_index] * atom_mass)

            return [sum(coord_list) / sum(masses) for coord_list in w_pos]

    def transform(self, transformation_matrix, sampling=1.):
        """ Geometrical transformation of a PDB structure

        :param entity: PDB biopython structure
        :param transformation_matrix -> 4x4 scipion matrix
        :param sampling: scipion transform matrix is applied to voxels so
              length must be multiplied by samplingRate

        internal variables:
             rotation matrix -> numpy.array(
                      [[ 0.5, -0.809017,  0.309017],
                       [ 0.809017,  0.30917,  -0.5     ],
                       [ 0.309017,  0.5,       0.809017]])
             translation: translation vector -> numpy.array([1., 0., 0.], 'd')
        :return: no return, new data overwrites entity
        """
        # bioPhython and Scipion conventions do not match
        rotation_matrix = numpy.transpose(transformation_matrix[:3, :3])
        # from geometry get euler angles and recreate matrix
        translation = translation_from_matrix(transformation_matrix)
        translation = [x * sampling for x in translation]
        self.structure.transform(rotation_matrix, translation)

    def _renameChainsIfNeed(self, struct2):
        """Rename chain, we assume that there is a single model per structure"""
        repeated = False

        def RepresentsInt(s):
            try:
                int(s)
                return True
            except ValueError:
                return False

        import uuid
        chainIDs1 = [chain.id for chain in self.structure.get_chains()]
        for chain in struct2.get_chains():
            counter = 2
            if chain.id in chainIDs1:
                repeated = True
                cId = chain.id
                l = len(cId)
                if l == 1:
                    while True:
                        try:
                            chain.id = "%s%03d" % (cId, counter)
                            break
                        except ValueError:
                            counter +=1
                            if counter > 1000:
                                raise ValueError('Error in _renameChainsIfNeed.')
                elif RepresentsInt(cId[1:]):  # try to fit a number and increase it by one
                    chain.id = "%s%03d" % (cId[0], int(cId[1:]) + 1)
                else:  # generate a 4 byte random string
                    chain.id = uuid.uuid4().hex[:4]
        if repeated:
            self._renameChainsIfNeed(struct2)

    def addStruct(self, secondPDBfileName, outPDBfileName=None, useModel=False):
        """ Join the second structure to the first one.
            If cheon numes are the same rename them.
            if outPDBfileName id provided then new
            struct is saved to a file"""
        # read new structure
        if outPDBfileName is not None:
            pdbID = (os.path.splitext(os.path.basename(outPDBfileName))[0])[:4]
        else:
            pdbID = (os.path.splitext(os.path.basename(secondPDBfileName))[0])[:4]

        if secondPDBfileName.endswith(".pdb") or secondPDBfileName.endswith(".ent"):
            parser = PDBParser(PERMISSIVE=self.permissive)
        else:
            parser = MMCIFParser()

        struct2 = parser.get_structure(pdbID, secondPDBfileName)

        if useModel:
            modelNumber = 0
            modelID = 0
            # model.id = model.serial_num = len(self.structure)?  # not sure this
            # is valid always
            for model in self.structure:
                pass
            modelNumber = model.serial_num
            modelID = model.id
            for model in struct2:
                modelNumber += 1
                modelID += 1
                model.detach_parent()
                model.serial_num = modelNumber
                model.id = modelID
                self.structure.add(model)
        else:
            self._renameChainsIfNeed(struct2)

            for model in struct2:
                for chain in model:
                    chain.detach_parent()
                    self.structure[0].add(chain)

        # create new output file
        if outPDBfileName is not None:
            self.write(outPDBfileName)

    def extractChain(self, chainID, start=0, end=-1,
                     modelID='0', filename="output.mmcif"):
        """Code for chopping  a structure.
        This module is used internally by the Bio.PDB.extract() function.
        """
        sel = ChainSelector(chain_id=chainID, start=start,
                            end=end, model_id=int(modelID))
        io = scipionMMCIFIO()
        io.set_structure(self.structure)
        io.save(filename, sel)

    def renameChain(self, chainID, newChainName, modelID='0',
                    filename="output.mmcif"):
        self.structure[modelID][chainID].id = newChainName
        self.write(filename)

    def renumberChain(self, chainID, offset=0, modelID='0',
                      filename="output.mmcif"):
        # get chain object
        chain = self.structure[modelID][chainID]
        # remove chain from model
        self.structure[modelID].detach_child(chainID)
        from Bio.PDB.Chain import Chain
        # create new chain
        newChain = Chain(chainID)
        for residue in chain:
            # remove residue, otherwise we cannot renumber it
            residue.detach_parent()
            rId = residue.id
            res_id = list(rId)
            res_id[1] = res_id[1] + offset
            if res_id[1] < 0:
                raise ValueError('Residue number cant be <= 0')
            residue.id = tuple(res_id)
            newChain.add(residue)
        self.structure[modelID].add(newChain)
        self.write(filename)


def cifToPdb(fnCif, fnPdb):
    h = AtomicStructHandler()
    h.read(fnCif)
    h.writeAsPdb(fnPdb)


def pdbToCif(fnPdb, fnCif):
    h = AtomicStructHandler()
    h.read(fnPdb)
    h.writeAsCif(fnCif)


def toPdb(inFileName, outPDBFile):
    if inFileName.endswith(".pdb") or inFileName.endswith(".ent"):
        return inFileName
    elif inFileName.endswith(".cif") or inFileName.endswith(".mmcif"):
        cifToPdb(inFileName, outPDBFile)
        return outPDBFile
    else:
        print("ERROR (toPdb), Unknown file type for file = %s" % inFileName)


def toCIF(inFileName, outCIFFile):
    if inFileName.endswith(".cif") or inFileName.endswith(".mmcif"):
        return inFileName
    elif inFileName.endswith(".pdb") or inFileName.endswith(".ent"):
        pdbToCif(inFileName, outCIFFile)
        return outCIFFile
    else:
        print("ERROR (toCIF), Unknown file type for file = %s" % inFileName)


def getEnviron():
    environ = pwutils.Environ(os.environ)
    environ.update({'RCSBROOT': os.path.join(Plugin.getMaxitHome())
                    }, position=pwutils.Environ.REPLACE)
    environ.update({'PATH': os.path.join(Plugin.getMaxitHome(), 'bin')
                    }, position=pwutils.Environ.BEGIN)
    return environ


def _frombase(inFileName, outFileName, log, oParam=1):
    # check if maxit exists,
    # if it does not then complain
    # convert pdb to cif using maxit
    global maxitAvailable
    try:
        maxitAvailable
    except:
        if not os.path.exists(Plugin.getMaxitBin()):
            maxitAvailable = False
            # show error message
        else:
            maxitAvailable = True

    if maxitAvailable:
        args = ' -input "' + inFileName + '" -output "' + outFileName + \
               '" -o %d' % oParam
        log.info('Launching: ' + Plugin.getMaxitBin() + args)
        # run in the background
        env = getEnviron()
        pwutils.runJob(None, Plugin.getMaxitBin(), args, env=env)
    else:
        # this is not the ideal conversion but it is better
        # than nothing
        aSH = AtomicStructHandler()
        aSH.read(inFileName)
        aSH.write(outFileName)
        # show error message
        print(pwutils.redStr("Please, install maxit with the command 'scipion installb maxit'"))
        print(pwutils.redStr("and restart scipion. Packages bison and flex are needed."))
        print(pwutils.redStr("If maxit is installed check %s in scipion.conf" % MAXIT_HOME))


def fromPDBToCIF(inFileName, outFileName, log):
    _frombase(inFileName, outFileName, log, 1)


def fromCIFToPDB(inFileName, outFileName, log):
    _frombase(inFileName, outFileName, log, 2)


def fromCIFTommCIF(inFileName, outFileName, log):
    _frombase(inFileName, outFileName, log, 8)


def retry(runEnvirom, program, args, cwd, listAtomStruct=[], log=None, clean_dir=None):
    try:
        runEnvirom(program, args, cwd=cwd)
    except:
        # first remove every files or directories in the extra folder,
        # except maps
        partiallyCleaningFolder(program, cwd)
        # something went wrong, may be bad atomStruct format
        log.info('retry with maxit conversion')

        for i, atomStructName in enumerate(listAtomStruct):
            if atomStructName.endswith(".pdb.cif"):
                aSH = AtomicStructHandler()
                aSH.read(atomStructName)
                aSH.write(atomStructName)

                runEnvirom(program, args, cwd=cwd)
            else:
                try:
                    if atomStructName.endswith(".pdb"):
                        newAtomStructName = os.path.abspath(
                            os.path.join(cwd, "retrypdb%d.cif" % i))
                        fromPDBToCIF(atomStructName, newAtomStructName, log)
                        _args = args.replace(atomStructName, newAtomStructName)
                        try:
                            runEnvirom(program, _args, cwd=cwd)
                        except:
                            fromCIFTommCIF(newAtomStructName, newAtomStructName, log)
                            runEnvirom(program, _args, cwd=cwd)
                    elif atomStructName.endswith(".cif"):
                        newAtomStructName = os.path.abspath(
                            os.path.join(cwd, "retrycif%d.cif" % i))
                        fromCIFTommCIF(atomStructName, newAtomStructName, log)
                        _args = args.replace(atomStructName, newAtomStructName)
                        try:
                            runEnvirom(program, _args, cwd=cwd)
                        except:
                            newAtomStructName = os.path.abspath(
                                os.path.join(cwd, "retrycif%d.pdb" % i))
                            fromCIFToPDB(atomStructName, newAtomStructName, log)
                            _args = args.replace(atomStructName, newAtomStructName)
                            runEnvirom(program, _args, cwd=cwd)
                except:
                    # first remove files or directories in the extra folder,
                    # except maps
                    partiallyCleaningFolder(program, cwd)
                    # biopython conversion
                    aSH = AtomicStructHandler()
                    if atomStructName.endswith(".pdb") or atomStructName.endswith(".ent"):
                        newAtomStructName = atomStructName.replace(".pdb", ".cif"). \
                            replace(".ent", ".cif")
                    else:
                        newAtomStructName = atomStructName
                    try:
                        aSH.read(atomStructName)
                        aSH.write(newAtomStructName)
                        _args = args.replace(atomStructName, newAtomStructName)
                        runEnvirom(program, _args, cwd=cwd)
                    except:
                        print("CIF file standardisation failed.")
                        # atomStructName = newAtomStructName


# TODO this should no be here, ROB
def partiallyCleaningFolder(program, cwd):
    l1 = os.listdir(cwd)
    l2 = ['molprobity.out', 'molprobity_probe.txt', 'molprobity_coot.py']
    l3 = ['validation_cryoem.pkl']
    l4 = ['placed_model.pdb', 'placed_model.cif']

    for item in l1:
        if (program.endswith('real_space_refine.py') and
            (item.endswith('geo') or
             item.endswith('real_space_refined.log') or
             item.endswith('real_space_refine.pdb') or
             item.endswith('real_space_refined.cif'))) or \
                (program.endswith('molprobity.py') and item in l2) or \
                (program.endswith('validation_cryoem.py') and item in l3) or \
                (program.endswith('emringer.py') and
                 (item.endswith("emringer.csv") or
                  item.endswith("emringer.pkl") or
                  item.endswith("emringer_plots") or
                  item.startswith('emringer_transfer'))) or \
                (program.endswith('dock_in_map.py') and item in l4):
            path1 = os.path.join(cwd, item)
            if os.path.exists(path1):
                if os.path.isfile(path1) or os.path.islink(path1):
                    os.remove(path1)
                elif os.path.isdir(path1):
                    shutil.rmtree(path1)
