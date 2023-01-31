# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
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
from urllib.error import URLError, HTTPError
from urllib.request import urlopen
from os.path import exists


# sequence related stuff


from Bio.Seq import Seq
from Bio import Entrez, SeqIO
import sys
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline
from Bio import pairwise2
from pyworkflow.utils import getExt
from pwem.objects.data import Alphabet


class SequenceHandler:
    def __init__(self, sequence=None,
                 iUPACAlphabet=Alphabet.DUMMY_ALPHABET, doClean=True):
        iUPACAlphabet = Alphabet.DUMMY_ALPHABET if iUPACAlphabet not in Alphabet.alphabets else iUPACAlphabet

        if iUPACAlphabet == Alphabet.DUMMY_ALPHABET:
            self.isAminoacid = None
        else:      
            self.isAminoacid =  (iUPACAlphabet <   Alphabet.AMBIGOUS_DNA_ALPHABET)

        self.alphabet = iUPACAlphabet ##  indexToAlphabet(isAminoacid, iUPACAlphabet)

        if sequence is not None:
            if doClean:
                self._sequence = cleanSequence(self.alphabet, sequence)
            else:
                self._sequence = sequence
        else:
            self._sequence = None

    def getSequence(self):
        """Returns the sequence as Bio.Seq.Seq object"""
        return Seq(self._sequence)

    def getTypeFromFile(self, fileName):
        '''Returns the expected BioPython file type according to the filename'''
        ext = getExt(fileName)
        if ext == '.fasta':
            type = 'fasta'
        elif ext == '.genbank' or ext == '.gb':
            type = 'genbank'
        else:
            type = 'fasta'
        return type

    def appendFile(self, fileName, seqID, sequence=None, name=None,
                 seqDescription=None, type="fasta"):
        '''Appends a sequence to the sequences in a file'''
        if sequence is not None:
            self._sequence = sequence

        if type is None:
            type = self.getTypeFromFile(fileName)

        if exists(fileName):
            records = list(SeqIO.parse(fileName, type))
        else:
            records = []
        records.append(SeqRecord(Seq(self._sequence), id=seqID, name=name,
                                 description=seqDescription))
        with open(fileName, "w") as output_handle:
            SeqIO.write(records, output_handle, type)

    def saveFile(self, fileName, seqID, sequence=None, name=None,
                 seqDescription=None, type="fasta"):
        '''Saves a single sequence into a specified file'''
        if sequence is not None:
            self._sequence = sequence

        self.saveFiles(fileName, [seqID], [self._sequence], [name], [seqDescription], type)

    def saveFiles(self, fileName, seqIDs, sequences, names,
                  seqDescriptions, type="fasta"):
        '''Saves multiple sequences into a specified file'''
        if type is None:
            type = self.getTypeFromFile(fileName)

        records = []
        for seq, seqID, name, seqDescription in zip(sequences, seqIDs, names, seqDescriptions):
            records.append(SeqRecord(Seq(seq), id=seqID, name=name,
                               description=seqDescription))
        # type(record): < class 'Bio.SeqRecord.SeqRecord'>
        with open(fileName, "w") as output_handle:
            SeqIO.write(records, output_handle, type)

    def readSequenceFromFile(self, fileName, type="fasta", isAmino=True):
        '''From a sequences file, returns a dictionary with ther FIRST sequence info.
        Dictionary: {'seqID': seqID1, 'sequence': sequence1, 'description': description1, 'alphabet': alphabet1}'''

        if type is None:
            type = self.getTypeFromFile(fileName)
        return self.readSequencesFromFile(fileName, isAmino=isAmino)[0]

    def readSequencesFromFile(self, fileName, type='fasta', isAmino=True):
        '''From a sequences file, returns a list of dictionaries with each sequence info.
        Dictionary: [{'seqID': seqID1, 'sequence': sequence1, 'description': description1, 'alphabet': alphabet1},
                     ...]'''

        if type is None:
            type = self.getTypeFromFile(fileName)

        sequences = []

        records = SeqIO.parse(fileName, type)

        for rec in records:
            sequence = str(rec.seq)

            if isAmino:
                self.alphabet = Alphabet.EXTENDED_PROTEIN_ALPHABET
            elif 'U' in sequence:
                self.alphabet = Alphabet.UNAMBIGOUS_RNA_ALPHABET
            elif 'T' in sequence:
                self.alphabet = Alphabet.UNAMBIGOUS_DNA_ALPHABET
            else:
                self.alphabet = Alphabet.NUCLEOTIDES_ALPHABET

            sequences.append({'seqID': rec.id, 'sequence': sequence,
                              'name': rec.name, 'description': rec.description,
                              'isAminoacids': isAmino, 'alphabet': self.alphabet})
        return sequences

    def downloadSeqFromDatabase(self, seqID, dataBase=None):
        # see http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
        # for format/databases

        seqID = str(seqID)
        sys.stdout.flush()
        counter = 1
        retries = 5
        record = None
        error = ""
        if dataBase is None:
            if self.isAminoacid:
                dataBase = 'UniProt'
            else:
                dataBase = 'GeneBank'

        while counter <= retries:  # retry up to 5 times if server busy
            try:
                if dataBase == 'UniProt':
                    url = "http://www.uniprot.org/uniprot/%s.xml"
                    format = "uniprot-xml"

                    handle = urlopen(url % seqID)
                    alphabet = Alphabet.EXTENDED_PROTEIN_ALPHABET
                else:
                    if self.isAminoacid:
                        db = "protein"
                        alphabet = Alphabet.EXTENDED_PROTEIN_ALPHABET
                    else:
                        db = "nucleotide"
                        alphabet = Alphabet.NUCLEOTIDES_ALPHABET
                    Entrez.email = "scipion@cnb.csic.es"
                    format = "fasta"
                    handle = Entrez.efetch(db=db, id=seqID,
                                           rettype=format, retmode="text")

                record = SeqIO.read(handle, format)
                break
            except HTTPError as e:
                error = "%s is a wrong sequence ID" % seqID
                print(e.code)
            except URLError as e:
                error = "Cannot connect to %s" % dataBase
                print(e.args)
            except Exception as ex:
                template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                message = template.format(type(ex).__name__, ex.args)
                error = message
            if counter == retries:
                break
            counter += 1

        seqDic = None
        if record is not None:
            seqDic = {'seqID': record.id, 'sequence': str(record.seq), 'description': record.description,
                      'alphabet': alphabet}
        return seqDic, error

    def alignSeq(self, referenceSeq):
        if self._sequence is not None:
            alignments = pairwise2.align.globalds(self._sequence.seq,
                                                  referenceSeq.seq)
            return alignments
        else:
            print("read the sequence first")
            exit(0)


def sequenceLength(filename, format='fasta'):
    handler = SequenceHandler()
    return len(handler.readSequenceFromFile(filename, format)['sequence'])


def cleanSequenceScipion(isAminoacid, iUPACAlphabet, sequence):
    return cleanSequence(iUPACAlphabet, sequence)


def cleanSequence(alphabet, sequence):
    """Remove all characters that are not in the alphabet
       :param alphabet: the alphabet to use as integer
       :param sequence: the sequence to clean
       """
    str_list = []
    for item in sequence.upper():
        if item in Alphabet.alphabets[alphabet]:
            str_list.append(item)
    value = ''.join(str_list)
    return ''.join(str_list)




def saveFileSequencesToAlign(SeqDic, inFile, type="fasta"):
    # Write my sequences to a fasta file
    with open(inFile, "w") as output_handle:
        for index, seq in SeqDic.items():
            record = SeqRecord(seq, id=str(index),
                               name="", description="")
            SeqIO.write(record, output_handle, type)


def alignClustalSequences(inFile, outFile):
    # Alignment of sequences with Clustal Omega program
    clustalomega_cline = ClustalOmegaCommandline(
            infile=inFile,
            outfile=outFile,
            verbose=True, auto=True)
    return clustalomega_cline


def alignMuscleSequences(inFile, outFile):
    # Alignment of sequences with Muscle program
    muscle_cline = MuscleCommandline(input=inFile, out=outFile)
    return muscle_cline


def alignBioPairwise2Sequences(structureSequenceId, structureSequence,
                               userSequenceId, userSequence,
                               outFileName):
    """aligns two sequences and saves them to disk using fasta format"""
    # see alignment_function for globalms parameters
    alignments = pairwise2.align.globalms(structureSequence,
                                          userSequence,  3, -1, -3, -2)
    align1, align2, score, begin, end = alignments[0]
    with open(outFileName, "w") as handle:
        handle.write(">%s\n%s\n>%s\n%s\n" % (structureSequenceId,
                                             align1,
                                             userSequenceId,
                                             align2))
