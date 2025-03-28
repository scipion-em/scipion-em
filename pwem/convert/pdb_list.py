#!/usr/bin/env python
# Copyright 2003, by Kristian Rother. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# PDBList.py
#
# A tool for tracking changes in the PDB Protein Structure Database.
#
# (c) 2003 Kristian Rother
# This work was supported by the German Ministry of Education
# and Research (BMBF). Project http://www.bcbio.de
#
# Contact the author
#    homepage : http://www.rubor.de/bioinf
#    email    : krother@genesilico.pl
#
#
# (c) 2016 Wiktoria Karwicka & Jacek Smietanski
#   - updated and Python 3.x compatible code
#   - new options to enable download PDBx/mmCif, PDBML and mmtf formatted
#       files as well as large PDB bundles
#   - unit tests for the module
#
# Contact the corresponding author
#   homepage : http://jaceksmietanski.net
#   email    : jacek.smietanski@ii.uj.edu.pl
#
# It may be distributed freely with respect to the original authors.
# Any maintainer of the Biopython code may change this notice
# when appropriate.

"""Access the PDB over the internet (e.g. to download structures)."""


import contextlib
import ftplib
import gzip
import os
import re
import shutil
import sys

from urllib.request import urlopen
from urllib.request import urlretrieve
from urllib.request import urlcleanup


class PDBList:
    """Quick access to the structure lists on the PDB or its mirrors.

    This class provides quick access to the structure lists on the
    PDB server or its mirrors. The structure lists contain
    four-letter PDB codes, indicating that structures are
    new, have been modified or are obsolete. The lists are released
    on a weekly basis.

    It also provides a function to retrieve PDB files from the server.
    To use it properly, prepare a directory /pdb or the like,
    where PDB files are stored.

    All available file formats (PDB, PDBx/mmCif, PDBML, mmtf) are supported.
    Please note that large structures (containing >62 chains
    and/or 99999 ATOM lines) are no longer stored as a single PDB file
    and by default (when PDB format selected) are not downloaded.

    Large structures can be downloaded in other formats, including PDBx/mmCif
    or as a .tar file (a collection of PDB-like formatted files for a given
    structure).

    If you want to use this module from inside a proxy, add
    the proxy variable to your environment, e.g. in Unix:
    export HTTP_PROXY='http://realproxy.charite.de:888'
    (This can also be added to ~/.bashrc)
    """

    PDB_REF = """
    The Protein Data Bank: a computer-based archival file for macromolecular structures.
    F.C.Bernstein, T.F.Koetzle, G.J.B.Williams, E.F.Meyer Jr, M.D.Brice, J.R.Rodgers, O.Kennard, T.Shimanouchi, M.Tasumi
    J. Mol. Biol. 112 pp. 535-542 (1977)
    http://www.pdb.org/.
    """

    # Code line modified from the original 
    # biopython code so now https is used instead of ftp
    # server="ftp://ftp.wwpdb.org" -> server="https://files.wwpdb.org"
    def __init__(
        self, server="https://files.wwpdb.org", 
        pdb=None, obsolete_pdb=None, verbose=True
    ):
        """Initialize the class with the default server or a custom one.

        Argument pdb is the local path to use, defaulting to the current
        directory at the moment of initialisation.
        """
        self.pdb_server = server  # remote pdb server
        if pdb:
            self.local_pdb = pdb  # local pdb file tree
        else:
            self.local_pdb = os.getcwd()

        # enable or disable verbose
        self._verbose = verbose

        # local file tree for obsolete pdb files
        if obsolete_pdb:
            self.obsolete_pdb = obsolete_pdb
        else:
            self.obsolete_pdb = os.path.join(self.local_pdb, "obsolete")
            if not os.access(self.obsolete_pdb, os.F_OK):
                os.makedirs(self.obsolete_pdb)

        # variable for command-line option
        self.flat_tree = False

    @staticmethod
    def _print_default_format_warning(file_format):
        """Print a warning to stdout (PRIVATE).

        Temporary warning (similar to a deprecation warning) that files
        are being downloaded in mmCIF.
        """
        if file_format is None:
            sys.stderr.write(
                "WARNING: The default download format has changed from PDB to PDBx/mmCif\n"
            )
            return "mmCif"
        return file_format

    @staticmethod
    def get_status_list(url):
        """Retrieve a list of pdb codes in the weekly pdb status file from given URL.

        Used by get_recent_changes. Typical contents of the list files parsed
        by this method is now very simply - one PDB name per line.
        """
        with contextlib.closing(urlopen(url)) as handle:
            answer = []
            for line in handle:
                pdb = line.strip()
                assert len(pdb) == 4
                answer.append(pdb.decode())
        return answer

    def get_recent_changes(self):
        """Return three lists of the newest weekly files (added,mod,obsolete).

        Reads the directories with changed entries from the PDB server and
        returns a tuple of three URL's to the files of new, modified and
        obsolete entries from the most recent list. The directory with the
        largest numerical name is used.
        Returns None if something goes wrong.

        Contents of the data/status dir (20031013 would be used);:

            drwxrwxr-x   2 1002     sysadmin     512 Oct  6 18:28 20031006
            drwxrwxr-x   2 1002     sysadmin     512 Oct 14 02:14 20031013
            -rw-r--r--   1 1002     sysadmin    1327 Mar 12  2001 README

        """
        path = self.pdb_server + "/pub/pdb/data/status/latest/"

        # Retrieve the lists
        added = self.get_status_list(path + "added.pdb")
        modified = self.get_status_list(path + "modified.pdb")
        obsolete = self.get_status_list(path + "obsolete.pdb")
        return [added, modified, obsolete]

    def get_all_entries(self):
        """Retrieve the big file containing all the PDB entries and some annotation.

        Returns a list of PDB codes in the index file.
        """
        url = self.pdb_server + "/pub/pdb/derived_data/index/entries.idx"
        if self._verbose:
            print("Retrieving index file. Takes about 27 MB.")
        with contextlib.closing(urlopen(url)) as handle:
            all_entries = [
                line[:4].decode() for line in handle.readlines()[2:] if len(line) > 4
            ]
        return all_entries

    def get_all_obsolete(self):
        """Return a list of all obsolete entries ever in the PDB.

        Returns a list of all obsolete pdb codes that have ever been
        in the PDB.

        Gets and parses the file from the PDB server in the format
        (the first pdb_code column is the one used). The file looks
        like this::

             LIST OF OBSOLETE COORDINATE ENTRIES AND SUCCESSORS
            OBSLTE    31-JUL-94 116L     216L
            ...
            OBSLTE    29-JAN-96 1HFT     2HFT
            OBSLTE    21-SEP-06 1HFV     2J5X
            OBSLTE    21-NOV-03 1HG6
            OBSLTE    18-JUL-84 1HHB     2HHB 3HHB
            OBSLTE    08-NOV-96 1HID     2HID
            OBSLTE    01-APR-97 1HIU     2HIU
            OBSLTE    14-JAN-04 1HKE     1UUZ
            ...

        """
        url = self.pdb_server + "/pub/pdb/data/status/obsolete.dat"
        with contextlib.closing(urlopen(url)) as handle:
            # Extract pdb codes. Could use a list comprehension, but I want
            # to include an assert to check for mis-reading the data.
            obsolete = []
            for line in handle:
                if not line.startswith(b"OBSLTE "):
                    continue
                pdb = line.split()[2]
                assert len(pdb) == 4
                obsolete.append(pdb.decode())
        return obsolete

    def retrieve_pdb_file(
        self, pdb_code, obsolete=False, pdir=None, file_format=None, overwrite=False
    ):
        """Fetch PDB structure file from PDB server, and store it locally.

        The PDB structure's file name is returned as a single string.
        If obsolete ``==`` True, the file will be saved in a special file tree.

        NOTE. The default download format has changed from PDB to PDBx/mmCif

        :param pdb_code: 4-symbols structure Id from PDB (e.g. 3J92).
        :type pdb_code: string

        :param file_format:
            File format. Available options:

            * "mmCif" (default, PDBx/mmCif file),
            * "pdb" (format PDB),
            * "xml" (PDBML/XML format),
            * "mmtf" (highly compressed),
            * "bundle" (PDB formatted archive for large structure)

        :type file_format: string

        :param overwrite: if set to True, existing structure files will be overwritten. Default: False
        :type overwrite: bool

        :param obsolete:
            Has a meaning only for obsolete structures. If True, download the obsolete structure
            to 'obsolete' folder, otherwise download won't be performed.
            This option doesn't work for mmtf format as obsoleted structures aren't stored in mmtf.
            Also doesn't have meaning when parameter pdir is specified.
            Note: make sure that you are about to download the really obsolete structure.
            Trying to download non-obsolete structure into obsolete folder will not work
            and you face the "structure doesn't exists" error.
            Default: False

        :type obsolete: bool

        :param pdir: put the file in this directory (default: create a PDB-style directory tree)
        :type pdir: string

        :return: filename
        :rtype: string
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)

        # Get the compressed PDB structure
        pdb_code = pdb_code.lower()
        archive = {
            "pdb": f"pdb{pdb_code}.ent.gz",
            "mmCif": f"{pdb_code}.cif.gz",
            "xml": f"{pdb_code}.xml.gz",
            "mmtf": f"{pdb_code}",
            "bundle": f"{pdb_code}-pdb-bundle.tar.gz",
        }
        archive_fn = archive[file_format]

        if file_format not in archive.keys():
            raise (
                f"Specified file_format {file_format} doesn't exists or is not supported. Maybe a "
                "typo. Please, use one of the following: mmCif, pdb, xml, mmtf, bundle"
            )

        if file_format in ("pdb", "mmCif", "xml"):
            pdb_dir = "divided" if not obsolete else "obsolete"
            file_type = (
                "pdb"
                if file_format == "pdb"
                else "mmCIF"
                if file_format == "mmCif"
                else "XML"
            )
            url = (
                self.pdb_server
                + f"/pub/pdb/data/structures/{pdb_dir}/{file_type}/{pdb_code[1:3]}/{archive_fn}"
            )
        elif file_format == "bundle":
            url = (
                self.pdb_server
                + f"/pub/pdb/compatible/pdb_bundle/{pdb_code[1:3]}/{pdb_code}/{archive_fn}"
            )
        else:
            url = f"https://mmtf.rcsb.org/v1.0/full/{pdb_code}"

        # Where does the final PDB file get saved?
        if pdir is None:
            path = self.local_pdb if not obsolete else self.obsolete_pdb
            if not self.flat_tree:  # Put in PDB-style directory tree
                path = os.path.join(path, pdb_code[1:3])
        else:  # Put in specified directory
            path = pdir
        if not os.access(path, os.F_OK):
            os.makedirs(path)
        filename = os.path.join(path, archive_fn)
        final = {
            "pdb": f"pdb{pdb_code}.ent",
            "mmCif": f"{pdb_code}.cif",
            "xml": f"{pdb_code}.xml",
            "mmtf": f"{pdb_code}.mmtf",
            "bundle": f"{pdb_code}-pdb-bundle.tar",
        }
        final_file = os.path.join(path, final[file_format])

        # Skip download if the file already exists
        if not overwrite:
            if os.path.exists(final_file):
                if self._verbose:
                    print(f"Structure exists: '{final_file}' ")
                return final_file

        # Retrieve the file(s)
        if self._verbose:
            print(f"Downloading PDB structure '{pdb_code}'...")
        try:
            urlcleanup()
            urlretrieve(url, filename)
        except OSError:
            print("Desired structure doesn't exist")
        else:
            with gzip.open(filename, "rb") as gz:
                with open(final_file, "wb") as out:
                    out.writelines(gz)
            os.remove(filename)
        return final_file

    def update_pdb(self, file_format=None, with_assemblies=False):
        """Update your local copy of the PDB files.

        I guess this is the 'most wanted' function from this module.
        It gets the weekly lists of new and modified pdb entries and
        automatically downloads the according PDB files.
        You can call this module as a weekly cron job.
        """
        assert os.path.isdir(self.local_pdb)
        assert os.path.isdir(self.obsolete_pdb)

        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)

        new, modified, obsolete = self.get_recent_changes()

        for pdb_code in new + modified:
            try:
                self.retrieve_pdb_file(pdb_code, file_format=file_format)
                if with_assemblies:
                    assemblies = self.get_all_assemblies(file_format)
                    for a_pdb_code, assembly_num in assemblies:
                        if a_pdb_code == pdb_code:
                            pl.retrieve_assembly_file(
                                pdb_code,
                                assembly_num,
                                file_format=file_format,
                                overwrite=True,
                            )

            except Exception as err:
                print(f"error {pdb_code}: {err}\n")
                # you can insert here some more log notes that
                # something has gone wrong.

        # Move the obsolete files to a special folder
        # NOTE: This should be updated to handle multiple file types and
        # assemblies. As of now, it only looks for PDB-formatted files.
        # Using pathlib will be probably the best approach here, to build
        # and index of which files we have efficiently (or glob them).
        for pdb_code in obsolete:
            if self.flat_tree:
                old_file = os.path.join(self.local_pdb, f"pdb{pdb_code}.ent")
                new_dir = self.obsolete_pdb
            else:
                old_file = os.path.join(
                    self.local_pdb, pdb_code[1:3], f"pdb{pdb_code}.ent"
                )
                new_dir = os.path.join(self.obsolete_pdb, pdb_code[1:3])
            new_file = os.path.join(new_dir, f"pdb{pdb_code}.ent")
            if os.path.isfile(old_file):
                if not os.path.isdir(new_dir):
                    os.mkdir(new_dir)
                try:
                    shutil.move(old_file, new_file)
                except Exception:
                    print(f"Could not move {old_file} to obsolete folder")
            elif os.path.isfile(new_file):
                if self._verbose:
                    print(f"Obsolete file {old_file} already moved")
            else:
                if self._verbose:
                    print(f"Obsolete file {old_file} is missing")

    def download_pdb_files(
        self, pdb_codes, obsolete=False, pdir=None, file_format=None, overwrite=False
    ):
        """Fetch set of PDB structure files from the PDB server and stores them locally.

        The PDB structure's file name is returned as a single string.
        If obsolete ``==`` True, the files will be saved in a special file tree.

        :param pdb_codes: a list of 4-symbols structure Ids from PDB
        :type pdb_codes: list of strings

        :param file_format:
            File format. Available options:

            * "mmCif" (default, PDBx/mmCif file),
            * "pdb" (format PDB),
            * "xml" (PMDML/XML format),
            * "mmtf" (highly compressed),
            * "bundle" (PDB formatted archive for large structure)

        :param overwrite: if set to True, existing structure files will be overwritten. Default: False
        :type overwrite: bool

        :param obsolete:
            Has a meaning only for obsolete structures.
            If True, download the obsolete structure
            to 'obsolete' folder, otherwise download won't be performed.
            This option doesn't work for mmtf format as obsoleted structures are not available as mmtf.
            (default: False)

        :type obsolete: bool

        :param pdir: put the file in this directory (default: create a PDB-style directory tree)
        :type pdir: string

        :return: filenames
        :rtype: string
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)
        for pdb_code in pdb_codes:
            self.retrieve_pdb_file(
                pdb_code,
                obsolete=obsolete,
                pdir=pdir,
                file_format=file_format,
                overwrite=overwrite,
            )
