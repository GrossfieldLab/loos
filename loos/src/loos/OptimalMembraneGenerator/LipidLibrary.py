#!/usr/bin/env python3

import os
import glob
import loos
import random
import sys

# @cond TOOLS_INTERNAL


class LipidLibrary:
    """
    Class to access a lipid library directory containing either PDB or
    CHARMM CRD files.
    """
    def __init__(self, path):
        self.path = path
        if not os.path.exists(path):
            print(self.error_message())
            sys.exit()

        all_files = os.listdir(self.path)

        # select files, either pdb or crd
        self.structures = []
        for file in all_files:
            if (file.endswith(".pdb") or file.endswith(".crd")):
                self.structures.append(file)

        if self.size() == 0:
            print(self.error_message())
            sys.exit()

    def size(self):
        return len(self.structures)

    def pick_structure(self):
        # pick a structure file at random, and build up the filename
        filename = os.path.join(self.path, random.choice(self.structures))

        # read the file
        lipid = loos.createSystem(filename)

        # return the AtomicGroup
        return lipid

    def __getitem__(self, index):
        return self.structures[index]

    def error_message(self):
        msg =  "Error reading lipid library (" + self.path + ")\n"
        msg += "Either directory doesn't exist or there are no pdb/crd files\n"
        msg += "Try downloading our lipid library from \n"
        msg += "http://sourceforge.net/projects/loos/files/loos%20material/lipid_library.tgz/download\n"
        msg += "Exiting...."

        return msg

# @endcond
