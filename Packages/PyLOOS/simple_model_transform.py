#!/usr/bin/env python
"""
simple_model_transform.py : simple skeleton of a program that reads in a 
    structure, makes a selection, performs some kind of transformation,
    then prints out pdb file

Alan Grossfield
University of Rochester Medical Center
"""

import loos
import sys

header = " ".join(sys.argv)
print "#", header

# parse the command line arguments -- in a more complex example,
# you'd use the argparse module
model_filename = sys.argv[1]
selection_string = sys.argv[2]

# Create the system 
model = loos.createSystem(model_filename)

# Select a subset of the system
subset = loos.selectAtoms(model, selection_string)


# Iterate over atoms and calculate something#
for atom in subset:
    # do something
    continue

# Convert the subset to a PDB file
pdb = loos.PDB.fromAtomicGroup(subset)

# add the command line to the pdb file's header
pdb.remarks().add(header)

print pdb


