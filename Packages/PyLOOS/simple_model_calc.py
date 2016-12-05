#!/usr/bin/env python
"""
simple_model_calc.py : simple skeleton of a program that reads in a 
    structure, makes a selection, then loops over the atoms.

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


