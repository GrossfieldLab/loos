#!/usr/bin/env python
"""
simple_traj_transform.py : simple skeleton of a program that reads in a 
    a structure and trajectory, selects a subset, does something, and
    writes out a new trajectory

Alan Grossfield
University of Rochester Medical Center
"""

import loos
import loos.pyloos
import sys

header = " ".join(sys.argv)
print "#", header

# parse the command line arguments -- in a more complex example,
# you'd use the argparse module
model_filename = sys.argv[1]
trajectory_filename = sys.argv[2]
selection_string = sys.argv[3]
output_trajectory_name = sys.argv[4]

# Create the system 
model = loos.createSystem(model_filename)

# Select a subset of the system
subset = loos.selectAtoms(model, selection_string)

# Create the trajectory
# You can automatically set skip and stride options here,
# so that when you loop over traj using a standard iterator
# those options are automatically applied
traj = loos.pyloos.Trajectory(trajectory_filename, model)

# Create the output trajectory
# We also have a DCDWriter, which (as the name suggests)
# writes DCD files
outtraj = loos.XTCWriter(output_trajectory_name)


# Iterate over trajectory
first_frame = True
for frame in traj:
    # Example of a trivial transformation
    subset.centerAtOrigin()

    # write out just the selection to the trajectory
    outtraj.writeFrame(subset)

    # if it's the first frame, write a pdb file too
    if first_frame:
        first_frame = False

        # Renumber the so that the atoms in the selection are
        # are sequential (using a copy so we don't screw up the 
        # original. 
        frame_copy = subset.copy()
        frame_copy.pruneBonds()
        frame_copy.renumber()

        # make a PDB object
        pdb = loos.PDB.fromAtomicGroup(frame_copy)
        pdb.remarks().add(header)
        print pdb


