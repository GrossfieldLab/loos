#!/usr/bin/env python
"""
simple_traj_calc.py : simple skeleton of a program that reads in a 
    structure and a rajectory, makes a selection, then loops over the frame.

Alan Grossfield
University of Rochester Medical Center
"""

import loos
import loos.pyloos
import sys

header = " ".join(sys.argv)
print("#", header)

# parse the command line arguments -- in a more complex example,
# you'd use the argparse module
model_filename = sys.argv[1]
trajectory_filename = sys.argv[2]
selection_string = sys.argv[3]

# Create the system 
model = loos.createSystem(model_filename)

# Select a subset of the system
subset = loos.selectAtoms(model, selection_string)

# Create the trajectory
# You can automatically set skip and stride options here,
# so that when you loop over traj using a standard iterator
# those options are automatically applied
traj = loos.pyloos.Trajectory(trajectory_filename, model)


# Iterate over 
for frame in traj:
    # Example: compute the centroid of the selection
    centroid = subset.centroid()


# Write something out

