#!/usr/bin/env python3

import loos
import loos.pyloos
import sys
import numpy
import argparse



header = "#" + " ".join(sys.argv)
print(header)



parser = argparse.ArgumentParser(description="Compute correlation matrix")
parser.add_argument('system_file',
                    help="System file")
parser.add_argument('selection_string',
                    help="Selection string for atoms to use")
parser.add_argument('output_filename',
                    help="File to write correlation matrix")
parser.add_argument('--traj',
                    help="One or more trajectory files",
                    nargs='+'
                    )
parser.add_argument('--skip',
                    help="Number of frames to skip for each traj",
                    type=int,
                    default=0)
parser.add_argument('--stride',
                    help="Step size for striding through trajectories",
                    type=int,
                    default=1)
args = parser.parse_args()

system = loos.createSystem(args.system_file)
subsystem = loos.selectAtoms(system, args.selection_string)
num_beads = len(subsystem)

trajs = []
for t in args.traj:
    traj = loos.pyloos.Trajectory(t,
                                  system,
                                  skip=args.skip,
                                  stride=args.stride)
    trajs.append(traj)

# This does the alignment for us, and manages stepping through all of the
# trajectories consistently
# The "*" unpacks the list of trajs into individual arguments
vtraj = loos.pyloos.AlignedVirtualTrajectory(*trajs,
                                             alignwith=args.selection_string)

# get the average structure from the aligned trajectory
average = loos.pyloos.ensembles.averageStructure(vtraj)
average = loos.selectAtoms(average, args.selection_string)

# create empty matrix to store the correlation matrix
correl = numpy.zeros([num_beads, num_beads], dtype=float)
distances = numpy.zeros([num_beads], dtype=float)

for frame in vtraj:
    # the way loos works, when the trajectory gets read, it gets written into
    # system. Since subsystem was created from system, it updates too.
    diffs = subsystem.differenceVectors(average)
    for i in range(len(diffs)):
        distances[i] = diffs[i].length()

    # there's got to be a fast numpyish way to do this
    for i in range(num_beads):
        for j in range(i, num_beads):
            correl[i, j] += distances[i]*distances[j]

# symmetrize the matrix
for i in range(num_beads):
    for j in range(i, num_beads):
        correl[j, i] += correl[i, j]

# normalize the matrix
correl /= len(vtraj)

numpy.savetxt(args.output_filename, correl)

# You can put the networkx stuff here, or in a separate script that
# reads in the correlation matrix using
# correl = numpy.loadtxt("filename")
# The call will look something like
# graph = networkx.from_numpy_array(correl)
