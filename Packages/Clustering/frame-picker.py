#!/usr/bin/env python3
"""
frame_picker.py: make trajectories out of the clusters.

The purpose of this tool is to read the output of cluster-kgs and create a
set of trajectory files, each containing 1 cluster.
"""

import loos
from loos import pyloos
import json
from sys import argv

# usage = 'python frame-picker.py results.json start:step:stop prefix model traj [traj] ...'
usage = """
    frame-picker.py results.json start:step: prefix model traj [traj] ...

        results.json:    the output from cluster-kgs
        start:step:stop: how to step through the trajectories. The same options
                         given to rmsds or multi-rmsds should be used here
                         Note: at the moment, stop should not be specified
        prefix:          the core of the name of the output trajectories, eg.
                         ../output/foo would produce
                         ../output/foo-frac-cluster.dcd where frac is the
                         fraction of population in this cluster
                         and cluster is the index of that cluster. It will also
                         create ../output/foo-frac-cluster.pdb, where that
                         structure is the center of the cluster
        model:           model file used in the distance calculation
        traj:            one or more trajectories. Should be the same
                         trajectories you fed to rmsds or multi-rmsds.
"""

if len(argv) == 1 or "-h" in argv:
    print(usage)
    exit(0)

for i, arg in enumerate(argv):
    print(i, arg)

cluster_results_fn = argv[1]
# parse clustering output here
with open(cluster_results_fn, "r") as fp:
    cluster_results = json.load(fp)

# trajectory stuff
prerange = argv[2].split(":")
if len(prerange) != 3:
    print("Stop feeding me weird things!")
    print(trajranges)
    print(usage)
    exit(0)

if prerange[-1]:
    trajranges = tuple(map(int, prerange))
    print('Warning: code does not currently handle "stop"')

else:
    print('Warning: code does not currently handle "stop"')
    trajranges = (int(prerange[0]), int(prerange[1]), -1)
skip, stride, stop = trajranges
prefix = argv[3]
model = loos.createSystem(argv[4])
traj = pyloos.Trajectory(argv[5], model, skip=skip, stride=stride)
vtraj = pyloos.VirtualTrajectory(traj)

# if there are more trajectories, add them to the vtraj
if len(argv) > 6:
    for tfile in argv[6:]:
        vtraj.append(pyloos.Trajectory(tfile, model, skip=skip, stride=stride))


total = len(vtraj)
print("Frames to process: ", total)

for cluster_ix, cluster in enumerate(cluster_results["clusters"]):
    # puts the fraction of the trajectory lumped into ith traj in the name
    frac = round(len(cluster) / float(total), 4)
    outprefix = prefix + "-" + str(frac) + "-" + str(cluster_ix)
    cluster_traj = loos.DCDWriter(outprefix + ".dcd")
    exemplar_ix = cluster_results["exemplars"][cluster_ix]
    frame = vtraj[exemplar_ix]
    exemplar = loos.PDB.fromAtomicGroup(frame)
    with open(outprefix + ".pdb", "w") as fp:
        fp.write(str(exemplar))

    # OK now write out the trajectory.
    for member_idx in cluster:
        member = vtraj[member_idx]
        cluster_traj.writeFrame(member)
