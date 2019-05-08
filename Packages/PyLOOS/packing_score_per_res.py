#!/usr/bin/env python3

import sys
import loos
import loos.pyloos
import numpy


header = "#" + " ".join(sys.argv)

system_file = sys.argv[1]
traj_file = sys.argv[2]
protein_selection = sys.argv[3]
probe_selection = sys.argv[4]
output_filename_core = sys.argv[5]


system = loos.createSystem(system_file)
# TODO: I'll want to support skip and stride, and probably virtual trajectories
traj = loos.pyloos.Trajectory(traj_file, system)

protein = loos.selectAtoms(system, protein_selection)
probe = loos.selectAtoms(system, probe_selection)

residues = protein.splitByResidue()
probes = probe.splitByMolecule()

# pre-allocate storage as a numpy array
scores = numpy.zeros([len(residues), len(probes), len(traj)], numpy.float)

frame_index = 0
for frame in traj:
    box = frame.periodicBox()
    for r in range(len(residues)):
        for p in range(len(probes)):
            s = residues[r].packing_score(probes[p], box, False)
            scores[r, p, frame_index] += s
    frame_index += 1
print("finished calc")

numpy.savetxt(output_filename_core + ".dat", scores.reshape(len(residues), len(probes)*len(traj)), header=header)
