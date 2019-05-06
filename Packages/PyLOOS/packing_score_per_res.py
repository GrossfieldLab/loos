#!/usr/bin/env python3

import sys
import loos
import loos.pyloos
import numpy


# TODO: this should probably be a method on AtomicGroup for performance
def packing_score(residue, probe, norm=False):
    box = residue.periodicBox()
    score = 0.0
    for a1 in residue:
        for a2 in probe:
            dist2 = a2.coords().distance2(a1.coords(), box)
            score += dist2 * dist2 * dist2

    if norm:
        score /= residue.size() * probe.size()

    return score


# Begin the actual program
header = "#" + " ".join(sys.argv)

system_file = sys.argv[1]
traj_file = sys.argv[2]
protein_selection = sys.argv[3]
probe_selection = sys.argv[4]
output_filename = sys.argv[5]


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
    for r in range(len(residues)):
        for p in range(len(probes)):
            s = packing_score(residues[r], probes[p])
            scores[r, p, frame_index] += s
    frame_index += 1

scores.reshape(len(residues), len(probes)*len(traj))

scores.savetxt(output_filename, header=header)
