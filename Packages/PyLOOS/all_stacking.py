#!/usr/bin/env python3

import sys
import loos
import loos.pyloos
import numpy as np

system_file = sys.argv[1]
traj_file = sys.argv[2]
selection_string = sys.argv[3]
outfile_name = sys.argv[4]

header = "#" + " ".join(sys.argv)

box = loos.GCoord(1000., 1000., 1000.)

system = loos.createSystem(system_file)
traj = loos.pyloos.Trajectory(traj_file, system)

selection = loos.selectAtoms(system, selection_string)
# TODO: make doing this optional
selection = loos.selectAtoms(selection, '!backbone && name =~"^C"')

residues = selection.splitByResidue()

scores = np.zeros((len(residues), len(residues)))

for frame in traj:
    for i in range(len(residues)-1):
        for j in range(i+1, len(residues)):
            p = residues[i].stacking(residues[j], box, 5.0)
            scores[i, j] += p
            scores[j, i] += p
scores /= len(traj)

np.savetxt(outfile_name, scores, header=header)
