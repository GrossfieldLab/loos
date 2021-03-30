#!/usr/bin/env python3

import loos
import loos.pyloos

import numpy as np
import sys

system_file = sys.argv[1]
selection_string = sys.argv[2]
outfile_name = sys.argv[3]
traj_file = sys.argv[4]

q_min = 0.5
q_max = 6
num_qvals = 12

system = loos.createSystem(system_file)
traj = loos.pyloos.Trajectory(traj_file, system)

subset = loos.selectAtoms(system, selection_string)

formFactors = loos.FormFactorSet()

total = np.zeros([num_qvals])
q_vals = np.arange(q_min, q_max, (q_max - q_min)/num_qvals)
for frame in traj:
    total += np.asarray(subset.scattering(q_min, q_max, num_qvals, formFactors))

total /= len(traj)

total = np.column_stack((q_vals, total))

np.savetxt(outfile_name, total)
