#!/usr/bin/env python3

import loos
import loos.pyloos

import numpy as np
import sys



system_file = sys.argv[1]
selection_string = sys.argv[2]
traj_file = sys.argv[3]

q_min = 0
q_max = 6
num_qvals = 60

system = loos.createSystem(system_file)
traj = loos.pyloos.Trajectory(traj_file, system)

subset = loos.selectAtoms(system, selection_string)

formFactors = loos.FormFactorSet()

total = np.zeros([num_qvals])
q_vals = np.arange(q_min, q_max, (q_max - q_min)/num_qvals)
for frame in traj:
    scattering = subset.scattering(q_min, q_max, num_qvals, formFactors)
