#!/usr/bin/env python3

import loos
import sys

system_file = sys.argv[1]
input_traj_file = sys.argv[2]
output_traj_file = sys.argv[3]

system = loos.createSystem(system_file)
input_traj = loos.createTrajectory(input_traj_file, system)

output_traj = loos.DCDWriter(output_traj_file)

while (input_traj.readFrame()):
    input_traj.updateGroupCoords(system)
    output_traj.writeFrame(system)