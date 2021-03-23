#!/usr/bin/env python3

import loos
import loos.pyloos
import sys

print("#", " ".join(sys.argv))

system = loos.createSystem(sys.argv[1])
traj = loos.pyloos.Trajectory(sys.argv[2], system)

res1 = loos.selectAtoms(system, sys.argv[3])
res2 = loos.selectAtoms(system, sys.argv[4])

for frame in traj:
    box = system.periodicBox()
    stacking = res1.stacking(res2, box, 5)
    packing = res1.packingScore(res2, box, norm=True)
    print(traj.index(), stacking, packing)
