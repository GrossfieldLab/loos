#!/usr/bin/env python3
#
#  domain.py: track the motions of 2 portions of a protein, printing out the
#             distance between the centroids, as well as the angle and torsion
#             between the two domains first principal axis
#             NOTE: This was written assuming the two chunks are part of the
#             same molecule, and so doesn't respect periodicity.
#  Alan Grossfield

import sys
import loos
import loos.pyloos
import math

header = " ".join(sys.argv)
print("# ", header)

system_file = sys.argv[1]
traj_file = sys.argv[2]
sel_string1 = sys.argv[3]
sel_string2 = sys.argv[4]

system = loos.createSystem(system_file)
traj = loos.pyloos.Trajectory(traj_file, system)

sel1 = loos.selectAtoms(system, sel_string1)
sel2 = loos.selectAtoms(system, sel_string2)

while traj.nextFrame():

    # compute distance
    centroid1 = sel1.centroid()
    centroid2 = sel2.centroid()

    diff = centroid2 - centroid1
    distance = diff.length()

    # Compute angle between principal axes
    vectors1 = sel1.principalAxes()
    axis1 = vectors1[0]

    vectors2 = sel2.principalAxes()
    axis2 = vectors2[0]
    angle = math.acos(axis1 * axis2) * 180/math.pi

    # Compute torsion between principal axes
    p1 = centroid1 + axis1
    p2 = centroid2 + axis2

    tors = loos.torsion(p1, centroid1, centroid2, p2)

    # write output
    print(traj.index(), distance, angle, tors)
