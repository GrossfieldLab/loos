#!/usr/bin/env python
"""
Program to detect lipid chains that have made their way inside a membrane
protein
"""

import sys
import loos
import ConvexHull

system_file = sys.argv[1]
traj_file = sys.argv[2]
protein_string = sys.argv[3]
zmin = float(sys.argv[4])
zmax = float(sys.argv[5])
zblocks = int(sys.argv[6])
target_string = sys.argv[7]
output_directory = sys.argv[8]

system = loos.createSystem(system_file)
traj = loos.createTrajectory(traj_file, system)

protein = loos.selectAtoms(system, protein_string)

target = loos.selectAtoms(system, target_string)
target = loos.selectAtoms(target, "!hydrogen")

chains = target.splitByResidue()

threshold = 7

if (zmin > zmax): 
    tmp = zmax
    zmax = zmin
    zmin = tmp
    
slicers = []
block_size = (zmax - zmin) / zblocks
for i in range(zblocks):
    z1 = zmin + i*block_size
    z2 = zmin + (i+1)*block_size

    ch = ConvexHull.ZSliceSelector(z1, z2)
    slicers.append(ch)


frame = 0

bound_lipids = {}

while (traj.readFrame()):
    traj.updateGroupCoords(system)

    # set up the convex hulls for this slice
    hulls = []
    for i in range(zblocks):
        h = ConvexHull.ConvexHull(slicers[i](protein))
        h.generate_hull()
        h.generate_vertices()
        hulls.append(h)

    for chain in chains:
        atoms_inside = 0
        for atom in chain:
            z = atom.coords().z()

            # skip atoms outside the z range
            if not (zmin < z < zmax):
                continue
            
            index = int((z-zmin)/block_size)

            if hulls[i].is_inside(atom.coords()):
                atoms_inside += 1
                #print frame, atom.resid(), atom.id(), atoms_inside
                if atoms_inside >= threshold:
                    #print frame, atom.resid()
                    key = atom.segid() + ":" + str(atom.resid())
                    if not bound_lipids.has_key(key):
                        bound_lipids[key] = []
                    bound_lipids[key].append(frame)
                    break
    frame += 1

for lipid in bound_lipids.keys():
    frames = bound_lipids[lipid]
    occ = float(len(frames)) / frame

    file = open(output_directory + "/" + lipid + ".dat", "w")
    file.write("# Occupancy = %f\n" % (occ))
    file.write("#Frame\tBound\n")
    for i in range(frame):
        if i in frames:
            file.write("%d\t%d\n" %(i, 1))
        else:
            file.write("%d\t%d\n" %(i, 0))
    file.close()
