#!/usr/bin/env python
"""
Program to detect lipid chains that have made their way inside a membrane
protein
"""

import sys
import loos
import ConvexHull

if len(sys.argv)==1 or sys.argv[1] == "-h" or sys.argv[1] == "--fullhelp":
    print """
Usage: %s system trajectory protein_selection zmin zmax zblocks target_selection output_directory helix_ranges

    system: system description file, e.g. PDB, PSF, parmtop, etc
    trajectory: trajectory file, e.g. DCD, XTC
    protein_selection: string to select the protein, e.g. 'segid == "PROT" && backbone'
    zmin, zmax: z-range to be considered for lipid penetration events
    zblocks: how many slices along the z-axis the system should be divided into
    target_selection: selection string to pick out the lipids you want to look at
    output_directory: where to write the output files (assumed to already exist)
    helix_ranges: list of colon-separated integers specifying residue ranges for the helices

Example command line:
inside_helices.py rhodopsin.psf sim.dcd 'segid == "RHOD" && backbone' -17 17 10 'segid =~ "[TB]PE" && name =~ "C2\d+"' dha_helices 35:64 71:100 107:139 151:173 200:233 246:277 286:308

What's actually going on:
The purpose of this program is to detect cases where lipid chains or 
headgroups move "inside" a helical protein.  The user specifies what the 
protein is and the residue ranges of the helices, and tells us how to 
identify the lipids (in the given example, we're looking at the C2 lipid
chain from PE lipids).

At each time point, the code breaks the system into slices along z, and
within each slice computes the centroid for each helix; this lets the
program correctly capture helix tilts and kinks, which can alter where the
"surface" of the protein is.  Note: if you set zblocks too big, you won't
have enough atoms for each helix in that block to compute a sane centroid.
It's up to you whether you want to include side chains.

The centroids of the helices are used to compute a convex hull for each
z-slice.  Then, for each lipid chain (or headgroup, etc) we check to
see if there are any atoms that are inside the hull that's in its z range.
If the number of atoms is greater than the threshold (currently hardwired
to 7, TODO: fix this), then that chain is considered to be "inside" the
protein at that time.  

The output of the program is a set of time series: each lipid that is ever
identified as being inside the protein gets its own file (named as
SEGID:resid.dat).  The first column is the frame number, and the second is
1 if the lipid was inside on that frame and 0 otherwise.  

If you just want the distribution of occupancies (e.g. what fraction of
time the lipid was present), you can say something like

grep Occ *.dat | awk '{print $4}' > vals.dat

and then histogram vals.dat.
    

    """ % (sys.argv[0])


system_file = sys.argv[1]
traj_file = sys.argv[2]
protein_string = sys.argv[3]
zmin = float(sys.argv[4])
zmax = float(sys.argv[5])
zblocks = int(sys.argv[6])
target_string = sys.argv[7]
output_directory = sys.argv[8]
helix_ranges = sys.argv[9:]

system = loos.createSystem(system_file)
traj = loos.createTrajectory(traj_file, system)

protein = loos.selectAtoms(system, protein_string)

helices = []
helix_centroids = loos.AtomicGroup()
for h in helix_ranges:
    first, last = h.split(":")
    helix_string ='(resid >= ' + first + ') ' + '&& (resid <= ' + last  + ')'
    helix = loos.selectAtoms(protein, helix_string)
    helices.append(helix)
    helix_centroids.append(loos.Atom())


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
        current_list = []
        for h in range(len(helices)):
            helix = slicers[i](helices[h])
            if len(helix) > 0:
                centroid = helix.centroid()
                helix_centroids[h].coords(centroid)
                current_list.append(helix_centroids[h])

        if len(current_list) < 3:
           print "Warning: only %d helices represented in frame %d and slice %d" % (len(current_list), frame, i)
           hulls.append(None)
        else:
            hull = ConvexHull.ConvexHull(helix_centroids)
            hull.generate_hull()
            hull.generate_vertices()
            hulls.append(hull)

    for chain in chains:
        atoms_inside = 0
        for atom in chain:
            z = atom.coords().z()

            # skip atoms outside the z range
            if not (zmin < z < zmax):
                continue
            
            index = int((z-zmin)/block_size)

            if hulls[index] and hulls[index].is_inside(atom.coords()):
                atoms_inside += 1
                if atoms_inside >= threshold:
                    key = atom.segid() + ":" + str(atom.resid())
                    if not bound_lipids.has_key(key):
                        bound_lipids[key] = []
                    bound_lipids[key].append(frame)
                    break
    frame += 1
    if frame % 20 == 0: print frame

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
