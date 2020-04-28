#!/usr/bin/env python3
"""
Program to detect lipid chains that have made their way inside a membrane
protein
"""

import sys
import loos
from loos.pyloos import ConvexHull
import argparse
import os

fullhelp= """

Example command line:
inside_helices.py rhodopsin.psf sim.dcd 'segid == "RHOD" && backbone' -17 17 10 'segid =~ "[TB]PE" && name =~ "C2\d+"' 35:64 71:100 107:139 151:173 200:233 246:277 286:308 -d dha_helices

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
z-slice.  Then, for each lipid chain (or headgroup, etc) we check to see if
there are any atoms that are inside the hull that's in its z range.  If the
number of atoms is greater than the threshold (default value is 7), then that
chain is considered to be "inside" the protein at that time.   Note: you must
specify at least 3 helices for this program to make any sense, since you can't
have a convex hull with 2 points.

The output of the program is a set of time series: each lipid that is ever
identified as being inside the protein gets its own file (named as
SEGID:resid.dat).  The first column is the frame number, and the second is
1 if the lipid was inside on that frame and 0 otherwise.

If you just want the distribution of occupancies (e.g. what fraction of
time the lipid was present), you can say something like

grep Occ *.dat | awk '{print $4}' > vals.dat

and then histogram vals.dat.
    """


parser = argparse.ArgumentParser(prog=sys.argv[0], epilog=fullhelp, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("system_file",
                    help="File describing system contents, e.g. PDB or PSF")
parser.add_argument("traj_file", help="Trajectory file")
parser.add_argument("protein_string", help="Selection for the protein")
parser.add_argument("zmin",
                    type=float,
                    help="Minimum z value to consider for target atoms")
parser.add_argument("zmax",
                    type=float,
                    help="Maxmimum z value to consider for target atoms")
parser.add_argument("zblocks",
                    type=int,
                    help="Number of blocks to slice the system into along the z-axis")
parser.add_argument("target_string",
                    help="Selection string for the stuff that may penetrate")
parser.add_argument("helix_ranges",
                    nargs="+",
                    help="list of colon-delimited residue number ranges specifying where the helices are")
parser.add_argument("-d", "--directory",
                    help="Directory to create output files",
                    default=".")
parser.add_argument("-t", "--threshold",
                    default=7,
                    type=int,
                    help="Number of atoms inside for the chain to be considered inside")

args = parser.parse_args()

system = loos.createSystem(args.system_file)
traj = loos.createTrajectory(args.traj_file, system)

protein = loos.selectAtoms(system, args.protein_string)

output_directory = args.directory
if not os.path.exists(output_directory):
    try:
        os.mkdir(output_directory)
    except OSError as inst:
        print('Error creating output directory %s : ' % output_directory)
        print(inst)
        sys.exit(1)
if not os.access(output_directory, os.W_OK):
    print("Error: no permission to write to output directory ", output_directory)
    sys.exit(1)


helices = []
helix_centroids = loos.AtomicGroup()
for h in args.helix_ranges:
    first, last = h.split(":")
    helix_string ='(resid >= ' + first + ') ' + '&& (resid <= ' + last + ')'
    helix = loos.selectAtoms(protein, helix_string)
    helices.append(helix)
    helix_centroids.append(loos.Atom())


target = loos.selectAtoms(system, args.target_string)
target = loos.selectAtoms(target, "!hydrogen")

chains = target.splitByResidue()

if (args.zmin > args.zmax):
    tmp = args.zmax
    args.zmax = args.zmin
    args.zmin = tmp

slicers = []
block_size = (args.zmax - args.zmin) / args.zblocks
for i in range(args.zblocks):
    z1 = args.zmin + i*block_size
    z2 = args.zmin + (i+1)*block_size

    ch = ConvexHull.ZSliceSelector(z1, z2)
    slicers.append(ch)


frame = 0

bound_lipids = {}

while (traj.readFrame()):
    traj.updateGroupCoords(system)

    # set up the convex hulls for this slice
    hulls = []
    for i in range(args.zblocks):
        current_list = []
        for h in range(len(helices)):
            helix = slicers[i](helices[h])
            if len(helix) > 0:
                centroid = helix.centroid()
                helix_centroids[h].coords(centroid)
                current_list.append(helix_centroids[h])

        if len(current_list) < 3:
            print("Warning: only %d helices represented in frame %d and slice %d" % (len(current_list), frame, i))
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
            if not (args.zmin < z < args.zmax):
                continue

            index = int((z-args.zmin)/block_size)

            if hulls[index] and hulls[index].is_inside(atom.coords()):
                atoms_inside += 1
                if atoms_inside >= args.threshold:
                    key = atom.segid() + ":" + str(atom.resid())
                    if key not in bound_lipids:
                        bound_lipids[key] = []
                    bound_lipids[key].append(frame)
                    break
    frame += 1
    if frame % 20 == 0: print(frame)

for lipid in list(bound_lipids.keys()):
    frames = bound_lipids[lipid]
    occ = float(len(frames)) / frame

    file = open(output_directory + "/" + lipid + ".dat", "w")
    file.write("# Occupancy = %f\n" % (occ))
    file.write("#Frame\tBound\n")
    for i in range(frame):
        if i in frames:
            file.write("%d\t%d\n" % (i, 1))
        else:
            file.write("%d\t%d\n" % (i, 0))
    file.close()
