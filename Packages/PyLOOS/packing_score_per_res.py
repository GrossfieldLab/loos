#!/usr/bin/env python3

import sys
import loos
import loos.pyloos
import numpy
import argparse
from sklearn import decomposition


header = "#" + " ".join(sys.argv)
print(header)

parser = argparse.ArgumentParser(description="Compute per-residue packing score for small molecules")
parser.add_argument('system_file',
                    help="System file")
parser.add_argument('protein_selection',
                    help="Selection string identifying the protein")
parser.add_argument('probe_selection',
                    help="Selection string identifying the small molecules")
parser.add_argument('output_core',
                    help="Core name for the output files")
parser.add_argument('--traj',
                    help="One or more trajectory files",
                    nargs='+'
                    )

parser.add_argument('--skip',
                    help="Number of frames to skip for each traj",
                    type=int,
                    default=0)
parser.add_argument('--stride',
                    help="Step size for striding through trajectories",
                    type=int,
                    default=1)
parser.add_argument('--no_h',
                    help="Exclude hydrogens from calculation",
                    action='store_true')
parser.add_argument('--pca',
                    help="Perform PCA on the per-residue profiles",
                    action='store_true')
args = parser.parse_args()


system = loos.createSystem(args.system_file)

protein = loos.selectAtoms(system, args.protein_selection)
probe = loos.selectAtoms(system, args.probe_selection)
if args.no_h:
    protein = loos.selectAtoms(protein, '!hydrogen')
    probe = loos.selectAtoms(probe, '!hydrogen')

residues = protein.splitByResidue()
probes = probe.splitByMolecule()

traj = loos.pyloos.Trajectory(args.traj[0],
                              system,
                              skip=args.skip,
                              stride=args.stride)
vtraj = loos.pyloos.VirtualTrajectory(traj)
for t in args.traj[1:]:
    traj = loos.pyloos.Trajectory(t,
                                  system,
                                  skip=args.skip,
                                  stride=args.stride)
    vtraj.append(traj)

# pre-allocate storage as a numpy array
scores = numpy.zeros([len(residues), len(probes), len(vtraj)], numpy.float)

frame_index = 0
for frame in vtraj:
    box = frame.periodicBox()
    for r in range(len(residues)):
        for p in range(len(probes)):
            s = residues[r].packing_score(probes[p], box, False)
            scores[r, p, frame_index] += s
    frame_index += 1
print("finished calculating packing scores")

scores = scores.reshape(len(residues), len(probes)*len(vtraj)).transpose()
numpy.savetxt(args.output_core + "_raw_score.dat", scores, header=header)

ave = numpy.add.reduce(scores, axis=0)
numpy.savetxt(args.output_core + "_ave.dat", ave)

if args.pca:
    pca = decomposition.PCA()
    pca.fit(scores)
    numpy.savetxt(args.output_core + "_var.dat",
                  pca.explained_variance_ratio_)
    numpy.savetxt(args.output_core + "_var_cum.dat",
                  numpy.add.accumulate(pca.explained_variance_ratio_))
    numpy.savetxt(args.output_core + "_comp.dat",
                  numpy.transpose(pca.components_))
