#!/usr/bin/env python3
"""
Compute packing score between each residue of a protein and a number of
small molecules in solution.

Alan Grossfield, University of Rochester Medical Center
"""

"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2019 Tod Romo, Grossfield Lab
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
import sys
import loos
import loos.pyloos
import numpy
import argparse
from sklearn import decomposition

def fullhelp():
    print("""
    Compute packing score between each residue of a protein and a number of
    small molecules in solution.

    Summary of arguments
    system_file: file describing the full system contents. For this program
                 to work correctly, it must contain connectivity information, so either a PSF or a PDB with CONECT records is best
    protein_selection: selection string for the protein
    probe_selection: selection string for the small molecules (will be
                     split into individual molecules)
    output_core: root name of the output files
    --traj: one or more trajectory files

    Options
    --skip value : number of frames to skip at the front of each trajectory
    --stride value: how to step through the trajectory (default=1)
    --no_h: remove hydrogens from the selection (recommended)
    --pca: perform principal component analysis on the time series

    Packing score is a measure of packing between two groups or moieties
    inspired by the Lennard-Jones function.  It is the sum over all
    inter-moiety pairs of atoms of 1/r^6.  To my knowledge, this term
    was first introduced in Grossfield et al, Proc. Nat. Acad. Sci.
    USA, 2006, 103, 4888-4893.

    Output:
    If output_core is set to "foo", then the following files will be
    produced:

    foo_raw_score.dat: matrix of packing scores.  Each row is the packing
                       score between 1 specific probe molecule against each
                       protein residue from 1 frame.
    foo_ave.dat: column of per-residue summed packing scores

    If --pca is specified, principal component analysis is performed on the
    matrix from foo_raw_score, producing 3 additional files:

      foo_var.dat: the fraction of the signal for each eigenmode
                   (essentially the eigenvalues)
      foo_var_cum.dat: the cumulative variance, so you can assess how many
                       eigenmodes you should analyze)
      foo_comp.dat: each column is an eigenvector,
                    representing an intensity for each residue.  This is
                    useful for identifying residues that simultaneously
                    interact with an individual probe molecule, and thus are
                    candidate binding sites.

    Example command
    packing_score_per_res.py system.psf 'segid == "PROT"' 'resname == "LIG"' score --no_h --pca --traj traj_1.dcd traj_2.dcd

    -- the protein is segid "PROT", while the ligand has residue name "LIG"
    -- analyzes 2 trajectories (traj_1.dcd, traj_2.dcd)
    -- automatically removes hydrogens from the selection
    -- performs pca
    -- output files are score_raw_score.dat, score_ave.dat, etc

    """)

class FullHelp(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        kwargs['nargs']=0
        super(FullHelp, self).__init__(option_strings, dest, **kwargs)
    def __call__(self, parser, namespace, values, option_string = None):
        fullhelp()
        parser.print_help()
        setattr(namespace, self.dest, True)
        parser.exit()



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
parser.add_argument('--fullhelp', help="Print detailed description of all options", action=FullHelp)

args = parser.parse_args()


system = loos.createSystem(args.system_file)

protein = loos.selectAtoms(system, args.protein_selection)
probe = loos.selectAtoms(system, args.probe_selection)
if args.no_h:
    protein = loos.selectAtoms(protein, '!hydrogen')
    probe = loos.selectAtoms(probe, '!hydrogen')

residues = protein.splitByResidue()
probes = probe.splitByMolecule()

frame_index = 0
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
scores = numpy.zeros([len(residues), len(probes), len(vtraj)], float)

frame_index = 0
for frame in vtraj:
    box = frame.periodicBox()
    for r in range(len(residues)):
        for p in range(len(probes)):
            s = residues[r].packingScore(probes[p], box, False)
            scores[r, p, frame_index] += s
    frame_index += 1

print("finished calculating packing scores")

scores = scores.reshape(len(residues), len(probes)*len(vtraj)).transpose()
numpy.savetxt(args.output_core + "_raw_score.dat",
              scores,
              header=header,
              fmt='%.6e')

ave = numpy.add.reduce(scores, axis=0)
ave /= len(vtraj)
resids = numpy.arange(1, len(residues)+1)

numpy.savetxt(args.output_core + "_ave.dat",
              numpy.column_stack((resids, ave)),
              fmt='%.6e',
              header="Residue\tAverage signal")

if args.pca:
    pca = decomposition.PCA()
    pca.fit(scores)
    numpy.savetxt(args.output_core + "_var.dat",
                  numpy.column_stack((resids, pca.explained_variance_ratio_)),
                  fmt='%.6e',
                  header="Mode\tFraction variance")
    numpy.savetxt(args.output_core + "_var_cum.dat",
                  numpy.column_stack((resids,
                                      numpy.add.accumulate(
                                          pca.explained_variance_ratio_))),
                  fmt='%.6e',
                  header="Mode\tCum fraction variance")
    numpy.savetxt(args.output_core + "_comp.dat",
                  numpy.column_stack((resids,
                                      numpy.transpose(pca.components_))),
                  fmt='%.6e',
                  header="Residue\tMode1\tMode2\t...")
