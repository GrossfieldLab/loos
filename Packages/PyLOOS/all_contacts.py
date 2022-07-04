#!/usr/bin/env python3
"""
Track a set of contacts through a trajectory.  Intended for use with a protein
or RNA, to track all residue-residue contacts within the trajectory.
"""

"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2013 Tod Romo, Grossfield Lab
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

import loos
import loos.pyloos
import loos.pyloos.options as options
from loos.pyloos.SymmMatrix import SymmMatrix
import sys
import numpy as np
from os.path import basename, splitext
from sklearn import decomposition


fullhelp = """
  all_contacts.py: compute the probability of residue-residue contact
  over the course of a trajectory

  Mandatory arguments:
  -m  model: file describing system contents, e.g. a psf or pdb
  --selection  'selection': selection string for which residues to look at
  traj: 1 or more trajectory files

  Options
  --cutoff: distance for atom-atom contacts, defaults to 4.0 Ang
  --no-hydrogens: ignore hydrogens when looking for contacts.
  --no-backbone: ignore the protein or RNA backbone.  If this flag isn't
                 given, you'll have a contact probability of 1 between
                 consecutive residues in a chain
  --individual: if more than 1 trajectory was given, write probability
                files for each one as well as the total.  The files will
                be written in the working directory, with the file name
                the same as the trajectory file, replacing the extension with
                ".dat", so "foo/bar/baz.dcd" would produce "./baz.dat".
                This flag has no effect if only one trajectory was given
                on the command line.
  --fullhelp: produce this message
  --pca:        perform PCA in contact space
  --ncomp:      specify how many PCA modes are computed (default=10)
  --outfile:    if --pca is not specified, the name of the contact map file.
                if --pca is specified, the prefix for the 3 files written.

  This program does not explicitly handle periodicity; it assumes you've
  already fixed any periodicity issues before you ran it.

  If the --pca flag is not specified, then the average contact map is output.
  The file will contain a symmetric matrix indexed by the residues used in the
  calculation (e.g. the entry [6, 10] will be the fraction of frames in which
  the 7th and 11th residues in the selection are in contact).

  PCA calculation

  The idea behind the --pca flag is to look for sets of pairs of residues that
  co-vary (tend to be formed at the same time). This should reveal residual
  structures in highly variable proteins that might not show up in
  clustering due to issues superposing flexible chains. We implement this
  by taking the contact map from each indivual frame, flattening it to 1
  dimension, and using singular value decomposition to perform PCA. Since this
  can be an expensive operation and contact maps change slowly, we recommend
  using the '-s' option to stride through the trajectory.

  If the --pca flag is specified, the contact map is not written. Instead, 3
  files are produced: (assume you specified "--outfile foo"):
        foo_var.dat: 2 column text file. First column is the mode number,
                     second is the fraction of the variance in the data
                     explained by that mode
        foo_comp.dat: each column of the file is a single eigenmode, each row
                      designates a specific pair of residues. The high values
                      in a column indicate pairs of residues that tend to be
                      in contact at the same time
        foo.index: maps rows of foo_comp.dat to pairs of residues. It is
                   0-based and reflects the residues selected.  For example,
                   if you specified 'resid >=5 && resid <=9' the entry
                   2       0       3
                   would say that the 3rd entry of foo_comp.dat would reflect
                   contact between the 1st and 4th residues specified (in this
                   case, residues 5 and 8)

  The "--ncomp X" says that only the first X modes of the eigendecomposition
  should be computed. The default value is 10 -- computing the full svd would
  be costly in time and memory, and only the first few modes are likely to be
  meaningful. This flag is meaningless if the --pca option isn't also supplied.
  """


if __name__ == '__main__':

    lo = options.LoosOptions("Compute probability of residue-residue contacts",
                             fullhelp)
    lo.modelSelectionOptions()
    lo.trajOptions()

    lo.add_argument('--out_file',
                    default='outfile',
                    help="File with the average contact occupancies")
    lo.add_argument('--cutoff', type=float,
                    help="Cutoff distance for contact", default=4.0)
    # TODO: add a number of contacts option
    lo.add_argument('--no_hydrogens', action='store_true',
                    help="Don't include hydrogens")
    lo.add_argument('--no_backbone', action='store_true',
                    help="Don't include the backbone")
    lo.add_argument('--individual', action='store_true',
                    help="Write contact maps for each trajectory")
    lo.add_argument('--pca', action='store_true',
                    help="Perform PCA on the residue-residue maps")
    lo.add_argument('--ncomp',
                    default=10,
                    type=int,
                    help="Number of PCA components to compute")
    args = lo.parse_args()
    header = lo.header()

    system = loos.createSystem(args.model)
    all_trajs = []
    out_names = []
    num_trajs = len(args.traj)
    for t in args.traj:
        traj = loos.pyloos.Trajectory(t, system)
        all_trajs.append(traj)
        if (num_trajs > 1) and args.individual:
            t_base = basename(t)
            core, ext = splitext(t_base)
            out_names.append(core + ".dat")

    if args.no_hydrogens:
        no_hydrogens = loos.selectAtoms(system, "!hydrogen")
        target = loos.selectAtoms(no_hydrogens, args.selection)
    else:
        target = loos.selectAtoms(system, args.selection)

    residues = target.splitByResidue()
    num_res = len(residues)
    # now remove the backbone -- doing before the split loses the glycines
    if args.no_backbone:
        residues = list([loos.selectAtoms(r, "!backbone") for r in residues])

    if args.pca:
        total_frames = 0
        for traj in all_trajs:
            total_frames += len(traj)
        num_pairs = int((num_res * (num_res-1))/2)
        frac_contacts_frame = np.zeros([num_pairs, total_frames],
                                       np.float64)
        symm_indexer = SymmMatrix(num_res)
    else:
        frac_contacts = np.zeros([len(residues), len(residues), num_trajs],
                                 np.float64)

    for traj_id in range(num_trajs):
        traj = all_trajs[traj_id]
        current_frame = 0
        for frame in traj:
            for i in range(len(residues)):
                for j in range(i+1, len(residues)):
                    if residues[i].contactWith(args.cutoff, residues[j]):
                        if args.pca:
                            index = symm_indexer.toFlat(i, j)
                            frac_contacts_frame[index, current_frame] = 1
                        else:
                            frac_contacts[i, j, traj_id] += 1.0
                            frac_contacts[j, i, traj_id] += 1.0
        if not args.pca:
            frac_contacts[:, :, traj_id] /= len(traj)
            if (num_trajs > 1) and args.individual:
                np.savetxt(out_names[traj_id], frac_contacts[:, :, traj_id],
                           header=header)

    # either output the averages or output the pca but not both
    # TODO: write new code to get the average from the larger data set
    # so we can do both
    if not args.pca:
        average = np.add.reduce(frac_contacts, axis=2)
        average /= len(args.traj)

        np.savetxt(args.out_file, average, header=header)

    # do pca if requested
    else:
        pca = decomposition.PCA(n_components=args.ncomp)

        pca.fit(frac_contacts_frame.transpose())

        labels = np.arange(pca.n_components)
        np.savetxt(args.out_file + "_var.dat",
                   np.column_stack((labels,
                                    pca.explained_variance_ratio_)),
                   fmt='%.6e',
                   header=header + "\n" + "Mode\tFraction variance")

        np.savetxt(args.out_file + "_comp.dat",
                   pca.components_.transpose(),
                   fmt='%.6e',
                   header=header + "\n" + "Mode1\tMode2\t...")

        # write out a mapping of indices in the pca to residue pairs
        with open(args.out_file + ".index", "w") as index_file:

            index_file.write("#Index\tRes1\tRes2\n")
            for index in range(pca.components_.shape[1]):
                i, j = symm_indexer.toSymm(index)
                index_file.write(f"{index}\t{i}\t{j}\n")
