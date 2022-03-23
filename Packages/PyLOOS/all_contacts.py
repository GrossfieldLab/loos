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
import sys
import numpy as np
from os.path import basename, splitext
from sklearn import decomposition


fullhelp = """
  all_contacts.py: compute the probability of residue-residue contact
  over the course of a trajectory

  Mandatory arguments:
  model: file describing system contents, e.g. a psf or pdb
  selection: selection string for which residues to look at
  out_file: name for the average contact map, written in matlab format
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

  This program does not explicitly handle periodicity; it assumes you've
  already fixed any periodicity issues before you ran it.

  """


def make_index(i, j, num_res):
    """ Utility function to flatten symmetric matrix
    Assumes the matrix is 0-based
    """
    # canonicalize order
    if i > j:
        i, j = j, i

    index = 0
    for k in range(1, i+1):
        index += num_res - k
    index += j - i - 1
    return index


def get_residues(index, num_res):
    """ Utility function to get symmetric indices from flattened index
    Assumes the matrix is 0-based
    """
    i = 0
    j = 1
    val = 0
    next = 0
    while (val + next <= index):
        val += next
        i += 1
        next = num_res - i
    j = index - val + i
    return i-1, j

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
    # now remove the backbone -- doing before the split loses the glycines
    if args.no_backbone:
        residues = list([loos.selectAtoms(r, "!backbone") for r in residues])

    if args.pca:
        total_frames = 0
        for traj in all_trajs:
            total_frames += len(traj)
        num_pairs = int((len(residues) * (len(residues)-1))/2)
        frac_contacts_frame = np.zeros([num_pairs, total_frames],
                                       np.float64)
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
                            index = make_index(i, j, len(residues))
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
        pca = decomposition.PCA()
        pca.fit()

        # hardwired file names for now
        pairs = np.arange(num_pairs)
        np.savetxt('pca' + "_var.dat",
                   np.column_stack((pairs, pca.explained_variance_ratio_)),
                   fmt='%.6e',
                   header="Mode\tFraction variance")

        np.savetxt('pca' + "_comp.dat",
                   np.column_stack((pairs,
                                   np.transpose(pca.components_))),
                   fmt='%.6e',
                   header="Pair\tMode1\tMode2\t...")

        # write out a mapping of indices in the pca to residue pairs
        with open("index_file", "w") as index_file:
            index_file.write("Index\tRes1\tRes2")
            for index in pairs:
                i, j = get_residues(index, num_pairs)
                index_file.write(index, i, j)
