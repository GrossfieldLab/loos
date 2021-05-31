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
import numpy
from os.path import basename, splitext


fullhelp = """
  all_contacts.py: compute the probability of residue-residue contact
  over the course of a trajectory

  Mandatory arguments:
  system_file: file describing system contents, e.g. a psf or pdb
  selection: selection string for which residues to look at
  out_file: name for the average contact map, written in matlab format
  traj_files: 1 or more trajectory files

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

lo = options.LoosOptions(fullhelp)
lo.modelSelectionOptions()
lo.trajOptions()


lo.parser.add_argument('--out_file',
                       required=True,
                       help="File with the average contact occupancies")
lo.parser.add_argument('--cutoff', type=float,
                       help="Cutoff distance for contact", default=4.0)
# TODO: add a number of contacts option
lo.parser.add_argument('--no_hydrogens', action='store_true',
                       help="Don't include hydrogens")
lo.parser.add_argument('--no_backbone', action='store_true',
                       help="Don't include the backbone")
lo.parser.add_argument('--individual', action='store_true',
                       help="Write contact maps for each trajectory")
args = lo.parser.parse_args()


header = lo.header()


system = loos.createSystem(args.system_file)
all_trajs = []
out_names = []
num_trajs = len(args.traj_files)
for t in args.traj_files:
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

frac_contacts = numpy.zeros([len(residues), len(residues), num_trajs],
                            numpy.float)


for traj_id in range(num_trajs):
    traj = all_trajs[traj_id]
    for frame in traj:
        for i in range(len(residues)):
            for j in range(i+1, len(residues)):
                if residues[i].contactWith(args.cutoff, residues[j]):
                    frac_contacts[i, j, traj_id] += 1.0
                    frac_contacts[j, i, traj_id] += 1.0
    frac_contacts[:, :, traj_id] /= len(traj)
    if (num_trajs > 1) and args.individual:
        numpy.savetxt(out_names[traj_id], frac_contacts[:, :, traj_id],
                      header=header)

average = numpy.add.reduce(frac_contacts, axis=2)
average /= len(args.traj_files)

numpy.savetxt(args.out_file, average, header=header)
