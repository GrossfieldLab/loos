#!/usr/bin/env python3
"""
Given a system file, write a new PDB file with the segid field written in
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
import argparse


def fullhelp():
    print("""
    The purpose of this program is, given a system file that does not use
    segment IDs, to create a PDB file that does have them.

    Example: set_segids.py old.pdb new.pdb --selection FIR 'resid < 100' --selection SEC 'resid >= 100 && resid < 200' --selection THI 'resid >= 200'

    This would read old.pdb, create 3 selections, and set their segids to FIR,
    SEC, and THI, respectively. It would then write the new file out as
    new.pdb.

    Atom order and coordinates remain unchanged. The only purpose of this
    program is put add segids, which can be useful for selection purposes. You
    can specify as many selections as you want.

    Having a mix where some atoms have a segid and others don't seems like a
    bad idea, so atoms not matched by any of the selections will be assigned a
    default segid, "OTH". You can override that choice with the --default
    option.

    The segid field is supposed to be 4 characters in length, so any segids
    specified will be truncated at 4 characters. The code does not check
    whether truncation will cause 2 segids to be identical.

    Note: segids are generally expected to be contiguous blocks of atoms, but
    this program doesn't make any attempt to enforce that. If you want to break
    that assumption, you can, but I don't know if that would break other things
    downstream.
        """)


class FullHelp(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        kwargs['nargs'] = 0
        super(FullHelp, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        fullhelp()
        parser.print_help()
        setattr(namespace, self.dest, True)
        parser.exit()


# Begin main program #########
parser = argparse.ArgumentParser(description="Add segids to system file")
parser.add_argument('system_file', help="File describing the system")
parser.add_argument('output_file', help="Name of PDB file to be created")
parser.add_argument('--selection',
                    nargs=2,
                    action='append',
                    help="segid selection-string"
                    )
parser.add_argument('--default',
                    nargs=1,
                    default="OTH",
                    help="default segid to apply to unmatched atoms"
                    )
parser.add_argument('--fullhelp',
                    help="Print detailed description of all options",
                    action=FullHelp)
args = parser.parse_args()

system = loos.createSystem(args.system_file)

all_selections = loos.AtomicGroup()
for seg, sel in args.selection:
    selection = loos.selectAtoms(system, sel)
    # SEGIDs are supposed to be 4 characters
    if len(seg) > 4:
        seg = seg[:4]
    for a in selection:
        a.segid(seg)

    all_selections.append(selection)

remaining_atoms = system.clone()
remaining_atoms.remove(all_selections)
print(len(remaining_atoms), " atoms unlabeled, applying default segid ",
      args.default)
for a in remaining_atoms:
    a.segid(args.default)

pdb = loos.PDB.fromAtomicGroup(system)
with open(args.output_file, "w") as pdb_file:
    pdb_file.write(str(pdb))
