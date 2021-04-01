#!/usr/bin/env python3
"""
Compute X-ray scattering using the spherical scattering approximation. Performs
both regular intensity and Kratky.

Structure factors from Szalóki, I. Empirical Equations for Atomic Form Factor
and Incoherent Scattering Functions. X Ray Spectrom 1996, 25 (1), 21–28
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

import argparse

def fullhelp():
    print("""
    This program computes the X-ray scattering from a trajectory using the
    spherical scattering approximation:

    I(q) = (sum over atom pairs) Fi(q) Fj(q) sin qd / qd
    where F(q) are the atomic scaterring profiles and d is the distance
    between the atoms.

    Command line:
        scattering.py system_file selection out_file traj_file
            system_file: the file describing system content (eg PDB, PSF)
            selection: loos selection string out atoms to do the scattering
                       calculation on
            out_file: file containing the output (see below)
            traj_file: trajectory file

    Other options:
        --qmin: minimum q value used. Default=0 The intensity curve is
                normalized by the I at this value
        --qmax: maximum q value. Default=6.0
        --nq: number of q values to use. Default=20
        --deduce_element: If specified, attempt to deduce the element of each
                atom from the atom's name. THIS IS NOT VERY ROBUST, and will
                give incorrect guesses for atoms named "CAL" (calcium) or
                "SOD" (sodium). If possible, avoid using this option in favor
                of specifying the system with a file format that gives us the
                element or the mass (eg PSF). It should get proteins, RNA, lipids,
                and water correct.
        --skip: skip this number of frames from the front of the trajectory,
                default=0
        --stride: step through the trajectory by this number, default=1

    Output: The resulting data file is field-delimited text. The first column
        is the q value. Second is I(q)/I(0), the normalized intensity. Third
        is q*Rg, where Rg is the radius of gyration. Fourth column is
        (q*Rg)^2 I(q)/I(0).  Plotting 1 vs 2 gives a standard intensity
        plot (I recommend log scale), while 3 vs 4 gives a Kratky plot

    Note: although the guts of this code are in C++, it's still doing a double
    loop over all pairs of atoms (including self), times the number of q values,
    for each frame. Think carefully about what you want to include in the
    calculation, since this will probably be pretty slow for an explicit atom
    system.

    Also, this program treats each atom essentially independently, and makes
    no effort to account for solvent structure unless the water is explicitly
    included in the calculation. If you want a program that does that, consider
    using something like CRYSOL.

    For more information on the implementation on the atomic form factors, see
    FormFactor.cpp and FormFactor.hpp in the loos distribution. The actual
    scattering calculation is done in the scattering() method to AtomicGroup,
    implemented in AG_numerical.cpp
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



def deduceAtomicNumber(system):
    """
    Given an AtomicGroup system, guess the atomic number based on the atom's
    name. NOT a particularly robust method, but useable if your system file
    doesn't give you the element.
    Only knows hydrogen, carbon, nitrogen, oxygen, sulfur. This covers
    protein, nuclei acid, and lipid, but will probably fail on salt (calcium and
    sodium will probably clash)
    """

    for atom in system:
        n = atom.name()
        if n.startswith("H"):
            atom.atomic_number(1)
        elif n.startswith("C"):
            atom.atomic_number(6)
        elif n.startswith("N"):
            atom.atomic_number(7)
        elif n.startswith("O"):
            atom.atomic_number(8)
        elif n.startswith("P"):
            atom.atomic_number(15)
        elif n.startswith("S"):
            atom.atomic_number(16)


if __name__ == '__main__':

    import loos
    import loos.pyloos
    import numpy as np
    import sys

    header = " ".join(sys.argv)
    parser = argparse.ArgumentParser(description="Compute X-ray scattering")
    parser.add_argument('system_file', help="File describing system contents")
    parser.add_argument('selection', help='Selection of atoms to use in calc')
    parser.add_argument('out_file', help="File to write scattering data to")
    parser.add_argument('traj_file', help="Trajectory file")
    parser.add_argument('--fullhelp',
                        help="Print detailed description of all options",
                        action=FullHelp)
    parser.add_argument('--qmin', default=0.0, type=float,
                        help="Minimum q value")
    parser.add_argument('--qmax', default=6.0, type=float,
                        help="Maximum q value")
    parser.add_argument('--nq', default=20, type=int,
                        help="Number of q values")
    parser.add_argument('--deduce_element', action="store_true",
                        help="Attempt to deduce elements from name")
    parser.add_argument('--skip',
                        help='Skip frames the start of the combined trajectory',
                        type=int,
                        default=0)
    parser.add_argument('--stride', help='Do every nth frame',
                        type=int,
                        default=1)


    args = parser.parse_args()


    q_min = args.qmin
    q_max = args.qmax
    num_qvals = args.nq

    system = loos.createSystem(args.system_file)
    traj = loos.pyloos.Trajectory(args.traj_file, system,
                                  skip=args.skip, stride=args.stride)

    subset = loos.selectAtoms(system, args.selection)

    if args.deduce_element:
        deduceAtomicNumber(subset)

    formFactors = loos.FormFactorSet()

    total = np.zeros([num_qvals])
    q_vals = np.arange(q_min, q_max, (q_max - q_min)/num_qvals)
    rgyr = 0.0
    for frame in traj:
        total += np.asarray(subset.scattering(q_min, q_max, num_qvals,
                                              formFactors))
        rgyr += subset.radiusOfGyration()

    total /= total[0]  # output I/I(0)
    rgyr /= len(traj)

    qrg = q_vals * rgyr
    kratky = qrg * qrg * total

    total = np.column_stack((q_vals, total, qrg, kratky))

    np.savetxt(args.out_file, total, header=header + "\n" + "Q\tI/I0\tQ*Rg\t(Q*Rg)^2 I/I(0)")
