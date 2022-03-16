#!/usr/bin/env python3
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

import sys
import loos
import loos.pyloos
import loos.pyloos.options as options
import numpy
import math


if __name__ == '__main__':

    fullhelp = """
cylindrical-density.py : compute the density of atoms ("targets") around
                        a centering selection, output a 2D histogram in
                        r,z
The expected use-case for this is looking at lipid packing around a membrane
protein, when for whatever reason you want to average out the "polar"
degrees of freedom.

For example, a command line could look like:

cylindrical-density.py --model sim.psf --traj sim.dcd --sel 'segid == "PROT"' \
    --target_selection 'segid =~ "PE" && \!hydrogen' --zmin -25 --zmanx 25  \
    --zbins 50 --rmin 10 --rmax 30 --rebins 20

(the line breaks marked by "\" are for readability purposes only)

This would read the system info from sim.psf, and use the trajectory sim.dcd.
The system would be translated such that "PROT" is at the origin, and then
the density of all heavy atoms with segment names containing "PE" would be
computed.  The z-range of the histogram would go from -25:25 with 50 bins,
and r-range would be 10:30 with 20 bins.

I highly suggest including the "&& \!hydrogen" part of the target selection,
since that will make the program run significantly faster without
substantially changing the information content (the slash in front of the "!"
may or may not be necessary, depending on which shell you use).
"""

    lo = options.LoosOptions("Compute 2D histogram in r,z about a selection",
                             fullhelp)

    lo.modelSelectionOptions()
    lo.trajOptions()
    lo.parser.add_argument('--target_selection',
                           help="Atoms whose density we're tracking")
    lo.parser.add_argument('--zmax', type=float,
                           help="Upper z limit of density")
    lo.parser.add_argument('--zmin', type=float,
                           help="Lower z limit of density")
    lo.parser.add_argument('--rmin', type=float,
                           help="Lower r limit of density")
    lo.parser.add_argument('--rmax', type=float,
                           help="Upper r limit of density")
    lo.parser.add_argument('--zbins', type=int,
                           help="Number of bins in z")
    lo.parser.add_argument('--rbins', type=int,
                           help="Number of bins in r")

    args = lo.parse_args()

    system = loos.createSystem(args.model)
    traj = loos.pyloos.Trajectory(args.traj, system,
                                  skip=args.skip, stride=args.stride)

    centering = loos.selectAtoms(system, args.selection)

    target = loos.selectAtoms(system, args.target_selection)


    zbin_width = (args.zmax - args.zmin) / arg.zbins

    rbin_width = (args.rmax - args.rmin) / args.rbins
    rmin2 = args.rmin * args.rmin
    rmax2 = args.rmax * args.rmax

    hist = numpy.zeros([args.rbins, args.zbins])

    for frame in traj:

        centroid = centering.centroid()
        target.translate(-centroid)
        target.reimageByAtom()

        for atom in target:
            x = atom.coords().x()
            y = atom.coords().y()
            z = atom.coords().z()

            r2 = x*x + y*y

            if (args.zmin < z < args.zmax) and (rmin2 < r2 < rmax2):
                r = math.sqrt(r2)

                rbin = int((r - rmin) / rbin_width)
                zbin = int((z - zmin) / zbin_width)

                if (0 <= rbin < args.rbins) and (0 <= zbin < args.zbins):
                    hist[rbin, zbin] += 1.0

    hist /= len(traj)

    for i in range(args.rbins):
        rinner = args.rmin + i*rbin_width
        router = rinner + rbin_width
        rval = args.rmin + (i+0.5)*rbin_width

        norm = math.pi * (router*router - rinner*rinner)
        for j in range(args.zbins):
            zval = args.zmin + (j+0.5)*zbin_width
            print(rval, zval, hist[i, j] / norm)
        print()
