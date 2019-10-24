#!/usr/bin/env python3
"""
   cylindrical-density.py : compute the density of atoms ("targets") around
                            a centering selection, output a 2D histogram in
                            r,z
   The expected use-case for this is looking at lipid packing around a membrane
   protein, when for whatever region you want to average out the "polar"
   degrees of freedom.

   For example, a command line could look like:

   cylindrical-density.py sim.psf sim.dcd 'segid == "PROT"' 'segid =~ "PE" && \!hydrogen' -25 25 50 10 30 20

   This would read the system info from sim.psf, and use the trajectory sim.dcd.
   The system would be translated such that "PROT" is at the origin, and then
   the density of all heavy atoms with segment names containing "PE" would be
   computed.  The z-range of the histogram would go from -25:25 with 50 bins,
   and r-range would be 10:30 with 20 bins.

   I highly suggest including the "&& \!hydrogen" part of the target selection,
   since that will make the program run significantly faster without
   substantially changing the information content (the slash in front of the "!"
   may or may not be necessary, depending on which shell you use).

   Alan Grossfield
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


import sys
import loos
import loos.pyloos
import numpy
import math

header = " ".join(sys.argv)
print("# ", header)

system_file = sys.argv[1]
traj_file = sys.argv[2]
centering_selection_string = sys.argv[3]
target_selection_string = sys.argv[4]
zmin = float(sys.argv[5])
zmax = float(sys.argv[6])
znum_bins = int(sys.argv[7])
rmin = float(sys.argv[8])
rmax = float(sys.argv[9])
rnum_bins = int(sys.argv[10])

system = loos.createSystem(system_file)
traj = loos.pyloos.Trajectory(traj_file, system)

centering = loos.selectAtoms(system, centering_selection_string)
target = loos.selectAtoms(system, target_selection_string)


zbin_width = (zmax - zmin) / znum_bins

rbin_width = (rmax - rmin) / rnum_bins
rmin2 = rmin*rmin
rmax2 = rmax*rmax

hist = numpy.zeros([rnum_bins, znum_bins])

for frame in traj:

    centroid = centering.centroid()
    target.translate(-centroid)
    target.reimageByAtom()

    for atom in target:
        x = atom.coords().x()
        y = atom.coords().y()
        z = atom.coords().z()

        r2 = x*x + y*y

        if (zmin < z < zmax) and (rmin2 < r2 < rmax2):
            r = math.sqrt(r2)

            rbin = int((r - rmin) / rbin_width)
            zbin = int((z - zmin) / zbin_width)

            if (0 <= rbin < rnum_bins) and (0 <= zbin < znum_bins):
                hist[rbin, zbin] += 1.0

hist /= len(traj)

for i in range(rnum_bins):
    rinner = rmin + i*rbin_width
    router = rinner + rbin_width
    rval = rmin + (i+0.5)*rbin_width

    norm = math.pi * (router*router - rinner*rinner)
    for j in range(znum_bins):
        zval = zmin + (j+0.5)*zbin_width
        print(rval, zval, hist[i, j] / norm)
    print()
