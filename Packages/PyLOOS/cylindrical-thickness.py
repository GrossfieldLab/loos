#!/usr/bin/env python
"""
   cylindrical-thickness.py : compute the thickness of the membrane around
                            a centering selection, output the average as 
                            a function of lateral distance from origin
   The expected use-case for this is looking at lipid packing around a membrane
   protein, when for whatever region you want to average out the "polar" 
   degrees of freedom.  

   For example, a command line could look like:

   cylindrical-thickness.py sim.psf sim.dcd 'segid == "PROT"' 'segid =~ "PE" && name == "P"' 10 30 20 

   This would read the system info from sim.psf, and use the trajectory sim.dcd.
   The system would be translated such that "PROT" is at the origin, and then
   the heights of phosphate with segment names containing "PE" would be 
   computed. 

   NOTE: This code assumes the membrane is centered at z = 0. Atoms with 
   z>0 are assigned to the upper leaflet, z<0 the lower leaflet, and the 
   difference between the averages is computed as a function of r.

   
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
print "# ", header

system_file = sys.argv[1]
traj_file = sys.argv[2]
centering_selection_string = sys.argv[3]
target_selection_string = sys.argv[4]
rmin = float(sys.argv[5])
rmax = float(sys.argv[6])
rnum_bins = int(sys.argv[7])

system = loos.createSystem(system_file)
traj = loos.pyloos.Trajectory(traj_file, system)

centering = loos.selectAtoms(system, centering_selection_string)
target = loos.selectAtoms(system, target_selection_string)


rbin_width = (rmax - rmin) / rnum_bins
rmin2 = rmin*rmin
rmax2 = rmax*rmax

upper_sum = numpy.zeros(rnum_bins)
upper_count = numpy.zeros(rnum_bins)
lower_sum = numpy.zeros(rnum_bins)
lower_count = numpy.zeros(rnum_bins)

for frame in traj:

    centroid = centering.centroid()
    target.translate(-centroid)
    target.reimageByAtom()

    for atom in target:
        x = atom.coords().x()
        y = atom.coords().y()
        z = atom.coords().z()

        r2 = x*x + y*y

        if (rmin2 < r2 < rmax2):
            r = math.sqrt(r2)
    
            rbin = int((r - rmin)/ rbin_width)
            if z > 0:
                upper_sum[rbin] += z
                upper_count[rbin] += 1
            else:
                lower_sum[rbin] += z
                lower_count[rbin] += 1


print "#r   Thick   UpperHeight LowerHeight"
for i in range(rnum_bins):
    r = rmin + (i+0.5)*rbin_width
    upper = lower = 0.0
    if upper_count[rbin] > 0:
        upper = upper_sum[i] / upper_count[i]
    if lower_count[rbin] > 0:
        lower = lower_sum[i] / lower_count[i]
    print r, upper-lower, upper, lower
