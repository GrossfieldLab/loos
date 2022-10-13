#!/usr/bin/env python3
"""
Compute the voronoi cross-sectional area for something (e.g. a protein) through the membrane.

Alan Grossfield, University of Rochester Medical Center
Copyright 2014
"""

import sys
import loos
import loos.pyloos
import numpy
from loos.Voronoi import *

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--fullhelp":
        print("""
Usage:
area_profile.py system trajectory skip stride zmin zmax num_slices padding all-selection-string target-selection-string

system: system file (e.g. pdb, psf, parmtop, gro)
trajectory: trajectory file (periodic boundary information required)
skip, stride: how to move through the trajectory, skipping the first "skip" values and stepping by
        "stride"
zmin, zmax: range of z-values to consider for area calculation
num_slices: number of slices to break the vertical span into
padding: floating point number specifying how many extra layers of atoms are generated
         (see below)
all-selection-string: the set of atoms used to compute the voronoi decomposition
target-selection-string: set of atoms whose areas we report

Notes
    1) all selections are forced to be subsets of the initial selection.  This is
       necessary for the mapping of areas to work correctly.
    2) this program assumes that the system has already been centered such that
       the z location of the membrane isn't drifting (z-slices are absolute, not
       relative to the membrane center) and such that the periodic box
       is centered at x=y=0.

Padding:
    This program uses scipy Voronoi implementation, which is itself a wrapper around
    Qhull.  Qhull doesn't know anything about periodic boundary conditions as far as I
    can tell, so we fake it by generating extra atoms so provide an environment around
    the real ones; otherwise, you get insane areas for the atoms at the "edge" of the
    box.  You need to pick this value such that there are ~2-3 layers of atoms around
    the real ones; in my experience, a value of 15 ang works well if you're selecting
    all atoms or all heavy atoms.  If your selection is more sparse (e.g. just phosphate
    atoms), you will need to make it significantly larger, perhaps 30 ang).  The best
    test is to run with different pad values on a short trajectory and find the value
    where the answer doesn't change.  Since a larger pad value increases the number
    of atoms used in the Voronoi decomposition, it does get more expensive, but it's
    never too bad, and there's always a faster way to get a wrong answer.



Example selection choice:
    '\!hydrogen' 'segid == "RHOD"'
        use all heavy atoms for the voronoi decomposition, then pick out the rhodopsin
        molecule to get its area
        (as a rule, you should get virtually identical answers with all atoms and all
        heavy atoms, but the latter will be dramatically faster)


        """)
        sys.exit(0)
    elif len(sys.argv) != 11 or sys.argv[1] == "-h" or sys.argv[1] == "--h":
        print("Usage: area_profile.py system trajectory skip stride zmin zmax num_slices padding all-selection-string target-selection-string")
        sys.exit(0)

    system_filename = sys.argv[1]
    traj_filename = sys.argv[2]
    skip = int(sys.argv[3])
    stride = int(sys.argv[4])
    zmin = float(sys.argv[5])
    zmax = float(sys.argv[6])
    num_slices = int(sys.argv[7])
    padding = float(sys.argv[8])

    all_selection = sys.argv[9]
    target_selection = sys.argv[10]

    print("# ", " ".join(sys.argv))

    system = loos.createSystem(system_filename)
    pytraj = loos.pyloos.Trajectory(traj_filename, system, skip=skip, stride=stride)


    all_atoms = loos.selectAtoms(system, all_selection)
    target_atoms = loos.selectAtoms(all_atoms, target_selection)

    slicers = []
    slice_width = (zmax - zmin) / num_slices
    for i in range(num_slices):
        low = zmin + i*slice_width
        high = zmin + (i+1)*slice_width
        slicers.append(ZSliceSelector(low, high))

    areas = numpy.zeros([num_slices])
    areas2 = numpy.zeros([num_slices])

    for snap in pytraj:
        system.reimageByAtom()

        for i in range(num_slices):
            all_slice_atoms = slicers[i](all_atoms)
            target_slice_atoms = slicers[i](target_atoms)

            # run voronoi
            v = VoronoiWrapper(all_slice_atoms, padding)
            v.generate_voronoi()

            sr = SuperRegion()
            sr.buildFromAtoms(target_slice_atoms, v)
            a = sr.area()
            areas[i] += a
            areas2[i] += a*a

    # normalize
    areas /= len(pytraj)
    areas2 /= len(pytraj)
    areas2 = numpy.sqrt(areas2 - areas*areas)
    print("# ZSlice\tArea\tDev")
    for i in range(num_slices):
        z = zmin + (i+0.5)*slice_width
        print(z, areas[i], areas2[i])

if __name__ == "__main__":
    main()
    