#!/usr/bin/env python3
"""
Compute areas for different sets of atoms within a particular slice along the membrane
normal.

Alan Grossfield, University of Rochester Medical Center
Copyright 2014
"""

import sys
import loos
from loos.Voronoi import *

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--fullhelp":
        print("""

Usage:
program system trajectory skip zmin zmax padding selection-string1 [selection-string2 ...]

system: system file (e.g. pdb, psf, parmtop, gro)
trajectory: trajectory file (Periodic boundary information required)
zmin, zmax: floating point numbers used to select a particular slice of the system
padding: floating point number specifying how many extra layers of atoms are generated
         (see below)
selection-string1: the set of atoms used to compute the voronoi decomposition
selection-string2, etc: sets of atoms for which areas are reported

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



Example selection choices:
    To map the areas for different membrane constituents, one could do something like
    '\!hydrogen' 'resname == "PCGL"' 'resname == "PALM"' 'resname == "OLEO"' 'resname == "TIP3"'
    The first selection is all heavy atoms, which we then break into contributions from
    the headgroup, palmitate chain, oleate chain, and water (this is the old charmm27
    format).  Note that skipping hydrogens is a good idea -- it won't change the answer
    noticeably, but it greatly reduces the number of atoms, which will make things
    run much faster.

    To track "lipid areas" in a bilayer with POPE and POPG lipids, one could do
    'name == "P"' 'resname == "POPE"' 'resname == "POPG"'



        """)
        sys.exit(0)
    elif len(sys.argv) < 8 or sys.argv[1] == "-h" or sys.argv[1] == "--h":
        print(sys.argv[0], " system trajectory skip zmin zmax padding selection-string1 [selection-string2 ...]")
        sys.exit(0)


    system_filename = sys.argv[1]
    traj_filename = sys.argv[2]
    skip = int(sys.argv[3])
    zmin = float(sys.argv[4])
    zmax = float(sys.argv[5])
    padding = float(sys.argv[6])
    selection_strings = sys.argv[7:] # first selection gives the all atoms to use
                                     # in area calculations, all others tell you
                                     # how to group the areas

    print("# ", " ".join(sys.argv))

    system = loos.createSystem(system_filename)
    traj = loos.createTrajectory(traj_filename, system)

    traj.readFrame(skip)
    traj.updateGroupCoords(system)

    slicer = ZSliceSelector(zmin, zmax)

    selections = []
    selections.append(loos.selectAtoms(system, selection_strings[0]))
    for s in selection_strings[1:]:
        selections.append(loos.selectAtoms(selections[0], s))

    frame = skip
    string = ""
    for i in range(len(selections)):
        string += "\tArea" + str(i)

    print("# Frame", string)
    while (traj.readFrame()):
        traj.updateGroupCoords(system)
        system.reimageByAtom()

        slice_atoms = slicer(selections[0])

        # run voronoi
        v = VoronoiWrapper(slice_atoms, padding)
        v.generate_voronoi()

        areas = []
        # generate the areas for the selections
        for s in selections:  # really should be selections[1:]
            sr = SuperRegion()
            sr.buildFromAtoms(slicer(s), v)
            areas.append(sr.area())
        print(frame, "\t".join(map(str, areas)))
        frame += 1

if __name__ == '__main__':
    main()