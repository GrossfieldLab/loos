#!/usr/bin/env python3
"""
Compute distribution of areas/molecule for a z-slice

Alan Grossfield, University of Rochester Medical Center
Copyright 2014
"""

import sys
import loos
import loos.pyloos
import loos.Voronoi
import numpy
from loos.Voronoi import *

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--fullhelp":
        print("""

The purpose of this program is to calculate histograms of areas for different
components of a membrane.  You might use this if you were looking at a
multicomponent bilayer, and wanted to know how much area is taken up
by PC lipids vs. PE lipids. Despite its name, this program does not split
your selections into individual molecules!

Usage:

program system trajectory skip zmin zmax padding min_area max_area num_area_bins selection-string1 selection-string2 [selection-string3 ...]

system: system file (e.g. psf) -- MUST CONTAIN CONNECTIVITY INFORMATION
trajectory: trajectory file (must have periodic boundary information)
zmin, zmax: floating point numbers used to select a particular slice of the system
padding: floating point number specifying how many extra layers of atoms are generated
         (see below)
min_area, max_area, num_area_bins: specifications for histograms of area
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

    This program uses the scipy Voronoi implementation, which is itself a
    wrapper around Qhull.  Qhull doesn't know anything about periodic boundary
    conditions as far as I can tell, so we fake it by generating extra atoms so
    provide an environment around the real ones; otherwise, you get insane
    areas for the atoms at the "edge" of the box.  You need to pick this value
    such that there are ~2-3 layers of atoms around the real ones; in my
    experience, a value of 15 ang works well if you're selecting all atoms or
    all heavy atoms.  If your selection is more sparse (e.g. just phosphate
    atoms), you will need to make it significantly larger, perhaps 30 ang).
    The best test is to run with different pad values on a short trajectory and
    find the value where the answer doesn't change.  Since a larger pad value
    increases the number of atoms used in the Voronoi decomposition, it does
    get more expensive, but it's never too bad, and there's always a faster way
    to get a wrong answer.



Example selection choices:
    To compute the distribution of areas/lipid molecule, where there are
    2 kinds of lipids, one with segids L1, L2, ... L120, and the other with
    segids L121, L121, ... you could do 3 selections like

    'segid =~ "^L\d+" && \!hydrogen' 'segid -> "^L(\d+)" <= 120' 'segid -> "^L(\d+)" > 120'

Note that for a bilayer you have to think carefully about the z-range you use.
If you use a z-range of [0, 20], you'll get the upper leaflet accurately, but
you'll also pick up some junk at very lower area due to the other leaflet
"leaking" across

If you see lines that look like "#Area outside range" followed by some numbers,
it means there was a molecule that had an area outside the range you set for
the histogram.  Either you need to adjust your histogram bounds, or (if the
area is absurdly large) it could suggest your padding value is too small.
        """)
        sys.exit(0)
    elif len(sys.argv) < 11 or sys.argv[1] == "-h" or sys.argv[1] == "--h":
        print(sys.argv[0], " system trajectory skip zmin zmax padding min_area max_area num_area_bins selection-string1 selection-string2 [selection-string3 ...]")
        sys.exit(0)

    system_filename = sys.argv[1]
    traj_filename = sys.argv[2]
    skip = int(sys.argv[3])
    zmin = float(sys.argv[4])
    zmax = float(sys.argv[5])
    padding = float(sys.argv[6])
    min_area = float(sys.argv[7])
    max_area = float(sys.argv[8])
    num_bins = int(sys.argv[9])
    selection_strings = sys.argv[10:] # first selection gives the all atoms to use
                                      # in area calculations, all others tell you
                                      # how to group the areas
    bin_width = (max_area - min_area) / num_bins

    #  Validate that there are multiple selection strings -- realistically, you
    #  need 2
    if len(selection_strings) < 2:
        print("""
        You need at least 2 selection strings -- the first specifies
        print("the atoms used in the Voronoi decomposition, and the others
        are components whose area gets reported (eg different lipid types)
        """)
        sys.exit(0)

    histograms = numpy.zeros([len(selection_strings[1:]), num_bins], float)

    print("# ", " ".join(sys.argv))

    system = loos.createSystem(system_filename)
    pytraj = loos.pyloos.Trajectory(traj_filename, system, skip=skip)

    slicer = ZSliceSelector(zmin, zmax)

    selections = []
    selections.append(loos.selectAtoms(system, selection_strings[0]))
    for s in selection_strings[1:]:
        selections.append(loos.selectAtoms(selections[0], s))
    for i in range(1, len(selections)):
        selections[i] = selections[i].splitByMolecule()

    string = ""
    for i in range(len(selections)):
        string += "\tArea" + str(i)

    for snap in pytraj:
        system.reimageByAtom()

        slice_atoms = slicer(selections[0])

        # run voronoi
        v = VoronoiWrapper(slice_atoms, padding)
        v.generate_voronoi()

        # generate the areas for the selections
        for i in range(len(selections[1:])):
            s = selections[i+1]
            for j in range(len(s)):
                sr = SuperRegion()
                gr = slicer(s[j])
                if (len(gr) == 0):  # skip if there are no atoms selected from mol
                    continue
                sr.buildFromAtoms(gr, v)
                a = sr.area()
                index = int((a-min_area) / bin_width)
                try:
                    histograms[i][index] += 1
                except IndexError:
                    print("#Area outside range:  ", pytraj.realIndex(), i, j,
                          a, index)

    # normalize the histograms
    for i in range(len(histograms)):
        histograms[i] /= numpy.add.reduce(histograms[i])

    # final output
    string = ""
    for i in range(len(selections[1:])):
        string += "\tSel" + str(i)
    print("# Area", string)
    for i in range(num_bins):
        a = min_area + (i+0.5)*bin_width
        string = list(map(str, histograms[:, i]))
        print(a, "\t".join(string))

if __name__ == '__main__':
    main()