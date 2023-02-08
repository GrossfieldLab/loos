#!/usr/bin/env python3
"""
Use Voronoi polygons to compute the lifetime of lipids near the surface of
a protein


Alan Grossfield, University of Rochester Medical Center
Copyright 2015
"""

import sys
import loos
import loos.pyloos
import numpy

from loos.Voronoi import *


def autocorrel(timeseries, expected=None):
    """
    Computes the autocorrelation function for a set of timeseries.

    Destroys the timeseries matrix.
    Does not trap for timeseries with zero std dev (e.g. all constant)
    """
    ave = numpy.mean(timeseries, axis=1, keepdims=True)
    dev = numpy.std(timeseries, axis=1, keepdims=True)
    timeseries -= ave
    timeseries /= dev

    length = timeseries.shape[1]
    if expected is None:
        expected = length

    corr = numpy.zeros([expected, len(timeseries)], float)
    for i in range(expected):
        num_vals = length - i
        num = timeseries[:, i:] * timeseries[:, :num_vals]
        corr[i] = numpy.average(num, axis=1)

    ave_corr = numpy.average(corr, axis=1)
    dev_corr = numpy.std(corr, axis=1)

    return ave_corr, dev_corr


def Usage():
    return """
Usage: lipid_lifetime.py system trajectory zmin zmax protein-selection all-lipid-selection target-lipids

This program uses 2D Voronoi decomposition to identify lipids in contact
with the protein surface, and returns the time-correlation function
for this contact (averaged over all lipids considered).  This can be fit
to a sum of exponentials to extract a "lifetime".

zmin and zmax define the "slice" of the membrane and protein used in the
Voronoi decomposition.  The details of the choice shouldn't matter too much,
as long as you pick out one leaflet or the other; beyond that, your choice
should be guided by the shape of the protein and the variation of whatever atom
you use to define the lipid.  At the moment, a lipid that leaves the z-slice
is treated the same as one that leaves the protein, so I'd say handling of
flip-flop is not very robust.

protein-selection is the set of atoms used to define the "protein", or the
thing that you're calculating the lifetime around.  It is assumed to be 1 thing
(e.g. all the atoms are treated together, and we actually just use the
centroid), so you couldn't put a selection for something like a bunch of
cholesterol molecules there and expect to get the lifetime for lipids around
cholesterol.  A fairly sparse selection is probably best, e.g. just alpha
carbons.

all-lipids is a selection that described all of the lipids that are found in
the leaflet you're examining.  The code is written assuming that you'll pick
one atom per lipid, e.g. the phosphorus.  The important thing is that all
lipids (and cholesterols, etc) in the leaflet in question must be represented,
or else the resulting "gaps" may cause some lipids to be incorrectly be
marked as neigboring the protein.

target-lipids is a selection of those lipids for which you want to compute the
lifetime.  This group is by definition a subset of the all-lipids selection;
since we actually apply this selection to the group created by the all-lipids
selection, you can use 'all' as this selection if you want to compute the
correlation function for every lipid.  Alternatively, if you have a mixture
of POPE and POPC lipids, but you want just the lifetime for just the POPE
lipids, you might use 'resname =~ "POP[CE]" and name == "P"' for all-lipids
and 'resname == "POPE"' for target-lipids.


Note: The third column in the output is the standard deviation in the
correlation function, where you've averaged all lipids.  If you believe
that the lipids can be treated as independent measurements, you would divide
this value by sqrt(N_lipids_considered) to get the standard error in the
correlation function, which is related to the statistical uncertainty.  I
didn't do this because it's not clear to me that that is a good assumption.
    """


def main():
    if len(sys.argv) <= 1 or sys.argv[1] == '--fullhelp':
        print(Usage())
        sys.exit()

    print("#", " ".join(sys.argv))

    system_file = sys.argv[1]
    traj_file = sys.argv[2]
    zmin = float(sys.argv[3])
    zmax = float(sys.argv[4])
    protein_selection = sys.argv[5]
    all_lipid_selection = sys.argv[6]
    target_lipid_selection = sys.argv[7]

    system = loos.createSystem(system_file)
    # TODO: go back and support skip/stride later
    traj = loos.pyloos.Trajectory(traj_file, system)

    protein = loos.selectAtoms(system, protein_selection)
    all_lipids = loos.selectAtoms(system, all_lipid_selection)
    # NOTE: target_lipids must be a subset of all_lipids
    target_lipids = loos.selectAtoms(all_lipids, target_lipid_selection)

    slicer = ZSliceSelector(zmin, zmax)

    # TODO: this should probably be a command line option
    # Note: normally, this small a padding might be a problem for Voronoi
    #       decomposition with 1 atom/lipid. However, we're not using the areas
    #       (which is what gets screwed up by lack of padding), and 15 ought
    #       to be big enough to make sure we've got 1 layer of lipid around
    #       the protein.
    padding = 15.

    protein_centroid = loos.AtomicGroup()
    protein_centroid.append(loos.Atom())

    # set up space to hold the neigbor time series
    neighbor_timeseries = numpy.zeros([len(target_lipids), len(traj)], float)

    for frame in traj:
        # Use the centroid of the protein slice to represent the protein
        protein_slice = slicer(protein)
        centroid = protein_slice.centroid()
        protein_centroid[0].coords(centroid)

        # We assume you're using 1 atom/lipid.
        # Applying slice operation here allows you to be sloppy
        # with your selections.
        all_lipids_slice = slicer(all_lipids)
        target_lipid_slice = slicer(target_lipids)

        combined = protein_centroid + all_lipids_slice

        # translate so that the protein is in the center, then reimage by
        # atom
        combined.translate(-protein_centroid[0].coords())
        combined.periodicBox(frame.periodicBox())
        combined.reimageByAtom()

        # run voronoi
        v = VoronoiWrapper(combined, padding)
        v.generate_voronoi()

        protein_region = v.get_region_from_atomid(protein_centroid[0].id())
        # loop over all target lipids
        for i in range(len(target_lipid_slice)):
            region = v.get_region_from_atomid(target_lipid_slice[i].id())
            if protein_region.is_neighbor(region):
                neighbor_timeseries[i][traj.index()] = 1

    # remove entries that are all-zero.
    ave = numpy.mean(neighbor_timeseries, axis=1)
    is_nonzero = numpy.select([abs(ave) > 1e-6], [ave])
    neighbor_timeseries = numpy.compress(is_nonzero,
                                         neighbor_timeseries, axis=0)

    ave_corr, dev_corr = autocorrel(neighbor_timeseries, neighbor_timeseries.shape[1]//2)

    print("#Frame\tCorrel\tDev")
    for i in range(len(ave_corr)):
        print(i, ave_corr[i], dev_corr[i])

if __name__ == '__main__':
    main()