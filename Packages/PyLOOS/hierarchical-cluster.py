#!/usr/bin/env python
"""
Given a distance matrix, compute hierarchical clustering of simulation 
structures.

Alan Grossfield, University of Rochester Medical Center
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
import numpy
import scipy
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from os.path import basename, splitext
import argparse

def fullhelp():
    print """
    Perform hierarchical clustering given a matrix of distances

    Summary of Arguments
      rmsd_file :     matlab-format text file containing "distances" between all
                      observations.  For example, this could be the output of 
                      the LOOS program rmsds
      num_clusters:   number of clusters to truncate the clustering

    Optional arguments
      --index_file    text file containing a list of filenames and which rows 
                      of the rmsd_file they correspond to.
                      The format is: 
                      START_FRAME    END_FRAME    FILENAME
                      using 0-based indexing.  This could be useful if you 
                      passed multiple trajectories to rmsds, but now want to 
                      pull out the cluster populations for each individual 
                      trajectory.
    --prefix          name to use for the outputted time series and histogram.  
                      --prefix "foo" produces "foo.dat" and "foo.hist"
    --link            writes out the linkage matrix to "all.link" or the 
                      equivalent name if you supplied a --prefix.

    Output
    If you did not pass in an index file or a prefix, this program will produce 
    2 output files, "all.dat" and "all.hist".  "all.dat" contains the time 
    series of cluster assignments, while "all.hist" is a normalized histogram of 
    cluster populations.  If you supplied a "--prefix" argument, that string 
    will replace "all" in the .dat and .hist filenames.

    If you did supply an index file, it will additionally supply a .dat and 
    .hist for each file listed in the index file.  It will strip off any 
    leading path and anything that looks like a filename extension, so that
    "/foo/bar/baz.dcd" would produce "baz.dat" and "baz.hist".

    Technical Details

    This program is mostly a thin wrapper around scipy's hierarchical clustering
    routines.  It uses the UPGMA algorithm: when trying to decide whether to 
    merge 2 clusters, it takes the average of the distances between all pairs of 
    points in the 2 clusters.  Because we supply a distance matrix instead of 
    the actual data points (which would be much bigger, and would make this 
    program much more complicated), we can't use some of the other algorithms, 
    like Ward's algorithm or the median or centroid algorithms.  

    """

cmd_args = " ".join(sys.argv)

parser = argparse.ArgumentParser()
parser.add_argument('rmsd_file', help="File containing a distance matrix in matlab format, e.g. output from the LOOS program rmsds")
parser.add_argument('num_clusters', help='Number of clusters', type = int) 

parser.add_argument('--index_file', help="File identifying which rows come from which trajectory", default=None, type = str)
parser.add_argument('--prefix', help="Core of the output filenames", default = 'all', type =str)
parser.add_argument('--fullhelp', action='store_true')
parser.add_argument('--link', action = 'store_true')

args = parser.parse_args()
if args.fullhelp:
    print fullhelp()
    sys.exit(1)
print "#", cmd_args

num_clusters = args.num_clusters
if args.index_file is not None:
    index_filename = args.index_file
else:
    index_filename = None

# read the RMSD file and
# convert to condensed upper triangular
rmsd = numpy.loadtxt(args.rmsd_file)
upper = squareform(rmsd)


link = sch.average(upper)

if args.link:
    numpy.savetxt(args.prefix + ".link", link, header=cmd_args)

assignments = sch.fcluster(link, num_clusters, criterion='maxclust')

# Read the index file if one was supplied
indices = {}
if index_filename:
    with(open(index_filename)) as index_file:
        for line in index_file.readlines():
            (first, last, filename) = line.split()
            first = int(first)
            last = int(last)
            t = basename(filename)
            sim, ext = splitext(t)
            indices[sim] = (first, last)

# always produce a histogram of all of the data
indices[args.prefix] = (0, len(rmsd)-1)

for sim in indices.keys():
    filename = sim + "_" + str(num_clusters) + ".dat"
    histname = sim + "_" + str(num_clusters) + ".hist"
    first, last = indices[sim]
    sim_slice = assignments[first:last+1]
    numpy.savetxt(filename, sim_slice.astype(int), fmt="%d", header=cmd_args)
    sim_hist = numpy.bincount(sim_slice, minlength=num_clusters+1)
    sim_hist = sim_hist.astype(float) 
    sim_hist /= numpy.add.reduce(sim_hist)
    numpy.savetxt(histname, sim_hist, header=cmd_args)
