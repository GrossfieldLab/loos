#!/usr/bin/env python3
"""
Cluster structures based from a simulation, using the K-means
algorithm applied to RMSD
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
import argparse
from scipy.cluster.vq import kmeans, vq
from sklearn.manifold import TSNE




parser = argparse.ArgumentParser()
parser.add_argument('model', help='Structure to use')
parser.add_argument('selection', help='Selection to use for clustering')
parser.add_argument('num_means', help='# of clusters to make', type=int)
parser.add_argument('prefix', help='Prefix output files with this')
parser.add_argument('traj', help='Trajectory to use', nargs='+')
parser.add_argument('--align',
                    help='Align trajectory using this LOOS selection')
parser.add_argument('--allskip',
                    help='Skip this amount from the start of each trajectory',
                    type=int,
                    default=0)
parser.add_argument('--allstride',
                    help='Step through each trajectory by this many frames',
                    type=int,
                    default=1)
parser.add_argument('--skip',
                    help='Skip frames the start of the combined trajectory',
                    type=int,
                    default=0)
parser.add_argument('--stride', help='Do every nth frame',
                    type=int,
                    default=1)
parser.add_argument('--tsne',
                    help="Perform tSNE with perplexity",
                    type=bool,
                    default=False)
parser.add_argument('--perplexity',
                    help="Set perplexity for tSNE",
                    type=float,
                    default=30.
                    )
args = parser.parse_args()
cmd_string = sys.argv[0]
for i in range(1, len(sys.argv)):
    arg = sys.argv[i].replace('\n', '\\n')
    cmd_string += " '" + arg + "'"
print('# ', cmd_string)


# Create the model & read in the trajectories
model = loos.createSystem(args.model)
allTrajs = loos.pyloos.VirtualTrajectory(skip=args.skip, stride=args.stride)
if args.align:
    allTrajs = loos.pyloos.AlignedVirtualTrajectory(alignwith=args.align,
                                                    skip=args.skip,
                                                    stride=args.stride)
for trajname in args.traj:
    allTrajs.append(loos.pyloos.Trajectory(trajname,
                                           model,
                                           subset=args.selection,
                                           skip=args.allskip,
                                           stride=args.allstride))


# Set up lists to hold the coordinates...

n = len(allTrajs.frame()) * 3
data = numpy.zeros((len(allTrajs), n))
for frame in allTrajs:
    coords = frame.getCoords()
    obs = numpy.reshape(coords, (1, n))
    data[allTrajs.index()] = obs


if args.align:
    print('# Iteratively aligned with %d iterations and final RMSD %g.' %
          (allTrajs._iters, allTrajs._rmsd))

if args.tsne:
    print("# Performing tSNE on aligned trajectories")
    t = TSNE(perplexity=args.perplexity,
             n_components=2,
             init="pca")
    embedded = t.fit_transform(data)


# Do the clustering...
# Computing K-Means with K = num_means clusters
# and assign each value to a cluster
# centroids  - the codebook of centroids
# distortion - total distortion
# Assign each sample to a cluster
# idx   - code (which cluster a point belongs to)
# dists - distance of each point from cluster
#         used for the distortion calculation
#         May want to output for determining
#         number of clusters in system
if not args.tsne:
    centroids, distortion = kmeans(data, args.num_means)
    idx, dists = vq(data, centroids)
    subset = allTrajs.frame()
    dists *= 1.0
    dists /= len(subset)
else:
    centroids, distortion = kmeans(embedded, args.num_means)
    idx, dists = vq(embedded, centroids)


# Write out the meta-data file
print("# Means\tDistortion: ")
print('# ', args.num_means, " \t", distortion)
print("# -------------------------------")
print("# Trajectory list:")
for i in range(len(args.traj)):
    print('# %5d = "%s"' % (i, args.traj[i]))

print('#\n# %8s %16s %8s %8s %8s' %
      ('Index', 'Trajectory', 'Frame', 'Cluster', 'Distance'))
print('# %8s %16s %8s %8s %8s' %
      ('--------', '----------------', '--------', '--------', '--------'))

minima = numpy.full_like(numpy.arange(args.num_means),
                         1.e10,
                         dtype=numpy.float)
minima_indices = numpy.zeros([args.num_means], dtype=numpy.int)

for i in range(len(idx)):
    loc = allTrajs.frameLocation(i)
    if (dists[i] < minima[idx[i]]):
        minima[idx[i]] = dists[i]
        minima_indices[idx[i]] = i
    print('%10d %16d %8d %8d %8f' % (i, loc[1], loc[3], idx[i], dists[i]))

print("\n#  Closest structure to each centroid")
print('# %8s %10s     %8s' % ("Cluster", "Index", "Distance"))
print('# %8s %10s     %8s' % ("-------", "-----", "--------"))
for i in range(args.num_means):
    print("# %8d %10d     %8f" % (i, minima_indices[i], minima[i]))


# Output centroids
cen_list = centroids.tolist()
subset = allTrajs.frame()
for j in range(len(cen_list)):
    troid = cen_list[j]
    centroid_structure = subset.copy()
    for i in range(0, len(troid), 3):
        centroid_structure[i//3].coords(loos.GCoord(troid[i],
                                                    troid[i+1],
                                                    troid[i+2]))
    pdb = loos.PDB.fromAtomicGroup(centroid_structure)
    pdb.remarks().add(cmd_string)
    pdb.remarks().add(">>> Means = %s, Distortion = %f" %
                      (args.num_means, distortion))

    filename = "%s-centroid-%d.pdb" % (args.prefix, j)
    file = open(filename, 'w')
    file.write(str(pdb))
    file.close()
