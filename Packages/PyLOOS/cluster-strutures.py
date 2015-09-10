#!/usr/bin/env python
"""
Cluster structures based from a simulation
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
from scipy.cluster.vq import kmeans,vq
from itertools import chain
import copy


verbose = 0

if (len(sys.argv) <= 5) or (sys.argv[1] == "-h"):
    print "Usage: cluster-structures.py model selection num_means prefix traj [traj...] > output-for-metadata"
    sys.exit(1)



## Command line arguments
model_name = sys.argv[1]
selection = sys.argv[2]
num_means = sys.argv[3]
outfile =  sys.argv[4]
cmd_string = sys.argv[0] 
for i in range(1, len(sys.argv)):
    arg = sys.argv[i].replace('\n', '\\n')
    cmd_string += " '" + arg + "'"
print '# ', cmd_string


trajList = sys.argv[5:len(sys.argv)]

# Create the model & read in the trajectories
model = loos.createSystem(model_name)
allTrajs = loos.pyloos.VirtualTrajectory()
for trajname in trajList:
    allTrajs.append(loos.pyloos.Trajectory(trajname, model, subset=selection))


# Set up lists to hold the coordinates...

n = len(allTrajs.currentFrame()) * 3
data = numpy.zeros( (len(allTrajs), n) )
for frame in allTrajs:
    coords = frame.getCoords()
    obs = numpy.reshape(coords, (1, n))
    data[allTrajs.currentIndex()] = obs


# Do the clustering...
# Computing K-Means with K = num_means clusters
centroids,distortion = kmeans(data, int(num_means))
# centroids  - the codebook of centroids 
# distortion - total distortion

# Assign each sample to a cluster
idx,dists = vq(data, centroids)
# idx   - code (which cluster a point belongs to)
# dists - distance of each point from cluster 
#         used for the distortion calculation
#         May want to output for determining 
#         number of clusters in system          


# Write out the meta-data file
print "# Means\tDistortion: "
print num_means," \t",distortion
print "# -------------------------------"
print "# Trajectory list:"
for i in range(len(trajList)):
    print '# %5d = "%s"' % (i, trajList[i])

print '#\n# %8s %16s %8s %8s' % ('Index', 'Trajectory', 'Frame', 'Cluster')
print '# %8s %16s %8s %8s' % ('--------', '----------------', '--------', '--------')

for i in range(len(idx)):
    loc = allTrajs.frameLocation(i)
    print '%10d %16d %8d %8d' % (i, loc[1], loc[3], idx[i])

    
# Output centroids
cen_list = centroids.tolist()
subset = allTrajs.currentFrame()
for j in range(len(cen_list)):
    troid = cen_list[j]
    centroid_structure = subset.copy()
    for i in range(0, len(troid), 3):
        centroid_structure[i/3].coords(loos.GCoord(troid[i], troid[i+1], troid[i+2]))
    pdb = loos.PDB.fromAtomicGroup(centroid_structure)
    pdb.remarks().add(cmd_string)
    pdb.remarks().add(">>> Means = %s, Distortion = %f" % (num_means, distortion))

    filename = "%s-centroid-%d.pdb" % (outfile, j)
    file = open(filename, 'w')
    file.write(str(pdb))
    file.close()
