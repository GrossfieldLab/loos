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
from loos import *
from numpy import vstack,array
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


trajList = []
for myargs in range(5, len(sys.argv)):
    trajList.append(sys.argv[myargs])


# Create the model & read in the trajectories
model = createSystem(model_name)
allTrajs = []
for trajs in range( len(trajList) ):
    allTrajs.append( createTrajectory(trajList[trajs], model) )
    if (verbose):
        print "# Trajectory", trajs, "contains", allTrajs[trajs].nframes(), "frames"

# Create a list of coords from each frame
fullTrajs = []
subset = selectAtoms(model, selection)
ensemble = AtomicGroupVector()
dataSet = 0 
for trajs in allTrajs:
    fullTrajs.append(readTrajectory(ensemble, subset, trajs))

print "# Subset has %d atoms." % (len(subset))
print "# Total of", len(ensemble), "structures."

# Set up lists to hold the coordinates...
thisElement = []
allElements = []

# Extract coords 
for eachElement in ensemble:
    for everyAtom in eachElement:
        thisElement.append(everyAtom.coords().x())
        thisElement.append(everyAtom.coords().y())
        thisElement.append(everyAtom.coords().z())
    allElements.append(thisElement)
    thisElement = []

# Format for numpy and scipy
fullArray = array(allElements)
dataStack = vstack((fullArray))

# Do the clustering...
# Computing K-Means with K = num_means clusters
centroids,distortion = kmeans(dataStack, int(num_means))
# centroids  - the codebook of centroids 
# distortion - total distortion

# Assign each sample to a cluster
idx,dists = vq(dataStack, centroids)
# idx   - code (which cluster a point belongs to)
# dists - distance of each point from cluster 
#         used for the distortion calculation
#         May want to output for determining 
#         number of clusters in system          


# Write out the meta-data file
print "# Means\tDistortion: "
print num_means," \t",distortion
print "# -------------------------------"
print "# index\t trajectory"
for item in range(5, len(sys.argv)):
    print "#", item-5,"\t ", sys.argv[item]

# Output centroids
cen_list = centroids.tolist()
for j in range(len(cen_list)):
    troid = cen_list[j]
    centroid_structure = subset.copy()
    for i in range(0, len(troid), 3):
        centroid_structure[i/3].coords(GCoord(troid[i], troid[i+1], troid[i+2]))
    pdb = PDB.fromAtomicGroup(centroid_structure)
    pdb.remarks().add(cmd_string)
    pdb.remarks().add(">>> Means = %s, Distortion = %f" % (num_means, distortion))

    filename = "%s-centroid-%d.pdb" % (outfile, j)
    file = open(filename, 'w')
    file.write(str(pdb))
    file.close()


    
# Output data to files named prefix.[traj-index]
currentTraj = 0;
counter = 0
current_outfile = outfile+"."+str(currentTraj)
filewr = open(current_outfile, 'w')
filewr.write("#Index Frame Cluster Trajectory\n")
filewr.write("#------------------------------\n")

for index in range(0, len(idx)):
    if (counter >= allTrajs[currentTraj].nframes() ):
        dataSet += 1
        currentTraj += 1
        counter = 0
        current_outfile = outfile+"."+str(currentTraj)
        filewr = open(current_outfile, 'w')
        filewr.write("#Index Frame Cluster Trajectory\n")
        filewr.write("#------------------------------\n")
    filewr.write( str(index) )
    filewr.write(" ")
    filewr.write( str(counter) )
    filewr.write(" ")
    filewr.write( str(idx[index]) )
    filewr.write(" ")
    filewr.write(str(dataSet))
    filewr.write("\n")
    counter += 1

