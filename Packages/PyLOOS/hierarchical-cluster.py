#!/usr/bin/env python
#
# See https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/


import sys
import numpy
import scipy
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from os.path import basename, splitext

#from matplotlib import pyplot as plt


rmsd = numpy.loadtxt(sys.argv[1])
distance = int(sys.argv[2])
index_filename = sys.argv[3]

# convert rmsd to condensed upper triangular
#indices = numpy.triu_indices_from(rmsd)
#upper = numpy.asarray( rmsd[indices] )
upper = squareform(rmsd)


#link = scipy.cluster.hierarchy.centroid(rmsd)
link = sch.average(upper)

#print 'Saving linkage to link.asc'
numpy.savetxt('link.asc', link)
#c, coph_dists = sch.cophenet(link, upper)
#print c


#plt.figure()
#plt.title('crap')
#plt.xlabel('Cluster Size')
#plt.ylabel('Distance')
#sch.dendrogram(link, leaf_rotation=90., leaf_font_size=12,
#                                   truncate_mode='lastp', p=25, show_contracted=True)
#
#plt.show()

#assignments = sch.fcluster(link, distance, criterion='distance')
assignments = sch.fcluster(link, distance, criterion='maxclust')
numpy.savetxt('clusters.asc', assignments.astype(int), fmt="%d")

# split up the assignments by trajectory name
indices = {}
index_file = open(index_filename)
for line in index_file.readlines():
    (first, last, filename) = line.split()
    first = int(first)
    last = int(last)
    t = basename(filename)
    sim, ext = splitext(t)
    indices[sim] = (first, last)

stim = numpy.zeros(distance+1)
stim2 = numpy.zeros(distance+1)
phos = numpy.zeros(distance+1)
phos2 = numpy.zeros(distance+1)
t389e = numpy.zeros(distance+1)
t389e2 = numpy.zeros(distance+1)

stim_count = 0
phos_count = 0
t389e_count =0 

for sim in indices.keys():
    filename = sim + "_" + str(distance) + ".dat"
    histname = sim + "_" + str(distance) + ".hist"
    first, last = indices[sim]
    sim_slice = assignments[first:last+1]
    numpy.savetxt(filename, sim_slice.astype(int), fmt="%d")
    sim_hist = numpy.bincount(sim_slice, minlength=distance+1)
    sim_hist = sim_hist.astype(float) 
    sim_hist /= numpy.add.reduce(sim_hist)
    numpy.savetxt(histname, sim_hist)

    if sim.startswith("stim"):
        stim += sim_hist
        stim2 += sim_hist*sim_hist
        stim_count += 1
    elif sim.startswith("phos"):
        phos += sim_hist
        phos2 += sim_hist*sim_hist
        phos_count += 1
    elif sim.startswith("t389e"):
        t389e += sim_hist
        t389e2 += sim_hist*sim_hist
        t389e_count += 1

stim /= stim_count
stim2 /= stim_count
stim_err = numpy.sqrt((stim2 - stim*stim)/stim_count)
stim_err =numpy.nan_to_num(stim_err)

phos /= phos_count
phos2 /= phos_count
phos_err = numpy.sqrt((phos2 - phos*phos)/phos_count)
phos_err =numpy.nan_to_num(phos_err)

t389e /= t389e_count
t389e2 /= t389e_count
t389e_err = numpy.sqrt((t389e2 - t389e*t389e)/t389e_count)
t389e_err =numpy.nan_to_num(t389e_err)
    
string = "Files Stim Stim Phos Phos T389E T389E"
print string

for i in range(distance+1):
    string = "%d\t" %(i)
    string += "%f\t%f\t" % (stim[i], stim_err[i])
    string += "%f\t%f\t" % (phos[i], phos_err[i])
    string += "%f\t%f\t" % (t389e[i], t389e_err[i])
    print string


# Alternatively, use maxd and distance as criterion
#maxk = input('How many clusters should I use? ')
#clusters = scipy.cluster.hierarchy.fcluster(link, maxk, criterion='maxclust')
#print 'Frame to cluster ID assignments saved to clusters.asc'
#numpy.savetxt('clusters.asc', clusters)


