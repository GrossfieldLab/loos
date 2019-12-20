#!/usr/bin/env python3

import loos
import loos.pyloos
import numpy
from sklearn.manifold import TSNE
import sys

command_line = " ".join(sys.argv)
system_file = sys.argv[1]
traj_file = sys.argv[2]
selection_string = sys.argv[3]
output_file = sys.argv[4]
perplexity = float(sys.argv[5])

system = loos.createSystem(system_file)
traj = loos.pyloos.Trajectory(traj_file, system)

sel = loos.selectAtoms(system, selection_string)

coord_list = []
for frame in traj:
    v = sel.getCoords()
    v = v.reshape([v.shape[0]*v.shape[1]])
    coord_list.append(v)

coord_list = numpy.array(coord_list)

# do the TSNE
t = TSNE(perplexity=perplexity,
         n_components=2,
         init="pca")
embedded = t.fit_transform(coord_list)
print(t.get_params(deep=True))

numpy.savetxt(output_file, embedded, header=command_line)
