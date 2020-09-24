
import loos
from loos import pyloos
import json
from sys import argv

usage = 'python frame-picker.py results.json start:step:stop prefix model traj [traj] ...'
for i, arg in enumerate(argv):
    print(i, arg)


if len(argv) == 1 or '-h' in argv:
    print(usage)
    exit(0)


cluster_results_fn = argv[1]
# parse clustering output here
with open(cluster_results_fn,'r') as fp:
    cluster_results = json.load(fp)
# trajectory stuff
prerange = argv[2].split(':')
if len(prerange) != 3:
    print('Stop feeding me weird things!')
    print(trajranges)
    print(usage)
    exit(0)
if prerange[-1]:
    trajranges = tuple(map(int, prerange))
    print('Warning: code does not currently handle "stop"')
else:
    print('Warning: code does not currently handle "stop"')
    trajranges = (int(prerange[0]), int(prerange[1]), -1)
skip, stride, stop = trajranges
prefix = argv[3]
model = loos.createSystem(argv[4])
vtraj = pyloos.Trajectory(argv[5], model, skip=skip, stride=stride)
# vtraj = pyloos.AlignedVirtualTrajectory(model, alignwith='!hydrogen')
# for t in argv[5:]:
#     traj = pyloos.Trajectory(t, model)[skip:stop:stride]
#     print(len(traj))
#     vtraj.append(traj)
total = len(vtraj)
print(total)

for cluster_ix, cluster in enumerate(cluster_results['clusters']):
    # puts the fraction of the trajectory lumped into ith traj in the name
    frac = round(len(cluster)/float(total), 4)
    outprefix = prefix+'-'+str(frac)+'-'+str(cluster_ix)
    cluster_traj = loos.DCDWriter(outprefix+'.dcd')
    exemplar_ix = cluster_results['exemplars'][cluster_ix]
    frame = vtraj[exemplar_ix]
    exemplar = loos.PDB.fromAtomicGroup(frame)
    with open(outprefix+'.pdb', 'w') as fp:
        fp.write(str(exemplar))

    # OK now write out the trajectory.
    for member_idx in cluster:
        member = vtraj[member_idx]
        cluster_traj.writeFrame(member)
        


