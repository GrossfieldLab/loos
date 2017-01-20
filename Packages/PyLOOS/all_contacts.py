#!/usr/bin/python

import loos
import loos.pyloos
import sys
import numpy
from os.path import basename, splitext

if len(sys.argv) < 5 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
    print "Usage: ", sys.argv[0], "system selection out_file trajectory_file [trajectory2 ...]"
    print """
        Description: all_contacts.py: tracks contacts between pairs of side
                     chains in a protein or nucleic acid.  Returns a matlab-format
                     matrix containing the fraction of frames containing each
                     residue-residue contact, as well as a matrix for each
                     trajectory individually.

                     The threshold for deciding if 2 residues are in contact is
                     whether any atoms from one residue are within 4.0 ang of
                     any atom from the other.  NOTE: PERIODICITY IS NOT RESPECTED,
                     since we're assuming they're both in the same molecule.
          """

system_file = sys.argv[1]
selection = sys.argv[2]
out_file = sys.argv[3]
traj_files = sys.argv[4:]

header =  " ".join(sys.argv) + "\n"


system = loos.createSystem(system_file)
all_trajs = []
out_names = []
num_trajs = len(traj_files)
for t in traj_files:
    traj = loos.pyloos.Trajectory(t, system)
    all_trajs.append(traj)
    t_base = basename(t)
    core,ext = splitext(t_base)
    out_names.append(core + ".dat")

# Used to remove hydrogens, but that messes up glycines
#no_hydrogens = loos.selectAtoms(system, "!hydrogen")
#target = loos.selectAtoms(no_hydrogens, selection)

target = loos.selectAtoms(system, selection)


residues = target.splitByResidue()
# now remove the backbone -- doing before the split loses the glycines
# Woohoo, look at me, I used a lambda!
residues = list(map(lambda r:loos.selectAtoms(r, "!backbone"), residues))

frac_contacts = numpy.zeros([len(residues), len(residues), num_trajs],
                            numpy.float)


for traj_id in range(num_trajs):
    traj = all_trajs[traj_id]
    for frame in traj:
        for i in range(len(residues)):
            for j in range(i+1, len(residues)):
                if residues[i].contactWith(4.0, residues[j]):
                    frac_contacts[i,j,traj_id] += 1.0
                    frac_contacts[j,i,traj_id] += 1.0
    frac_contacts[:,:,traj_id] /= len(traj)
    numpy.savetxt(out_names[traj_id], frac_contacts[:,:,traj_id],
                  header=header)



average = numpy.add.reduce(frac_contacts, axis=2)
average /= len(traj_files)


numpy.savetxt(out_file, average, header = header)
