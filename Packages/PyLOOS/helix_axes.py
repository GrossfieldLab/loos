#!/usr/bin/env python3

import loos
import loos.pyloos
import sys
import math

if len(sys.argv) < 5 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
    print("Usage: ", sys.argv[0], " system trajectory out_prefix selection1 [selection2...]")
    print("""
    Description: helix_axes.py: Writes a new trajectory replacing each 
                 selection with 3 particles, one at the centroid of the
                 selection and two stepped 10 angstroms up and down the 
                 first principal axis from it.  It's a quick way to 
                 simply view a helix's motion, and useful as input for
                 a simplified PCA.  
    Command line arguments:
                 system and trajectory have their usual meanings.
                 output_prefix is used to name the resulting PDB file
                 and trajectory (e.g. "foo" would produce "foo.pdb" and 
                 "foo.dcd".  
                 You need to supply at least one selection string, but 
                 you can supply as many as you want (e.g. for rhodopsin
                 it might make sense to give 7 selections, one for each
                 transmembrane helix).  The resulting file will have
                 3x (number of selections) atoms.
    Despite the name of the program, the selections do not have to be 
    helical.  However, I don't know that the approach makes much sense 
    for anything else.
          """)
    sys.exit()

header = " ".join(sys.argv)
system_filename = sys.argv[1]
traj_filename = sys.argv[2]
out_prefix = sys.argv[3]
selections = sys.argv[4:]


system = loos.createSystem(system_filename)
traj = loos.pyloos.Trajectory(traj_filename, system)


helices = []
for s in selections:
    helices.append(loos.selectAtoms(system, s))

# Make a fake atomicgroup to hold the placeholder atoms
new_group = loos.AtomicGroup()
resnum = 1
for i in range(len(helices)):
    a = loos.Atom()
    b = loos.Atom()
    c = loos.Atom()
    a.resid(resnum)
    b.resid(resnum)
    c.resid(resnum)
    a.name("CENT")
    b.name("PLUS")
    c.name("MIN")

    new_group.append(a)
    new_group.append(b)
    new_group.append(c)

    resnum+=1 

new_group.renumber()

first = True
for frame in traj:

    vec = loos.GCoord(0.,0.,0.)
    for i in range(len(helices)):
        h = helices[i]
        pca = h.principalAxes()
        v = pca[0]
        if v.z() < 0:
            v *= -1.

        centroid = h.centroid()
        plus = centroid + 10.*v
        minus = centroid - 10.*v

        helix_atoms = loos.selectAtoms(new_group, 'resid == ' + str(i+1))
        helix_atoms[0].coords(centroid)
        helix_atoms[1].coords(plus)
        helix_atoms[2].coords(minus)

    # set up output trajectory on the first frame
    if first:
        first = False
        out_pdb = loos.PDB_fromAtomicGroup(new_group)
        pdb_name = out_prefix + ".pdb"
        with open(pdb_name, 'w') as pdb_file:
            pdb_file.write(str(out_pdb))
    
        out_dcd = out_prefix + ".dcd"
        dcd_writer = loos.DCDWriter(out_dcd)
        dcd_writer.setTitle(header)


    dcd_writer.writeFrame(new_group)

