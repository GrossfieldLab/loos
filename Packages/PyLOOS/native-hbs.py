#!/usr/bin/env python

import loos
import loos.pyloos
import numpy as np
import argparse

"""
    native-hbs.py : find the hydrogen bonds present in some reference
                    selection, then output their average occupancy
                    over a trajectory (and optionally a timeseries of 
                    those occupancies and a PDB of the atoms in the HBs)

    The expected use-case is the commonly asked question, 'does this Hydrogen
    bonded structure persist in the simulation?' The expectation is that
    the user would have some molecular model (such as a PDB) where the 
    H-bonded structure is present (say, a Watson-Crick duplex) and a 
    trajectory of that molecule which may no longer retain the structure.

    You may optionally specify files to get not just the average HB occupancy
    in your selection, but also the per-frame occupancy of each hydrogen bond
    in the selection in an output file. You can also optionally provide an 
    outfile name for a PDB of the hydrogen bonded atoms the script finds.
    This makes for easy visualisation: load that file and your ref up in a 
    visualizer, then view the model in sticks (or something smaller) and the
    H-bonded atoms as spheres (or something bigger). Alternatively, you could
    have your structure viewer draw HB for one 'molecule' then the other
    and see whether the script found the ones you were hoping it would

    Example commandline with no optional arguments:
    native-hbs.py -m xtal_molA.pdb -x traj_molA.dcd

    This will print the following values separated by spaces to stdout:
    Number_frames Number_HBs Average_HBs

    Example commandline with all the optional arguments:
    native-hbs.py -m xtal_mol_A.pdb -x traj_mol_A.dcd \
    -n native_hbs_found_A.pdb -o native_hbs_timeseries_A.dat \
    -b 4 -a 20 -s (resid < 4 || resid > 8) && !backbone'

    Note that with all the options present the three values above are
    still written to standard out, but the outfile specified by -o will
    contain the occupancy of each found HB by frame so they can be readily
    recalculated from those files if you chose not to redirect the output.

    Louis G. Smith

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
# takes list of DHA triples, obtains their atomnames, 
def format_header(dhas):
    retlist = ['# frame']
    for x in dhas:
        names = ['_'.join([i.resname(), str(i.resid()), i.name()]) for i in x]
        retlist.append('-'.join(names))
    retlist.append('total_inframe')
    return '\t'.join(retlist) + '\n'

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model', help='filename of Structure to use')
parser.add_argument('-c', '--coords', help='filename of reference coordinates, if they are not contained in supplied model', default='')
parser.add_argument('-x', '--traj', help='Filename of trajectory or trajectories', nargs='+')
parser.add_argument('-k', '--skip', help='Skip this amount from the start of the trajectory', type=int, default=0)
parser.add_argument('-i', '--stride', help='Step through the trajectory by this many frames', type=int, default=1)
parser.add_argument('-s', '--selection', help='Use this selection for computing the native hbonds', default='all')
parser.add_argument('-n', '--nativeHBs', help='a filename to write PDB with just the native HBs used', type=str, default='')
parser.add_argument('-b', '--bondlength', help='maximum length for HB.', type=float, default = 3.0)
parser.add_argument('-o', '--outfile', help='file to write timeseries to.', type=str, default='')
parser.add_argument('-a', '--angle', help='maximum valid D-H-A angle.', type=float, default = 30.0)
args = parser.parse_args()

# takes list of DHA triples, obtains their atomnames, 
def format_header(dhas):
    retlist = ['# frame']
    for x in dhas:
        names = ['_'.join([i.resname(), str(i.resid()), i.name()]) for i in x]
        retlist.append('-'.join(names))
    retlist.append('total_inframe')
    return '\t'.join(retlist) + '\n'


# define the model, subset to selection at this point
if args.coords:
    model = loos.selectAtoms(loos.loadStructureWithCoords(args.model, args.coords), args.selection)
else:
    model = loos.selectAtoms(loos.createSystem(args.model), args.selection) 

# add connectivity, anything within 1.65 A (loos default) is assumed covalently bound
if model.hasCoords():
    if not model.hasBonds():
        model.findBonds()
else: # if model doesn't have any coordinates, this analysis is nonsensical; exit with error.
    print 'Error: You must supply a model structure with coordinates.'
    exit(-1)

# read in the traj
traj = loos.pyloos.VirtualTrajectory(skip=args.skip, stride=args.stride)
for trajname in args.traj:
    traj.append(loos.pyloos.Trajectory(trajname, model, skip = args.skip, stride = args.stride))




# make hb detector object
model_hbs = loos.HBondDetector(args.bondlength, args.angle, model)

# select all nitrogens and oxygens in the model as putative donors/acceptors
putative = loos.selectAtoms(model, 'name =~ "^N" || name =~ "^O"')

# select all hydrogens to test which are DA pairs
hydrogens = loos.selectAtoms(model, 'name =~ "^H"')

# initialize donor and acceptor groups for appending
donors = []
acceptors = loos.AtomicGroup()

# build list of donor and acceptor groups from putative groups.
for i in putative:
    bound_to_hydrogen = False
    for h in hydrogens:
        if i.isBoundTo(h):
            donors.append((h,i)) # Associate donor with its hydrogen
            bound_to_hydrogen = True
    if not bound_to_hydrogen:
        acceptors.append(i)

# list of donor-H-acceptor triples present in reference crds
ref_dhas = []

# fill out the donor-H-acceptor list
for i in acceptors:
    for j in donors:
        if model_hbs.hBonded(i, j[0], j[1]):
            ref_dhas.append(np.array([i, j[0], j[1]])) 
            donors.remove(j) # assumes D-H pairs can only be used once
            break # don't check any extra D-Pairs

# can't append to np arrays so array this after constructiion.
ref_dhas = np.array(ref_dhas)

# build an atomic group of just these D-H---A triples
dha_group = loos.AtomicGroup()
for dha in ref_dhas:
    for atom in dha:
        dha_group.append(atom)

# make counters to track number of frames and hb tallies
frame_count = 0
total_hbs = 0
# iterate over all the frames in the traj, checking each frame for each HB.
if args.outfile:
    with open(args.outfile, 'w') as outfile:
        outfile.write(format_header(ref_dhas))
        for frame in traj:
            frame_count += 1
            total_inframe = 0
            outform = [frame_count] # list that will aggregate a row of data in outfile
            frame_hbs = loos.HBondDetector(args.bondlength, args.angle, frame)
            for dha in ref_dhas:
                if frame_hbs.hBonded(dha[0], dha[1], dha[2]): # if there's an HB, add to counters
                    total_hbs += 1 # add to total hbs
                    total_inframe += 1 # add to inframe total
                    outform.append(1) # add a 1 to the row and column you're on in outfile
                else:
                    outform.append(0) # didn't find an HB so add a zero to row in outfile
            outform.append(total_inframe) # add last column to outfile
            outfile.write('\t'.join(map(str,outform)) + '\n') # write tab delim row to outfile
else:
    for frame in traj:
        frame_count += 1
        frame_hbs = loos.HBondDetector(args.bondlength, args.angle, frame)
        for dha in ref_dhas:
            if frame_hbs.hBonded(dha[0], dha[1], dha[2]):
                total_hbs += 1 # if DHA observed in frame, add one to counter.

average_hbs = total_hbs/float(frame_count)
print frame_count, total_hbs, average_hbs


# write built atomic group of native D-H---A triples to PDB for Vis.
# if arg to give this flag is called for
if args.nativeHBs:
    with open(args.nativeHBs, 'w') as outfile:
        outfile.write(str(loos.PDB_fromAtomicGroup(dha_group)))