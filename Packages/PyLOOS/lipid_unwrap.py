# @Author: Kyle Billings <kbillings>
# @Date:   2021-09-22T18:00:20-04:00
# @Email:  krb0073@mix.wvu.edu
# @Filename: lipid_un.py
# @Last modified by:   kbillings
# @Last modified time: 2021-09-26T15:07:15-04:00


#!/usr/bin/env python3

"""
Track a base stacking through a trajectory.  Intended for use with
nucleic acids, to identify base stacking.
"""

"""
  This file is part of LOOS.
  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2021 Alan Grossfield, Grossfield Lab
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
import numpy as np
import argparse
# making the fullhelp for the program
def fullhelp():
    print("""lipid_unwrap.py takes the coordinates from a trajectory file to unwrap
    the periodic box(PB) for a selection .At the moment the only rectangular
    box based periodic box conditions are able to be used. The subset of atoms
    that you use could be any valid loos selection ,but the program is intended
    to work on lipids,tracking the geometric center of a molecule at the
    current frame and comparing that position to the previous one.
    If the difference between the two positions vary by more than half of the
    P.B  length (in x,y, or z) a correction factor is added to the current position
    to correct the positions it moves outside the boundaries.

     If desired, a pdb and DCD of the unwrapped trjectory can be. This is done by
     taking every molecule if the selection to (0,0,0) then translating it to the
     unwrapped. For making the trajecotry,the assumtion is made that the
     Lx,Ly or Lz of the P.B will not be greater than 1,000 Angstroums.

     The output txt file contrain the radial position of the atom comparing it to the
     center of mass to the bilayer.

    Usage: lipid_unwrap.py [OPTIONS] selection_string system_file trajectory_file
        --  selection_string identifies the atoms to be analyzed. They will be
            split by molecule number calculate position for every lipid adjusted
            without the wrapping
        --  system_file is a PDB, PSF, prmtop, etc
        --  trajectory_file eg DCD, xtc. Its contents must precisely
            match the system file

    Options:
        --skip: # of residues to exclude from the front of each trajectory
        --stride: how to step through the trajectory (eg --stride 10 will read
              every 10th frame)
        --fullhelp: prints this message
        --output_traj: output a dcd of the unwrapped trajectory defult(--output_traj 0)
        --output_prefix: prefix to use for the output defult is output

    Exmaple:
    python /media/bak12/Analysis/unwrap/lipid_unwrap.py --ouput_traj 1 --output_prefix='unwrapped_foo'
    foo.psf 'resname == "POPC"' foo.dcd

    We will obatian a dcd of the unwrapped trajectory with the name of unwrapped_foo.dcd
    a pdb of unwrapped_foo.pdb to be used to veiw the dcd, and unwrapped_foo.txt that
    has the radial distacne of the ceneter of the lipid to the COM of bilayer
    at that frame for the POPC lipids in the membrane""")
def findFix():
    """ find Fix is the takes the center befor it is unwrapped then checks the lipids
    until to see if the new postion has moved more than half the box then the box
    lenght is added to the postion of that atom to negate the wrapping. If the
     corrdinate of fixed atom still is lager than half the box a factor is multiped to
      the box lenght added to the system untill the distacne is within reason"""
    fix = loos.GCoord() # blank GCorrd() == (0,0,0)
    a ,b ,c, d,e,f = 1, 1,1,1,1,1 # scaling factors of each
    check = True # checks if the condtion is true
    while check: # while the ditstance propational to a jump
        diff = centers[index] -  prev_centers[index]  # final minus statrs
        tracker = 0 # tracks which if it enterhs
        # looking for points that jumpedd
        if diff.x() >  half_box.x():
            fix[0] -= box.x() * a
        elif diff.x() < -half_box.x():
            fix[0] += box.x() * b
            tracker =1
        if diff.y() > half_box.y():
            fix[1] -= box.y() * c
            tracker = 2
        elif diff.y() < -half_box.y():
            fix[1] += box.y() * d
            tracker = 3
        if diff.z() > half_box.z():
            fix[2] -=  box.z() * e
            tracker = 4
        elif diff.z() < -half_box.z():
            fix[2] += box.z() * f
            tracker = 5
        # the new postion is the wrapped postion plus the fix
        ans = centers[index] + fix
        # if the distacne from the old location to the new is large than the
        ## half_box
        if prev_centers[index].distance(ans) > half_box.x():
            # using the value of tracker we find which if was used to fix the
            ## bad corrdinate
            if tracker == 0:
                # if the ans - previous is less than half box the jump
                ## was from the postive side of the pbc to the negative
                if ans.x() - prev_centers[index].x() < half_box.x():
                    a += 1
            elif tracker == 1:
                # is the diffrecne is grater than half_box the jump was from
                ## negative to postive
                if ans.x() -prev_centers[index].x() > half_box.x():
                    b += 1
            elif tracker == 2:
                if ans.y() - prev_centers[index].y() < half_box.y():
                    c += 1
            elif tracker == 3:
                if ans.y() - prev_center[index].y() > half_box.y():
                    d += 1
            elif tracker == 4:
                if ans.z() - prev_centers[index].z() < half_box.z():
                    c += 1
            elif tracker == 5:
                if ans.z() - prev_center[index].y() > half_box.z():
                    d += 1
        else:
            # if there was no large jump we exit the loop
            check = False
    return fix
class FullHelp(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        kwargs['nargs'] = 0
        super(FullHelp, self).__init__(option_strings, dest, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        fullhelp()
        parser.print_help()
        setattr(namespace,self.dest,True)
        parser.exit()

parser = argparse.ArgumentParser(description="Unwrap PBC lipids")
parser.add_argument('system_file', help="File describing the system")
parser.add_argument('selection_string',help="Selection string describing which residues to use")
parser.add_argument('traj_file',help="File contraing the trajecotry")
parser.add_argument('--fullhelp',help="Print detailed description of all options",action=FullHelp)
parser.add_argument('--skip',type=int,default=0,help='# of frame to skip')
parser.add_argument('--stride' , type=int,default=1,help="Read every nth frame")
parser.add_argument('--output_traj',type=int,default=0,help='produce an unwrapped trajectory')
parser.add_argument('--output_prefix',default='unwrapped_output',type=str,help='name of the tractory file to write (DCD format)')
args = parser.parse_args()
pre = args.output_prefix
header = " ".join([f"{i}" for i in sys.argv])
system = loos.createSystem(args.system_file)
# load the trajectory
traj = loos.pyloos.Trajectory(args.traj_file,system,stride=args.stride , skip=args.skip)
# select the atoms for the lipids
## spilt into the invidaul molecules of lipid
lipids = loos.selectAtoms(system,args.selection_string).splitByMolecule()
# boolean to check if this is frame is frame 0
first = True
# loop into the frame
## this is the order of the lipid resname-resid-segname
### this will be usefull later
if args.output_traj: # fixes the issues where an emptry traj is made when
## not wanting one out 
    outtraj = loos.DCDWriter(pre +".dcd")
header += '\nresname-resid-segid'
big_list = " ".join([f"{lipids[i][0].resname()}-{lipids[i][0].resid()}-{lipids[i][0].segid()} " for i in range(len(lipids))])
header += f"\n{big_list}"
header += "\nframe Lipid_1_center_R Lipid_2_center_R Lipid_3_center_R ..."
for frame in traj:
    # loop into the lipids
    COM_R = loos.selectAtoms(system,args.selection_string).centerOfMass()
    centers = []
    for index in range(len(lipids)):
        centers.append(lipids[index].centroid()) # geometric center of the lipid
        # position 0 for frame 0 is always the starting position
    # in frame one nothing should have moved
    R_coord = [traj.index()]
    if not first:
        box = frame.periodicBox() # get a G cord of the PBC
        half_box = 0.5 * box  # a square box's cneter is at the halfb way between all the postions
        for index in range(len(lipids)):
            fix =findFix()
            centers[index]+= fix
            R_coord.append(np.sqrt((centers[index].x()-COM_R.x())**2 +(centers[index].y()-COM_R.y())**2+(centers[index].z()-COM_R.z())**2))
            # fix the werid jump issssues if the atoms move more than wanted
            # if the user has the flag on for the trajectory then this is True
    if args.output_traj:
        # loop through the list of center locations and
        moved_lipids = loos.AtomicGroup()
        for lipid , center in zip(lipids,centers):
            lipid.centerAtOrigin() # move tha lipid to the center
            # translate to the center
            lipid.translate(center)
            moved_lipids.append(lipid)
        # setting the pbc to be large
        moved_lipids.periodicBox(loos.GCoord(1000,1000,1000))
        if len(lipids) >= 5:
            moved_lipids.centerAtOrigin()
        outtraj.writeFrame(moved_lipids)
        if first == True:
            temp_lipids = loos.selectAtoms(system,args.selection_string)
            frame_copy = temp_lipids.copy()
            temp_lipids.periodicBox(loos.GCoord(1000,1000,1000))
            frame_copy.pruneBonds()
            frame_copy.renumber()
            pdb = loos.PDB.fromAtomicGroup(frame_copy)
            pdb.remarks().add(header)
            with open(f"{pre}.pdb","w") as out:
                out.write(str(pdb))
    all = [traj.index()]
    if first:
        for index in range(len(lipids)):
            R_coord.append(np.sqrt((centers[index].x()-COM_R.x())**2 +(centers[index].y()-COM_R.y())**2+(centers[index].z()-COM_R.z())**2))
        all_centers = np.array(R_coord)
    else:
        R_coord = np.array((R_coord))
        all_centers = np.row_stack((all_centers,R_coord))
    prev_centers = centers
    first = False
np.savetxt(f'{pre}.txt',all_centers,header=header)
