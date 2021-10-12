#!/usr/bin/env python3

"""
Unwraps a systems perdoic box of for each frame of a
given selection. This creates an output for calcuting the
means sqaured diffusion.
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
def fullhelp():
    print("""lipid_unwrap.py takes the coordinates from a trajectory file to unwrap
    the periodic box for a selection. At the moment the only rectangular
    box based periodic box conditions are able to be used. The subset of atoms
    that you use could be any valid loos selection, but the program is intended
    to work on lipids,tracking the geometric center of a molecule at the
    current frame and comparing that position to the previous one.
    If the difference between the two positions vary by more than half of the
    periodic box length (in x,y, or z) a correction factor is added to the
    current position to correct the positions it moves outside the boundaries.

    If desired, a pdb and DCD of the unwrapped trajectory can be. This is done by
    taking every molecule if the selection to (0,0,0) then translating it to the
    unwrapped. For making the trajectory, the assumption is made that the
    Lx, Ly or Lz of the periodic box will not be greater than 1000000 Angstroms.

    The output text file containing the radial position of the atom comparing
    it to the center of mass to the bilayer.

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
        --output_traj: output a dcd of the unwrapped trajectory default(--output_traj 0)
        --output_prefix: prefix to use for the output default is output
        --no_center: Variable to turn off the recentering of the trajectory default(--no_center 1)

    Example:
    python3 lipid_unwrap.py --ouput_traj 1 --output_prefix='unwrapped_foo' foo.psf 'resname == "POPC"' foo.dcd

    We will obtain a dcd of the unwrapped trajectory with the name of
    unwrapped_foo.dcd a pdb of unwrapped_foo.pdb to be used to view the dcd, and
    unwrapped_foo.txt that has the radial distantness of the center of the lipid
    to the COM of bilayer at that frame for the POPC lipids in the membrane
    """)


def findFix():
    """
    find Fix is the takes the center before it is unwrapped then checks the lipids
    until to see if the new position has moved more than half the box then the box
    length is added to the position of that atom to negate the wrapping. If the
    corrdinate of fixed atom still is lager than half the box a factor is multiple to
    the box length added to the system until the distance is within reason
    """
    fix = loos.GCoord()
    bad = {0:1,1:1,2:1,3:1,4:1,5:1}  # scaling factors of each
    check = True  # checks if the condtion is true
    while check:  # while the ditstance propational to a jump
        diff = centers[index] - prev_centers[index]  # final minus statrs
        tracker = 0  # tracks which if it enterhs
        # looking for points that jumped
        bad_coord = {}
        if diff.x() > half_box.x():
            fix[0] -= box.x() * bad[tracker]
            bad_coord[tracker] = bad[tracker]
        elif diff.x() < -half_box.x():
            fix[0] += box.x() * bad[tracker]
            tracker = 1
            bad_coord[tracker] = bad[tracker]
        if diff.y() > half_box.y():
            fix[1] -= box.y() * bad[tracker]
            tracker = 2
            bad_coord[tracker] = bad[tracker]
        elif diff.y() < -half_box.y():
            fix[1] += box.y() * bad[tracker]
            tracker = 3
            bad_coord[tracker] = bad[tracker]
        if diff.z() > half_box.z():
            fix[2] -= box.z() * bad[tracker]
            tracker = 4
            bad_coord[tracker] = bad[tracker]
        elif diff.z() < -half_box.z():
            fix[2] += box.z() * bad[tracker]
            tracker = 5
            bad_coord[tracker] = bad[tracker]

        # the new position is the wrapped position plus the fix
        ans = centers[index] + fix
        if prev_centers[index].distance(ans) > half_box.x():
            # can assume here that the bad cord is empty
            if len(bad_coord) != 0:
                for i in bad_coord.keys():
                    bad[i] += 1
        else:
            check = False
    return fix


class FullHelp(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        kwargs['nargs'] = 0
        super(FullHelp, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        fullhelp()
        parser.print_help()
        setattr(namespace, self.dest, True)
        parser.exit()



if __name__ in '__main__':
    parser = argparse.ArgumentParser(description="Unwrap PBC lipids")
    parser.add_argument('system_file',
                        help="File describing the system")
    parser.add_argument('selection_string',
                        help="Selection string describing which residues to use")
    parser.add_argument('traj_file',
                        help="File contraing the trajecotry")
    parser.add_argument('--fullhelp',
                        help="Print detailed description of all options",
                        action=FullHelp)
    parser.add_argument('--skip',
                        type=int,
                        default=0,
                        help='# of frame to skip')
    parser.add_argument('--stride',
                        type=int,
                        default=1,
                        help="Read every nth frame")
    parser.add_argument('--output_traj',
                        type=int,
                        default=0,
                        help='produce an unwrapped trajectory')
    parser.add_argument('--output_prefix',
                        default='unwrapped_output',
                        type=str,
                        help='name of the tractory file to write (DCD format)')
    parser.add_argument('--no_center',
                        default=1,
                        type=int,
                        help='Do not wrap the trajecotry default is to center the lipid selection')
    args = parser.parse_args()
    pre = args.output_prefix
    header = " ".join([f"{i}" for i in sys.argv])

    system = loos.createSystem(args.system_file)
    traj = loos.pyloos.Trajectory(args.traj_file,
                                  system,
                                  stride=args.stride,
                                  skip=args.skip)

    # select the atoms for the lipids
    # split into the individual molecules of lipid
    lipid_sel =  loos.selectAtoms(system, args.selection_string)
    lipids =lipid_sel.splitByMolecule()
    # boolean to check if this is frame is frame 0
    first = True
    # making a numpy array to save the coorinates
    num_lipids = len(lipids)
    num_frames= len(traj)
    all_centers = np.zeros((num_lipids,num_frames))
    # loop into the frame
    if args.output_traj:
        outtraj = loos.DCDWriter(pre + ".dcd")

    # this is the order of the lipid resname-resid-segname
    header += '\nresname-resid-segid'
    big_list = " ".join(
        [f"{lipids[i][0].resname()}-{lipids[i][0].resid()}-{lipids[i][0].segid()} " for i in range(len(lipids))])
    header += f"\n{big_list}"
    header += "\nLipid_1_center_R Lipid_2_center_R Lipid_3_center_R ..."
    for frame in traj:
        # shift the lipids such that their total centroid is at the origin
        COM_R = loos.selectAtoms(system, args.selection_string).centroid()
        system.translate(-COM_R)
        centers = []
        for index in range(len(lipids)):
            # geometric center of the lipid
            centers.append(lipids[index].centroid())
            # position 0 for frame 0 is always the starting position
        # in frame one nothing should have moved
        R_coord = [traj.index()]
        if not first:
            box = frame.periodicBox()  # get a G cord of the PBC
            half_box = 0.5 * box  # a square box's center is at the half box way between all the positions
            for index in range(len(lipids)):
                fix = findFix()
                centers[index] += fix
                R_coord.append(centers[index].length())

        if args.output_traj:
            # loop through the list of center locations and
            moved_lipids = loos.AtomicGroup() # blank Atoms group
            for lipid, center in zip(lipids, centers):
                lipid.centerAtOrigin()  # move that lipid to the center
                # translate to the center
                lipid.translate(center) # then move the lipids to the center of
                # the unwrapped point
                moved_lipids.append(lipid) # add the lipid to the blank Atom loos

            # setting the pbc to be large
            moved_lipids.periodicBox(loos.GCoord(1000000, 1000000, 1000000)) # set the
            if args.no_center:
                moved_lipids.centerAtOrigin()

            outtraj.writeFrame(moved_lipids)

            # write a pdb
            if first:
                moved_lipids.periodicBox(loos.GCoord(1000000, 1000000, 1000000))
                # change the box size
                frame_copy = lipid_sel.copy()
                frame_copy.renumber()
                pdb = loos.PDB.fromAtomicGroup(frame_copy)
                pdb.remarks().add(header)
                with open(f"{pre}.pdb", "w") as out:
                    out.write(str(pdb))
            # now loop into the array of R_coord input the values
            for index in range(len(lipids)):
                # find the lenght of the centers indexs
                all_centers[index,traj.index()] = centers[index].length()
        prev_centers = centers
        first = False
    np.savetxt(f'{pre}.txt', all_centers, header=header)
