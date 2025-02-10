#!/usr/bin/env python3
"""
Build a lipid library that can be used for an OMG run

Alan Grossfield,  University of Rochester Medical Center, 2019
"""

"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2019 Tod Romo, Grossfield Lab
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
import os
import loos
import loos.pyloos


def main():
    system_file = sys.argv[1]  # must have connectivity
    selection = sys.argv[2]
    stride = int(sys.argv[3])
    center_atom = sys.argv[4]  # usually the phosphorus
    library_dir = sys.argv[5]
    traj_files = sys.argv[6:]


    system = loos.createSystem(system_file)

    all_lipids = loos.selectAtoms(system, selection)

    # get the individial lipid molecules, and the centering atoms
    lipids = all_lipids.splitByMolecule()
    centers = loos.selectAtoms(all_lipids, 'name == "' + center_atom + '"')

    # make the directory to hold the library
    try:
        os.mkdir(library_dir)
    except FileExistsError:
        print("Library directory %s exists" % (library_dir))
        pass

    # set up the virtual trajectory
    first = loos.pyloos.Trajectory(traj_files[0], system, stride=stride)
    vtraj = loos.pyloos.VirtualTrajectory(first)

    for t in traj_files[1:]:
        traj = loos.pyloos.Trajectory(t, system, stride=stride)
        vtraj.append(traj)

    x_axis = loos.GCoord(1.0, 0., 0.)
    index = 0
    for frame in vtraj:
        for i in range(len(lipids)):

            # rotate the lower leafet "right side up"
            c = lipids[i].centroid()
            if c.z() < 0:
                lipids[i].rotate(x_axis, 180.)

            # put the centering atom at the origin
            center = centers[i].coords()
            lipids[i].translate(-center)

            # covert to a PDB and write it out
            pdb = loos.PDB.fromAtomicGroup(lipids[i])

            # clean up the numbering
            pdb.renumber()
            for a in pdb:
                a.resid(1)

            filename = library_dir + "/lipid_" + str(index) + ".pdb"
            with open(filename, "w") as outfile:
                outfile.write(str(pdb))

            index += 1

if __name__ == "__main__":
    main()
    