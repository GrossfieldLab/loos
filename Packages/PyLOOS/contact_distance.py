#!/usr/bin/env python3
"""
From a trajectory, create an all-to-all distance matrix analogous to rmsds, but
working in the space of residue-residue contacts
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
import loos.pyloos.options as options
import numpy
from scipy.spatial.distance import pdist

fullhelp = """
  Insert something useful here

  If more than 1 trajectory is specified, any --skip or --stride values will be
  applied to each.
"""

lo = options.LoosOptions(fullhelp)
lo.modelSelectionOptions()
lo.trajOptions()
lo.parser.add_argument('--sigma',
                       default=6.0,
                       type=float,
                       help="exponent for the logistic contact")
lo.parser.add_argument('--radius',
                       default=5.0,
                       type=float,
                       help="midpoint distance for the logistic contact")
# fix this so you don't need to give a value
lo.parser.add_argument('--skip_backbone',
                       default=False,
                       help="Consider only sidechains")
lo.parser.add_argument('--include_h',
                       default=False,
                       help="Include hydrogens")
lo.parser.add_argument('--outfile',
                       help="Name of outputted file ")

args = lo.parse_args()

system = loos.createSystem(args.model)

chain = loos.selectAtoms(system, args.selection)
if args.skip_backbone:
    chain = loos.selectAtoms(chain, '!backbone')
if not args.include_h:
    chain = loos.selectAtoms(chain, '!hydrogen')

residues = chain.splitByResidue()
num_residues = len(residues)

# if there's 1 trajectory, args.traj will be a string, otherwise
# it's a list
if type(args.traj) == str:
    t = loos.pyloos.Trajectory(args.traj, system,
                               skip=args.skip, stride=args.stride)
    traj = loos.pyloos.VirtualTrajectory(t)
else:
    t = loos.pyloos.Trajectory(args.traj[0], system,
                               skip=args.skip, stride=args.stride)
    traj = loos.pyloos.VirtualTrajectory(t)
    for t in args.traj[1:]:
        t = loos.pyloos.Trajectory(args.traj[0], system,
                                   skip=args.skip, stride=args.stride)
        traj.append(t)

# set up storage
contacts = numpy.zeros((len(residues) * (len(residues)-1) / 2, len(traj)),
                       type=float)

default_box = loos.GCoord(10000., 10000., 10000.)
frame_number = 0
for frame in traj:
    # checking each frame is wasteful, but shouldn't have a measurable
    # effect on performance
    if frame.isPeriodic():
        box = frame.periodicBox()
    else:
        box = default_box

    # compute residue-residue contacts and store
    for i in range(num_residues - 1):
        r1 = residues[i]
        for j in range(i+1, num_residues):
            contact = r1.logisticContact(r2, args.radius, args.sigma, box)
            index = (i*num_residues) + j
            contacts[frame_number, index] = contact

print("Finished computing contacts, starting distance calculation")

dists = pdist(contacts, 'euclidean')
square_dists = np.zeros([num_residues, num_residues], np.float)
for i in range(num_residues - 1):
    for j in range(i+1, num_residues):
        index = i*num_residues + j
        square_dists[i, j] = dists[index]
        square_dists[j, i] = dists[index]

if args.outfile:
    outfile = open(args.outfile, "w")
else:
    outfile = sys.stdout

square_dists.savetxt(outfile, header=lo.header())
