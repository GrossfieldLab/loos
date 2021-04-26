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
from scipy.spatial.distance import pdist, squareform

fullhelp = """

Compute the frame-frame distance matrix for a trajetory in contact space.

This program is analogous to rmsds, except the distance is computed using a
logistic contact function for each pair of residues. For each pair of residues,
the degree of contact is computed using a logistic switching function
C(r) = 1/(1 + (r/radius)^sigma

So, each structure with N residues is represented by a vector with N*(N-1)/2
entries, each storing CC = 1-C(r). The distance between 2 structures is
computed as sqrt(sum CC). We then write out the full square matrix, suitable
for consumption by a clustering program such as our cluster-kgs.

Note: the final matrix is MxM, where M is the number of frames in the trajectory
or trajetories, so it will get very big very fast (eg ~2.5 GB for 1000 frames).
So, if you've written a lot of frames, you'll probably want to downsample the
calculation using the --stride option.

You can control the nature of the contact calculation using these flags:
    --sigma: the exponent in the contact equation. Must be an integer,
             defaults to 6. Larger values make the 1->0 transition sharper.
    --radius: mid-point of the switching distance. The default is 6 ang,
            which seems reasonable for sidechain-sidechain distances.

The program uses the selection string supplied on the command line to decide
which residues to look at. For your convenience, we've added 2 additional
flags:
    --skip_backbone: adds "!backbone" to the selection to remove the backbone of
                    protein, RNA, or DNA
    --include_h: by default, the calculation uses only the heavy atoms, since
                hydrogens don't generally shift the location of the centroid
                appreciably. If you want the hydrogens, use this flag.

By default, the matrix of distances is written to stdout, so you could in
principle stream it into something like cluster-kgs without saving it.
If you want to save it to a file, use the --outfile flag.

TODO: add an option to write out the contact structures themselves, so someone
could input them into PCA, tICA, or something like that.

"""

lo = options.LoosOptions("Compute distance matrix using contacts", fullhelp)
lo.modelSelectionOptions()
lo.trajOptions()
lo.parser.add_argument('--sigma',
                       default=6,
                       type=int,
                       help="exponent for the logistic contact")
lo.parser.add_argument('--radius',
                       default=5.0,
                       type=float,
                       help="midpoint distance for the logistic contact")
lo.parser.add_argument('--skip_backbone',
                       default=False,
                       action='store_true',
                       help="Consider only sidechains")
lo.parser.add_argument('--include_h',
                       default=False,
                       action='store_true',
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
        t = loos.pyloos.Trajectory(t, system,
                                   skip=args.skip, stride=args.stride)
        traj.append(t)

num_pairs = len(residues) * (len(residues)-1) // 2
contacts = numpy.zeros((len(traj), num_pairs), numpy.float)

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
    index = 0
    for i in range(num_residues - 1):
        r1 = residues[i]
        for j in range(i+1, num_residues):
            contact = r1.logisticContact(residues[j],
                                         args.radius,
                                         args.sigma,
                                         box)
            # High contact -> low distance
            # However, since we're taking distances between structures, you
            # get pretty much the same answer either way.
            contacts[frame_number, index] = 1. - contact
            index += 1
    frame_number += 1


dists = pdist(contacts, 'euclidean')
square_dists = squareform(dists)

if args.outfile:
    outfile = open(args.outfile, "w")
else:
    outfile = sys.stdout

#trajectories = "Trajectories: " + " ".join(args.traj)
#frame_boundaries = "First frames: " + " ".join(str(x) for x in traj.frameBoundaries())
frame_boundaries = traj.frameBoundaries()

traj_header = [" traj  start   end     filename"]
if (type(args.traj) == str):
    s = ["0", "0", str(len(traj)-1), args.traj]
    traj_header.append("\t".join(s))
else:
    for i in range(len(args.traj)-1):
        s = [str(i), str(frame_boundaries[i]), str(frame_boundaries[i+1]-1),
             args.traj[i]]
        traj_header.append("\t".join(s))
    i = len(args.traj) - 1
    s = [str(i), str(frame_boundaries[i]), str(len(traj)-1),
         args.traj[i]]
    traj_header.append("\t".join(s))
traj_header = "\n".join(traj_header)

header = "\n".join(str(x) for x in [lo.header(), traj_header])

numpy.savetxt(outfile, square_dists, header=header)
