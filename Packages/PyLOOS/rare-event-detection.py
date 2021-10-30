#!/usr/bin/env python3
"""
Use Non-Negative Matrix Factorization to identify conformational transitions. 
Intended for use with GPCRs. Based on PW, Ligand-Dependent Conformational Transitions 
in Molecular Dynamics Trajectories of GPCRs Revealed by a New Machine Learning 
Rare Event Detection Protocol, journal, year, doi.

Grace Julien 2021
https://doi.org/10.3390/molecules26103059
"""

"""

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2021 Grace Julien, Grossfield Lab
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

import loos
import loos.pyloos
import loos.pyloos.options as options
import sklearn.decomposition
import sys
import numpy
#import matplotlib.pyplot as plt
from os.path import basename, splitext


fullhelp = """
  red.py: compute the makeup of components of residue-residue
  contacts and the weights of said components over the course 
  of a trajectory using Non-Negative Matrix Factorization
  as implemented by sk-learn.
  https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html
  Based on methodology outlined by Plante and Weinstein in 
  https://doi.org/10.3390/molecules26103059
  Window_length should be approximately the number of frames per nanosecond
  It is recommended to play with n_components, starting with "the expected number 
  of rare events in the trajectory based on mechanistic hypotheses of the
  conformational changes required for the molecular process," (reword or cite to Plante ^)
  and experimenting with +/- components until convergence is attained.

  Mandatory arguments:
  system_file: file describing system contents, e.g. a psf or pdb
  selection: selection string for which residues to look at
  spatial_out_file: name for the component composition file
  temporal_out_file: name for the component weight file
  traj_file: 1 trajectory file
  window_length: window length for smoothing function

  Options
  --cutoff: distance for atom-atom contacts, defaults to 4.0 Ang
  --no-hydrogens: ignore hydrogens when looking for contacts.
  --no-backbone: ignore the protein or RNA backbone.  If this flag isn't
                 given, you'll have a contact probability of 1 between
                 consecutive residues in a chain
  --n_components: number of components for NMF, defaults to 5
  --max_iterations: maximum number of iterations for NMF, defaults to 200
  --fullhelp: produce this message

  This program does not explicitly handle periodicity; it assumes you've
  already fixed any periodicity issues before you ran it.

  """

def moving_average3d(a, n = 30):
    ret = numpy.cumsum(a, dtype=float, axis = 2)
    ret[0:,0:,n:] = ret[0:,0:,n:] - ret[0:,0:,:-n]
    return ret[0:,0:,n-1:] / n

def moving_ave1d(a, n = 30):
    ret = numpy.cumsum(a, dtype=float, axis = 1) #1 for across, 0 for down
    ret[0:,n:] = ret[0:,n:] - ret[0:,:-n] #going across is going across the frames (1), going down is going through the contacts for one frame (0)
    return ret[0:,n-1:] / n
    

def LowerTriIndex(row, col, n):
    index = (row*(row-1)/2 + col)
    return int(index)

if __name__ == '__main__':
    lo = options.LoosOptions(fullhelp)
    lo.modelSelectionOptions()
    lo.trajOptions()


    lo.parser.add_argument('--spatial_out_file',
                            #required=True,
                            help="File with the component composition from NMF (W)")
    lo.parser.add_argument('--temporal_out_file',
                            #required=True,
                            help="File with the weights of each component from NMF (H)")
    lo.parser.add_argument('--cutoff', type=float,
                            help="Cutoff distance for contact", default=4.0)
    # TODO: add a number of contacts option
    lo.parser.add_argument('--no_hydrogens', action='store_true',
                            help="Don't include hydrogens")
    lo.parser.add_argument('--no_backbone', action='store_true',
                            help="Don't include the backbone")
    lo.parser.add_argument('--window_length', type=int,
                            required=True,
                            help="Window length for smoothing function, number of frames per 30 nanoseconds")
    lo.parser.add_argument('--n_components', type=int,
                            help="Number of components for NMF function", default = 5)
    lo.parser.add_argument('--max_iterations', type=int,
                            help="Maximum iterations for NMF function", default = 200)
    args = lo.parser.parse_args()


    header = lo.header()


    system = loos.createSystem(args.model)


    traj = loos.pyloos.Trajectory(args.traj[0], system)

    if args.no_hydrogens:
        no_hydrogens = loos.selectAtoms(system, "!hydrogen")
        target = loos.selectAtoms(no_hydrogens, args.selection)
    else:
        target = loos.selectAtoms(system, args.selection)

    residues = target.splitByResidue()
    # now remove the backbone -- doing before the split loses the glycines
    if args.no_backbone:
        residues = list([loos.selectAtoms(r, "!backbone") for r in residues])

    frac_contacts = numpy.zeros([len(residues), len(residues), len(traj)],
                                numpy.float64)
    fc = numpy.zeros([int((len(residues)-1)*len(residues)/2), len(traj)], numpy.float64)


    for (frame, frame_id) in zip(traj, range(len(traj))): #frac_contact (len res len res len traj)
        #frac = numpy.zeros([subsum(len(residues))], numpy.float64)
        for i in range(len(residues)):
            for j in range(i+1, len(residues)):
                if residues[i].contactWith(args.cutoff, residues[j]):
                    if j>i:
                        fc[LowerTriIndex(j, i, len(residues)), frame_id] += 1.0

    #print(frac_contacts)
    #print(frac_contacts.shape) #upper triangular matrix
    ave = moving_ave1d(fc, args.window_length)
    #ave = moving_average3d(frac_contacts, args.window_length)
    #flat = ave.reshape((ave.shape[0]*ave.shape[1]), ave.shape[2])
    #print(flat.shape)

    model = sklearn.decomposition.NMF(n_components = args.n_components, 
                                      init = 'nndsvd', 
                                      max_iter = args.max_iterations)
    W = model.fit_transform(ave)
    H = model.components_

    print(W.shape)
    #print(H[1].shape)

    #numpy.savetxt(args.spatial_out_file, W, header=header)
    #numpy.savetxt(args.temporal_out_file, H, header=header)


    #numpy.savetxt(args.out_file, average, header=header)
