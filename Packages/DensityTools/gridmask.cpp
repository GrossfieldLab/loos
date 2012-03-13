/*
  gridmask.cpp


  Given an int grid that represents picked blobs, use this as a mask
  against a double grid...
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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
*/


#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>
#include <algorithm>
#include <sstream>
#include <limits>

#include <DensityGrid.hpp>

using namespace std;
using namespace loos;
using namespace DensityTools;


/*
 * Applies a mask to an edm grid.  This can select a particular
 * region within the grid.  
 * 
 * Example:
 *    gridmask < foo_grid foo_picked > masked_foo
 * This will apply the mask called \"foo_picked\" to the grid
 * \"foo_grid\" and write out a new grid, \"masked_foo\". This
 * assumes foo_picked was created with the \"pick_blob\" tool.
 *
 * This output my be converted to an xplor map and then used for 
 * visualization within pymol.
 * 
*/



int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr <<
      "SYNOPSIS\n\tExtracts a region of density given a mask grid\n"
      "\nDESCRIPTION\n\tThis tool will zero out any unwanted density\n"
      "given a density grid and a grid mask.  The grid mask is an integer\n"
      "grid.  Any non-zero element of the grid mask means that the corresponding\n"
      "density from the density grid will be copied to the output grid.\n"
      "All other locations will have a zero density value.  (Think of an alpha-mask\n"
      "in gimp or photoshop).\n"
      "\nEXAMPLES\n\tgridmask <density.grid mask.grid >masked_density.grid\n"
      "This will apply the mask.grid mask to the density.grid, writing the output\n"
      "to masked_density.grid\n"
      "\n"
      "\tblobid --threshold 1 <foo.grid >foo_id.grid\n"
      "\tpick_blob --model foo.pdb --selection 'resid == 65' <foo_id.grid >foo_picked.grid\n"
      "\tgridmask <foo.grid foo_picked.grid >foo_masked.grid\n"
      "This example will first threshold the density at 1.0, then it will find the blob\n"
      "closest to residue 65.  This blob is then used as a mask for the original density\n"
      "grid.  foo_picked.grid therefore contains the actual density values, but with\n"
      "all extraneous density removed.\n";
      
    cerr << "Usage- gridmask <edm_grid mask_grid >masked_edm_grid\n";
    exit(-1);
  }

  DensityGrid<int> mask;
  ifstream ifs(argv[1]);
  if (!ifs) {
    cerr << "Error - cannot open " << argv[1] << " for reading.\n";
    exit(-10);
  }

  ifs >> mask;

  DensityGrid<double> data;
  cin >> data;

  DensityGridpoint dims = data.gridDims();
  DensityGridpoint ddims = mask.gridDims();
  if (dims != ddims) {
    cerr << "Error - differing dimensions between mask and density grids.\n";
    exit(-10);
  }

  long k = dims[0] * dims[1] * dims[2];
  for (long i=0; i<k; i++)
    if (!mask(i))
      data(i) = 0.0;

  cout << data;
}
