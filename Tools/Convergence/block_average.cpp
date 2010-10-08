/* 
  Compute block-averaged standard error for a time series
  using the Flyvbjerg & Petersen approach.  Outputs the block averaged
  standard error as a function of block size -- you'll need to plot it
  and estimate the plateau value.
 
  References: 
  
  Flyvbjerg, H. & Petersen, H. G. Error estimates on averages 
  of correlated data J. Chem. Phys., 1989, 91, 461-466

  Grossfield, A., and Zuckerman, D. M. Quantifying uncertainty and sampling 
  quality in biomolecular simulations, Ann. Reports in Comp. Chem., 2009, 
  5, 23-48

 
  Alan Grossfield
  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School
 
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Alan Grossfield
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

using namespace std;
using namespace loos;


void Usage()
    {
    cerr << "Usage: block_average TimeSeriesFile column max_blocks skip"
         << endl;
    cerr << endl;
    cerr << "TimeSeriesFile is a columnated text file.  Blank lines and "
         << "lines starting with \"#\" are ignored"
         << endl;
    }

int main(int argc, char*argv[])
{
if ( (argc <= 1) ||
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 5)
   )
    {
    Usage();
    exit(-1);
    }

// Print the command line arguments
cout << "# " << invocationHeader(argc, argv) << endl;

char *datafile = argv[1];
int column = atoi(argv[2]);
unsigned int max_blocks = atoi(argv[3]);
unsigned int skip = atoi(argv[4]);

// Read the TimeSeries file
TimeSeries<float> data = TimeSeries<float>(datafile, column);
unsigned int num_points = data.size();

// validate the arguments
if (skip < 0)
    {
    cerr << "Command line value for skip must be >= 0" << endl;
    exit(-1);
    }

if (skip > num_points)
    {
    cerr << "You set skip ( " << skip << " ) greater than the number "
         << "of points in the trajectory ( " << num_points << " )."
         << endl
         << "This doesn't work." 
         << endl;
    exit(-1);
    }

// Remove the equilibration time
data.set_skip(skip);
num_points -= skip;

if (max_blocks < 0)
    {
    cerr << "Command line value for max_blocks must be >= 0" << endl;
    exit(-1);
    }

if (max_blocks > num_points)
    {
    cerr << "You set max_blocks ( " << max_blocks << " ) greater than "
         << "the number of points in the trajectory minus the number skipped "
         << "( " << num_points << " )."
         << endl
         << "This doesn't work." 
         << endl;
    exit(-1);
    }

// Loop over the number of blocks, computing the variance of the averages for
// each number of blocks

cout << "# Num_Blocks\tBlockLen\tStdErr" << endl;

for (int i=max_blocks; i>=2; i--)
    {
    float time = (float)num_points / i;
    float variance = data.block_var(i);
    float std_err = sqrt(variance/i);
    cout << i << "\t\t"
         << time << "\t\t"
         << std_err
         << endl;
    }


}
