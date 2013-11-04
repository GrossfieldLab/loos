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

string fullHelpMessage(void)
    {
    string s = 
"\n"
"SYNOPSIS\n"
"\n"
"Apply block averaging to estimate standard error of timeseries data\n"
"\n"
"DESCRIPTION\n"
"\n"
"This tool performs a block averaging analysis in order to estimate the \n"
"standard error of a set of time series data.  It takes as input a text\n"
"file with white-space delimited data in columns (each time point is a \n"
"row), and returns the estimated standard error as a function of block \n"
"size.  \n"
"\n"
"The command line arguments are as follows:\n"
"\n"
"block_average TimeSeriesFile column max_blocks skip\n"
"\n"
"TimeSeriesFile      columnated text file (blank lines and lines starting \n"
"                    with \"#\" are ignored) containing the time series data\n"
"column              which column to use for analysis (1-based)\n"
"max_blocks          maximum number of blocks to use in the analysis\n"
"skip                number of frames to skip from the beginning of the \n"
"                    trajectory\n"
"  \n"
"The algorithm used is in essence that of Flyvbjerg and Petersen [Ref 1],\n"
"and is intended to estimate the standard error for a correlated time\n"
"series.  For uncorrelated data, the standard error can be estimated as\n"
"\n"
"SE = sqrt(var(a) / N) = stdev(a) / sqrt(N)\n"
"\n"
"where \"a\" is the quantity of interest and \"N\" is the number of points.\n"
"When the data has correlations, as is the case for nearly all molecular\n"
"dynamics or Monte Carlo simulations, this formula significantly \n"
"underestimates the statistical uncertainty.  \n"
"\n"
"The block averaging algorithm works by breaking the \"N\" data points\n"
"into \"M\" equal-sized contiguous blocks , computing the average within \n"
"each block, and then combining them to get the standard deviation in the \n"
"averages.  By tracking how that standard dev changes as a function of block \n"
"size, we can estimate the standard error in the limit of inifinite\n"
"block size, which is an estimate of the true standard error.  \n"
"\n"
"As the blocks get longer, there are fewer of them, and their variance\n"
"can get very noisy.  If you've got really good data, the at long block \n"
"time will be pronounced before the curve gets noisy.  If not, you can\n"
"estimate the plateau value by averaging the values at the last few block\n"
"sizes (in a plot of std err vs. block size).  If there is no plateau (e.g.\n"
"the curve is still systematically rising), your data is sufficiently \n"
"unconverged that the statistical error cannot be estimated.\n"
"\n"
"It is important to note that block averaging can significantly\n"
"underestimate the standard error for extremely undersampled systems, \n"
"because it is entirely based on what has been seen in the trajectory.  For\n"
"example, in the case of a 2-state system with different positions along\n"
"a reaction coordinate x, a very short simulation might only have population\n"
"in 1 state; block averaging this data would produce a small estimated \n"
"uncertainty, because the data looks very homogeneous.  Basically, the\n"
"analysis can't know what it hasn't seen.\n"
"\n"
"Note: If the number of blocks doesn't evenly divide the number of points,\n"
"then the remainder will be discarded from the end of the trajectory.\n"
"\n"
"\n"
"See references 1 and 2 for more discussion of the block averaging algorithm.\n"
"\n"
"1.  Flyvbjerg, H. & Petersen, H. G. Error estimates on averages \n"
"    of correlated data J. Chem. Phys., 1989, 91, 461-466\n"
"\n"
"2.  Grossfield, A., and Zuckerman, D. M. Quantifying uncertainty and \n"
"    sampling quality in biomolecular simulations, Ann. Reports in Comp. \n"
"    Chem., 2009, 5, 23-48\n"
"\n"
"\n"
"EXAMPLE\n"
"\n"
"block_average trj_1.dat 2 20 100\n"
"\n"
"In this case, trj_1.dat is the data file (I used NAMD's output of the box\n"
"dimensions), 2 means analyze the 2nd column (the dimension of the x \n"
"coordinate), 20 means use from 2--20 blocks, and 100 means skip the \n"
"first 100 time points in the file.\n"
"\n"
"The output will look like:\n"
"\n"
"# block_average 'trj_1.dat' '2' '20' '100' - alan (Fri Apr  6 14:22:45 2012) {/home/alan/projects/IBM/lipopeptides/analysis/box_area/pope_popg} [2.0.0 120406]\n"
"# Num_Blocks    BlockLen        StdErr\n"
"20              111             0.190357\n"
"19              117             0.194111\n"
"18              123             0.193382\n"
"17              131             0.19941\n"
"(more lines like this)\n"
"3               743             0.239527\n"
"2               1114            0.354243\n"
"\n"
"The first column is the number of blocks the data is broken into, the \n"
"second is the number of time points in each block, and the last column \n"
"is the standard error of the averages for each block.  As a rule, you'll\n"
"want to plot this data using column 2 as the x-axis, and column 3 as the \n"
"y-axis.\n"
"\n"
        ;
    return(s);
    }


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

if ( (argc >=2) && (strncmp(argv[1], "--fullhelp", 10) == 0) )
    {
    cout << fullHelpMessage() << endl;
    exit(-1);
    }

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
    int time = num_points / i;
    float variance = data.block_var(i);
    float std_err = sqrt(variance/i);
    cout << i << "\t\t"
         << time << "\t\t"
         << std_err
         << endl;
    }


}
