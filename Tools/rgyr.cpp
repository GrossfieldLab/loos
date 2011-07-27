/*
  Compute the distribution of radii of gyration for a selection of atoms.
  

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
    cerr << "Usage: rgyr SystemFile Trajectory selection "
         << "min max num_bins skip [by-molecule]" 
         << endl;
    cerr << "by-molecule should be one if you want the selection "
         << "broken up based on connectivity, and 0 or absent otherwise."
         << endl;
    }

int main (int argc, char *argv[])
{
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 8)
   )
    {
    Usage();
    exit(-1);
    }

// Print the command line arguments
cout << "# " << invocationHeader(argc, argv) << endl;

// Create the system and the trajectory file
AtomicGroup system = createSystem(argv[1]);
pTraj traj = createTrajectory(argv[2], system);

char *selection = argv[3];  // String describing the first selection

// Get the histogram parameters
greal hist_min = atof(argv[4]); // Lower edge of the histogram
greal hist_max = atof(argv[5]); // Upper edge of the histogram
int num_bins = atoi(argv[6]); // Number of bins in the histogram

// Get the number of frames to discard as equilibration
int skip = atoi(argv[7]);  

int split_by_molecule = 0;
if (argc >= 9)
    {
    split_by_molecule = atoi(argv[8]);
    }

greal bin_width = (hist_max - hist_min)/num_bins;



vector<AtomicGroup> molecules;
if (split_by_molecule)
    {
    molecules = system.splitByMolecule();
    }
else
    {
    molecules.push_back(system);
    }

// Set up the selector to define the selected group
Parser parser(selection);
KernelSelector parsed_sel(parser.kernel());


// Loop over the molecules and add them to selection
vector<AtomicGroup> molecule_groups;
vector<AtomicGroup>::iterator m;
for (m=molecules.begin(); m!=molecules.end(); m++)
    {
    AtomicGroup tmp = m->select(parsed_sel);
    if (tmp.size() > 0)
        {
        molecule_groups.push_back(tmp);
        }
    }


// Skip the initial frames as equilibration
traj->readFrame(skip); 

// read the initial coordinates into the system
traj->updateGroupCoords(system);

// Create the histogram and zero it out
vector<greal> hist;
hist.reserve(num_bins);
hist.insert(hist.begin(), num_bins, 0.0);

// loop over the frames of the trajectory
int frame = 0;
int count = 0;
while (traj->readFrame())
    {
    // update coordinates and periodic box
    traj->updateGroupCoords(system);

    vector<AtomicGroup>::iterator m;
    for (m=molecule_groups.begin(); m!=molecule_groups.end(); m++)
        {
        greal rad = m->radiusOfGyration();
        if ( (rad >=hist_min) && (rad <hist_max) )
            {
            int bin = int((rad-hist_min)/bin_width);
            hist[bin]++;
            count++;
            }
        }
    frame++;
    }


// Output the results
cout << "# Rgyr\tProb\tCum" << endl;
greal cum = 0.0;
for (int i = 0; i < num_bins; i++)
    {
    greal d = bin_width*(i + 0.5) + hist_min;


    greal prob = hist[i]/ count;
    cum += prob;
    cout << d << "\t" 
         << prob << "\t" 
         << cum << endl;

    }
}

