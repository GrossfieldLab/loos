
/*
  Determine which and how many molecules are "bound" to the lipid, defined
  by the position of the z center of mass relative to user provided
  boundaries of the lipid membrane.
  
  Joshua Horn
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
    cerr << "Usage: bound SystemFile Trajectory selection "
         << "skip lastframe boundary [by-molecule]" 
         << endl;
    cerr << "by-molecule should be one if you want the selection "
         << "broken up based on connectivity, and 0 or absent otherwise."
         << endl;
    cerr << "Set lastframe to 0 to include entire trajectory."
         << endl;
     }

int main (int argc, char *argv[])
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

// Create the system and the trajectory file
AtomicGroup system = createSystem(argv[1]);
pTraj traj = createTrajectory(argv[2], system);

// String describing the first selection
char *selection = argv[3];

// Get the number of frames to discard as equilibration
int skip = atoi(argv[4]);

// Get the last frame to analyze
int lastFrame = atoi(argv[5]);
if (lastFrame == 0)
    {

        lastFrame = traj->nframes();

    }

// Set boundaries that define bound state
double boundary = strtod(argv[6], NULL);

// Split by molecule
int split_by_molecule = 0;
if (argc >= 8)
    {

    split_by_molecule = atoi(argv[7]);
    
    }

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


// Print headers for output
cout << "#Time\tAvg     \tStdev    \t";

// Loop over the molecules and add them to selection
vector<AtomicGroup> molecule_groups;
vector<AtomicGroup>::iterator m;

// Build atomic groups and push them into molecule_groups
for (m=molecules.begin(); m!=molecules.end(); m++)
    {
    AtomicGroup tmp = m->select(parsed_sel);
    if (tmp.size() > 0)
        {
        molecule_groups.push_back(tmp);

        // Retrieve segid of each molecule and print to output
            //cout << (tmp.getAtom(1))->segid() << "\t";
        }

    }

// Print endline for formatting
cout << "Total" << endl;

// Skip the initial frames as equilibration
traj->readFrame(skip); 

// read the initial coordinates into the system
traj->updateGroupCoords(system);

// loop over the frames of the trajectory
int frame = 0;

// Create a vector and a timeseries to store data
vector<double> series;
TimeSeries<double> tseries;

// Counter to keep track of total bound molecules
int count = 0;

while (traj->readFrame() && frame < lastFrame)
    {

    // update coordinates and periodic box
    traj->updateGroupCoords(system);

    // iterate through all the molecules
    vector<AtomicGroup>::iterator m;
    for (m=molecule_groups.begin(); m!=molecule_groups.end(); m++)
        {

        // Add current molecule's center of mass to vector
        series.push_back(fabs((m->centerOfMass()).z()));
        }
    
    // Print the frame to output
    cout << frame << "\t";

    // Wrap the vector in a time series and print avg, stdev
    tseries = TimeSeries<double>(series);
    cout.setf(ios::fixed);
    cout << setprecision(6) << tseries.average() << "\t";
    cout << setprecision(6) << tseries.stdev() << "\t";

    // For each in the series, print 0 (unbound) if less than
    // given boundary value and print 1 (bound) if more
    for (int j = 0; j < (int)series.size(); j++)
        {

	    if (series.at(j) > boundary || series.at(j) < -boundary)
            {
        
            //cout << "0\t";
		
		    } else {
		
            //cout << "1\t";
		    count++;
		
            }

        }
    
    // Clear the vector, print the total and reset the counter
    series.clear();
    //cout << count << endl;
	cout << endl;
    count = 0;

    // Increments the frame in the trajectory
    frame++;

    }

}

