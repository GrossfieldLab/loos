/*
  Calculate the area per lipid if the selection is the specific lipid leaflet
  you wish to analyze.
  
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
    cerr << "Usage: lipid_area SystemFile Trajectory selection "
         << "skip lastframe" 
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

vector<AtomicGroup> molecules;
molecules = system.splitByMolecule();
  
// Set up the selector to define the selected group
Parser parser(selection);
KernelSelector parsed_sel(parser.kernel());


// Print headers for output
cout << "#Time\tArea per molecule";

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
        }

    }

// Skip the initial frames as equilibration
traj->readFrame(skip); 

// read the initial coordinates into the system
traj->updateGroupCoords(system);

// loop over the frames of the trajectory
int frame = 0;

int num_mol;
double X;
double Y;
double area;
double area_per_mol;

// Loop through trajectory, getting area per lipid
while (traj->readFrame() && frame < lastFrame)
    {

    // update coordinates and periodic box
    traj->updateGroupCoords(system);
    
    // Print the frame to output
    cout << frame << "\t";

    // Get the x and y coords, then get the area
    //
    X = (traj->periodicBox()).x();
    Y = (traj->periodicBox()).y();
    area = X * Y;

    // Get the number of molecules
    //
    num_mol = molecule_groups.size();

    // Divide area by number of molecules, thus the area per mol
    area_per_mol = area / num_mol;

    // PRINT AREA PER LIPID
    //
    cout << area_per_mol << endl;

    // Increments the frame in the trajectory
    frame++;

    }

}

