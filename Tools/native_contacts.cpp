/*
  native_contacts: compute the fraction of native contacts in a trajectory,
                   based on an initial structure.

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
    cerr << "Usage: native_contacts system traj cutoff selection "
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

cout << "# " << invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names
AtomicGroup system = createSystem(argv[1]);
pTraj traj = createTrajectory(argv[2], system);
double cutoff = atof(argv[3]);
char *selection = argv[4];

double cut2 = cutoff*cutoff;


AtomicGroup sel = selectAtoms(system, selection);
vector<AtomicGroup> residues = sel.splitByResidue();

// Check to see if the system has coordinates -- we'll use them if it does,
// otherwise we'll use the first frame of the trajectory as the reference 
// structure.
if ( !(sel[0]->checkProperty(Atom::coordsbit)) )
    {
    traj->readFrame(0);
    traj->updateGroupCoords(system);
    }

// Compute the centers of mass of the selections
uint num_residues = residues.size();
vector<GCoord> centers_of_mass(num_residues);
for (uint i=0; i<num_residues-1; i++)
    {
    centers_of_mass[i] = residues[i].centerOfMass();
    }

vector<vector<uint> > contacts;
// Find contacts within the threshold distance
for (uint i=0; i<num_residues-1; i++)
    {
    for (uint j=i+1; j< num_residues; j++)
        {
        GCoord diff = centers_of_mass[j] - centers_of_mass[i]; 
        if (diff.length2() <= cut2)
            {
            vector<uint> v(2);
            v[0] = i;
            v[1] = j;
            contacts.push_back(v);
            // Print out the contacts as we go
            cout << "# " << i << "\t" << j << endl;
            }
        }
    }

// Number of native contacts, as a float because we'll need
// to do floating point arithmatic with it later anyway
float num_native_contacts = (float) contacts.size();
cout << "# Total native contacts: " << num_native_contacts << endl;

// Loop over structures in the trajectory
vector<vector<uint> >::iterator p;
int frame = 0;
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);

    // Loop over contacts from the native structure
    int num_contacts = 0;
    for (p=contacts.begin(); p!= contacts.end(); p++)
        {
        uint r1 = p->at(0);
        uint r2 = p->at(1);
        GCoord c1 = residues[r1].centerOfMass();
        GCoord c2 = residues[r2].centerOfMass();
        GCoord diff = c2 - c1;
        if (diff.length2() <= cut2)
            {
            num_contacts++;
            }
        }
    float fraction = num_contacts / num_native_contacts;
    cout << frame << "\t" << fraction << endl;
    frame++;
    }

}

