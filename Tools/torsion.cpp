/*
   Given a selection of 4 atoms and a trajectory, compute the torsion time 
   series.

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
    cerr << "Usage: torsion system trajectory selection1 selection2 "
         << "selection3 selection4"
         << endl;
    }

int main (int argc, char *argv[])
{
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 7)
   )
    {
    Usage();
    exit(-1);
    }

// Print the command line arguments
cout << "# " << invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names
// Create the system and read the trajectory file
AtomicGroup system = createSystem(argv[1]);
pTraj traj = createTrajectory(argv[2], system);

char *selection1 = argv[3];  // String describing the first selection
char *selection2 = argv[4];  // String describing the second selection
char *selection3 = argv[5];  // String describing the first selection
char *selection4 = argv[6];  // String describing the second selection


AtomicGroup group1 = loos::selectAtoms(system, selection1);
AtomicGroup group2 = loos::selectAtoms(system, selection2);
AtomicGroup group3 = loos::selectAtoms(system, selection3);
AtomicGroup group4 = loos::selectAtoms(system, selection4);


if ( (group1.size() != 1) ||
     (group2.size() != 1) ||
     (group3.size() != 1) || 
     (group4.size() != 1) )
    {
    cerr << "Each selection must contain exactly 1 atom" << endl;
    exit(-1);
    }

pAtom a1 = group1[0];
pAtom a2 = group2[0];
pAtom a3 = group3[0];
pAtom a4 = group4[0];

// read the initial coordinates into the system
traj->updateGroupCoords(system);

// loop over the frames of the trajectory
int frame = 0;
double t;
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);
    t = Math::torsion(a1,a2,a3,a4);
    cout << frame << "\t" << t << endl;
    frame++;
    }
}

