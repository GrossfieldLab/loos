/*
  Given 4 selections, compute the torsion angle for their centroids.  This 
  program loops over a trajectory and writes the torsion angle time series.

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
"    SYNOPSIS\n"
"\n"
"    Compute a time series of the torsion formed by 4 selections\n"
"\n"
"    DESCRIPTION\n"
"\n"
"    This tool loops over a trajectory, computing the torsion angle formed \n"
"    by the centroids of four selections.  \n"
"    \n"
"    Note: the code does not make any attempt to ensure that the entirety \n"
"    of a given selection is found within the same periodic image; if a \n"
"    selection is split (e.g. some of it is at the +x edge of the box and \n"
"    some at the -x edge), then the centroid is not a good description of the \n"
"    position.  However, if you're working with pieces of a protein and\n"
"    you've run the system through merge-traj with fix-imaging and centering,\n"
"    you will probably be fine.\n"
"\n"
"    EXAMPLE\n"
"\n"
"    torsion model.psf trajectory.dcd 'resid == 5' 'resid == 6' 'resid == 7' 'resid == 8'\n"
        ;
    return(s);
    }

void Usage()
    {
    cerr << "Usage: torsion system trajectory selection1 selection2 "
         << "selection3 selection4"
         << endl;
    }

int main (int argc, char *argv[])
{

if ( (argc >= 2) && (string(argv[1]) == string("--fullhelp") ) )
    {
    cerr << fullHelpMessage() << endl;
    exit(-1);
    }

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
char *selection3 = argv[5];  // String describing the third selection
char *selection4 = argv[6];  // String describing the fourth selection


AtomicGroup group1 = loos::selectAtoms(system, selection1);
AtomicGroup group2 = loos::selectAtoms(system, selection2);
AtomicGroup group3 = loos::selectAtoms(system, selection3);
AtomicGroup group4 = loos::selectAtoms(system, selection4);


// read the initial coordinates into the system
traj->updateGroupCoords(system);

// loop over the frames of the trajectory
int frame = 0;
double t;
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);
    t = Math::torsion(group1.centroid(), group2.centroid(),
                      group3.centroid(), group4.centroid());
    cout << frame << "\t" << t << endl;
    frame++;
    }
}

