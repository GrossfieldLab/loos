/*

   recenter-trj

   Reads a LOOS-supported trajectory format and a selection, writes a 
   dcd with the selection centered at the origin and the rest of the system
   recentered by molecule.

   Usage: recenter-trj model-file trajectory-file selection-string dcd-name
 
*/

/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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

int main(int argc, char *argv[])
{

if (argc != 6)
    {
    cerr << "Usage: recenter-trj model-file trajectory-file selection-string [Z|XY|A] dcd-name" << endl;
    exit(-1);
    }

AtomicGroup model = createSystem(argv[1]);
pTraj traj = createTrajectory(argv[2], model);
AtomicGroup center = selectAtoms(model, argv[3]);
string flag = string(argv[4]);
bool just_z = false;
bool just_xy = false;
if ( (flag == "Z") || (flag == "z"))
    {
    just_z = true;
    }
else if ( (flag == "XY") || (flag == "xy"))
    {
    just_xy = true;
    }

DCDWriter dcd(argv[5]);
dcd.setTitle(invocationHeader(argc, argv));

if (!model.hasBonds())
    {
    cerr << "Error: " << argv[0]
         << " will only work if the system has connectivity information."
         << endl
         << "You'll need to use something like a PSF or PDB with conect records"
         << endl;
    exit(-1);
    }

vector<AtomicGroup> molecules= model.splitByMolecule();
vector<AtomicGroup>::iterator m;

while (traj->readFrame())
    {
    traj->updateGroupCoords(model); 

    // Simple approach won't work if the centering selection is split
    // across the periodic image.  In that case, the centroid may be near the 
    // middle even if none of the atoms are near there.

    // pick a single atom in the selection, and center based on it.
    // This will make sure the selection is now _not_ split acrosst the 
    // periodic image
    GCoord centroid = center[0]->coords();
    if (just_z)
        {
        centroid.x() = 0.0;
        centroid.y() = 0.0;
        }
    else if (just_xy)
        {
        centroid.z() = 0.0;
        }

    model.translate(-centroid);
    for (m=molecules.begin(); m!=molecules.end(); m++)
        {
        m->reimage();
        }
    
    // now, center as we did in the original algorithm:
    // Move the whole system such that selected region is at the origin and
    // reimage
    centroid = center.centroid();
    if (just_z)
        {
        centroid.x() = 0.0;
        centroid.y() = 0.0;
        }
    else if (just_xy)
        {
        centroid.z() = 0.0;
        }

    model.translate(-centroid);
    for (m=molecules.begin(); m!=molecules.end(); m++)
        {
        m->reimage();
        }
    
    dcd.writeFrame(model);
    }

}
