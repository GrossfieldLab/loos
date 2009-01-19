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

if (argc != 5)
    {
    cerr << "Usage: recenter-trj model-file trajectory-file selection-string dcd-name" << endl;
    exit(-1);
    }

AtomicGroup model = createSystem(argv[1]);
pTraj traj = createTrajectory(argv[2], model);
AtomicGroup center = selectAtoms(model, argv[3]);

vector<AtomicGroup> molecules= model.splitByMolecule();

DCDWriter dcd(argv[4]);
//dcd.setHeader(model.size(), traj->nframes(), 1e-3, traj->hasPeriodicBox());
//dcd.setTitle(invocationHeader(argc, argv));

vector<AtomicGroup>::iterator m;

while (traj->readFrame())
    {
    traj->updateGroupCoords(model); 
    GCoord centroid = center.centroid();
    model.translate(-centroid);
    for (m=molecules.begin(); m!=molecules.end(); m++)
        {
        m->reimage();
        }
    
    dcd.writeFrame(model);
    }

}
