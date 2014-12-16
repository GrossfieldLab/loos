/*

   recenter-trj

   Reads a LOOS-supported trajectory format and a selection, writes a 
   dcd with the selection centered at the origin and the rest of the system
   recentered by molecule.
 
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

string fullHelpMessage(void)
    {
    string s = 
"\n"
"SYNOPSIS\n"
"\n"
"Read a trajectory and produce a new trajectory with the selected set of\n"
"atoms recentered.\n"
"\n"
"DESCRIPTION\n"
"\n"
"This program translates and reimages a trajectory such that the set of atoms\n"
"selected is recentered.  Some of the capabilities are redundant with \n"
"the merge-traj tool, and at some point they may be merged.  However, the\n"
"main unique use of this tool is the ability to recenter just in the x-y\n"
"plane, just along the z-axis, or in all 3 dimensions at once. No rotations\n"
"are performed.\n"
"\n"
"Unlike merge-traj, recenter-trj always handles the case where the centering\n"
"selection might be split across the periodic boundary, and so does not\n"
"need a flag like --selection-is-split.  \n"
"\n"
"recenter-trj will only work when the system file specifies the system's \n"
"connectivity, as with a CHARMM/NAMD psf file, or a PDB file with CONECT \n"
"records.\n"
"\n"
"EXAMPLE\n"
"\n"
"recenter-trj model.psf traj.dcd 'segname == \"PROT\"' A output.dcd\n"
"\n"
"Here, model.psf is the system file, traj.dcd is the input trajectory file,\n"
"and the selection string specifies a segment called PROT, presumably a \n"
"protein molecule.  The \"A\" argument means that the selection\n"
"is centered in all 3 dimensions.  \n"
    ;
    return(s);
    }

int main(int argc, char *argv[])
{

if ((argc > 1) && (string(argv[1]) == string("--fullhelp")))
    {
    cerr << fullHelpMessage() << endl;
    exit(-1);
    }
else if (argc != 6)
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


pTrajectoryWriter traj_out = createOutputTrajectory(argv[5]);
traj_out->setComments(invocationHeader(argc, argv));

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
    
    traj_out->writeFrame(model);
    }

}
