/*
  reimage-by-molecule
  
  Read a model and trajectory, reimage each frame by molecule, and write a new
  dcd

  Alan Grossfield
 
  
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

const int update_frequency = 250;

using namespace std;
using namespace loos;


void Usage()
{
  cerr << "Usage: reimage-by-molecule model trajectory outdcd [xbox ybox zbox]"
       << endl;
} 

int main(int argc, char *argv[]) 
{
  if ( (argc <= 1) ||
       ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
       (argc < 4)
       )
    {
      Usage();
      exit(-1);
    }

  string hdr = invocationHeader(argc, argv);

  AtomicGroup model = createSystem(argv[1]);
  if (!model.hasBonds())
    {
      cerr << "***WARNING***\nThe model does not have connectivity,\nso your results may not be what you expect.\n";
    }

  pTraj traj = createTrajectory(argv[2], model);
  DCDWriter dcd_out(argv[3]);

  bool box_override = false;
  GCoord newbox;
    
  if (argc == 7) {
    // use the supplied box, set it in model
    double xbox = atof(argv[4]);
    double ybox = atof(argv[5]);
    double zbox = atof(argv[6]);
    newbox = GCoord(xbox, ybox, zbox);
    model.periodicBox(newbox);
    box_override = true;
    if (!traj->hasPeriodicBox())
      {
        cerr << "Adding box " << newbox << endl;
      }
    else
      {
        cerr << "WARNING - Overriding existing box(es) with " << newbox << endl;
      }
  }

  if (!traj->hasPeriodicBox() && !box_override)
    {
      // if the trajectory doesn't have periodic box info, we need it from
      // the user
      cerr << "ERROR - The trajectory has no box information.  You must add it or supply it on the command-line.\n";
      exit(-1);
    }

  // duplicate the info from dcd in dcd_out
  dcd_out.setHeader(traj->natoms(), traj->nframes(), traj->timestep(), 
                    true );
  dcd_out.setTitle(hdr);
  dcd_out.writeHeader();

  // split the system by molecule
  vector<AtomicGroup> molecules = model.splitByMolecule();
  cerr << "Found " << molecules.size() << " molecules.\n";

  // split the system by segid
  vector<AtomicGroup> segments = model.splitByUniqueSegid();
  cerr << "Found " << segments.size() << " segments.\n";

  cerr << "Trajectory has " << traj->nframes() << " total frames.\n";


  // Loop over the frames of the dcd and reimage each molecule
  vector<AtomicGroup>::iterator m;
  int frame_no = 0;
  cerr << "Frames processed - ";
  while (traj->readFrame())
    {
      if (++frame_no % update_frequency == 0)
        {
          cerr << frame_no << " ";
        }

      traj->updateGroupCoords(model);

      if (box_override)
        model.periodicBox(newbox);

      for (m = segments.begin(); m != segments.end(); m++)
        {
          m->reimage();
        }
      for (m = molecules.begin(); m != molecules.end(); m++)
        {
          m->reimage();
        }

      dcd_out.writeFrame(model);
    }

  cerr << " - done\n";

}

