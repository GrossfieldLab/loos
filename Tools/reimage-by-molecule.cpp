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

string fullHelpMessage(void)
    {
    string s =
"\n"
"    SYNOPSIS\n"
"\n"
"    Read a trajectory and reimage it such that each molecule has its\n"
"    centroid in the central box.\n"
"\n"
"    DESCRIPTION\n"
"\n"
"    This tool reads a trajectory and processes it to produce a new \n"
"    trajectory in DCD format where each molecule has its centroid in\n"
"    the central image.\n"
"    \n"
"    This operation does not make a lot of sense if the system file \n"
"    does not contain connectivity information; it will warn you \n"
"    if you invoke it without connectivity, but will run.\n"
"\n"
"    If the trajectory has information on box size built in to it, that\n"
"    box data is used for the reimaging.  If not, the periodicity information\n"
"    may be read from the model file (e.g. a CRYSTL line from a PDB file).\n"
"    Alternatively, the user can provide box size information on the command \n"
"    line by supplying 3 extra arguments.  If this is done, the information\n"
"    overrides anything supplied in the trajectory or model file.\n"
"\n"
"    Note: this tool is largely redundant with merge-traj and recenter-traj \n"
"          (which also have additional capabilities), and may at some point \n"
"          be deprecated.\n"
"\n"
"    EXAMPLE\n"
"\n"
"    reimage-by-molecule model.psf input_traj.dcd output_traj.dcd \n"
"\n"
"    This reads the system information from model.psf, operates on \n"
"    input_traj.dcd (which presumably has periodicity information), and \n"
"    writes output_traj.dcd, which does have periodicity information.\n"
"\n"
"    reimage-by-molecule model.psf input_traj.dcd output_traj.xtc 55 77 100\n"
"\n"
"    This does essentially the same thing, but asserts that the periodic\n"
"    box is constant with x-dimension 55 angstrom, y-dimension 77 angstroms,\n"
"    and z-dimension 100 angstroms.  The trajectory is also converted to\n"
"    the GROMACS XTC format."
"\n"
        ;
    return(s);
    }


void Usage()
{
  cerr << "Usage: reimage-by-molecule model trajectory output-trajectory [xbox ybox zbox]"
       << endl;
} 

int main(int argc, char *argv[]) 
{
  if ( (argc >=2) )
    {
    if (string(argv[1]) == string("-h"))
        {
        Usage();
        exit(-1);
        }
    else if (string(argv[1]) == string("--fullhelp"))
        {
        cerr << fullHelpMessage() << endl;
        Usage();
        exit(-1);
        }
    }
  if (argc < 4)
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
  pTrajectoryWriter traj_out = createOutputTrajectory(argv[3]);

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
  traj_out->setComments(hdr);

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

      traj_out->writeFrame(model);
    }

  cerr << " - done\n";

}

