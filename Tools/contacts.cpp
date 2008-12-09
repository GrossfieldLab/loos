/*
 *     Compute 2d radial distribution function for two selections
 */


/*
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
    cerr << "Usage: contacts model trajectory selection1 selection_2 max" 
         << endl;
    }

int main (int argc, char *argv[])
{
  if ( (argc <= 1) || 
       ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
       (argc < 6)
       )
    {
      Usage();
      exit(-1);
    }

  cout << "# " << invocationHeader(argc, argv) << endl;

  // copy the command line variables to real variable names
  string model_filename(argv[1]);
  string traj_filename(argv[2]);
  string selection1(argv[3]);
  string selection2(argv[4]);
  double max = strtod(argv[5], 0);
  double max2 = max*max;

  AtomicGroup model = createSystem(model_filename);
  pTraj traj = createTrajectory(traj_filename, model);


  // The assumption here is that selection1 will specify a bunch of molecules,
  // eg, all lipid headgroups of type foo.  g1 will be the group containing all
  // of those atoms.  However, what we'll actually want to do is work with the 
  // individual molecules' centers of mass, so we'll split g1 into a vector of
  // individual segids (corresponding to individual lipids) called group1
  AtomicGroup g1 = selectAtoms(model, selection1);
  vector<AtomicGroup> group1 = g1.splitByUniqueSegid();

  // g2 / group2 works the same way
  AtomicGroup g2 = selectAtoms(model, selection2);
  vector<AtomicGroup> group2 = g2.splitByUniqueSegid();



  cout << "#Frame\tPairs\tPerGroup1\tPerGroup2" << endl;

  // loop over the frames of the dcd file
  int frame = 0;
  while (traj->readFrame())
    {
      // get the new coordinates
      traj->updateGroupCoords(model);
      int count = 0;

      // compute the number of contacts between group1 center of mass 
      // and group2 center of mass
      vector<AtomicGroup>::iterator first;
      for (first=group1.begin(); first!=group1.end(); first++)
        {
      GCoord com1 = first->centerOfMass();

      vector<AtomicGroup>::iterator second;
      for (second=group2.begin(); second!=group2.end(); second++)
            {
          // exclude self pairs 
          if (*first == *second)
                continue;

          GCoord com2 = second->centerOfMass();

          double d2 = com1.distance2(com2, model.periodicBox());
          if ( (d2 <= max2) )
                {
          count++;
                }
            }
        }
    
      // Output the results
      double per_g1_atom = (double)count / group1.size();
      double per_g2_atom = (double)count / group2.size();
      cout << frame << "\t" 
       << count << "\t"
       << per_g1_atom << "\t"
       << per_g2_atom << endl;

      frame++;
    }
}

