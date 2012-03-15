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

string fullHelpMessage(void)
    {
    string s =
   "\n"
   "SYNOPSIS\n"
   "\n"
   "Count the number of contacts between the centers of mass of two sets\n"
   "of selections.\n"
   "\n"
   "DESCRIPTION\n"
   "\n"
   "This tool counts the number of contacts between two selections.\n"
   "Each selection is split by unique segment name, and the various \n"
   "segments are treated separately, using their centers of mass.  \n"
   "\n"
   "This tool provides a subset of the functionality supplied by rdf;\n"
   "if you need splitting by something other than segment, you're better off\n"
   "using rdf and looking at the cumulative columns, which have equivalent \n"
   "information.  The only advantage to using this tool is that it\n"
   "avoids taking the square root in the distance calculation, so \n"
   "it might be a little bit faster.\n"
   "\n"
   "EXAMPLE\n"
   "\n"
   "contacts model.pdb traj.dcd 'segname ==\"RHOD\"' 'segname =~\"^L[0-9]+\"' 18\n"
   "\n"
   "This command line reads model.pdb, loops over the trajectory traj.dcd, \n"
   "and looks at 2 selections.  The first is segment RHOD, which is \n"
   "the protein rhodopsin, while the second is a set of lipid molecules \n"
   "with segment names L1, L2, etc.  It'll report the time series of the \n"
   "number of lipids with centers of mass within 18 angstroms of the center \n"
   "of mass of the protein.  It will also report the same data normalized \n"
   "by the number of groups in the first and second selection, respectively.\n"
   "In this case, that means that since there's 1 protein, the second and \n"
   "third columns will be the same, while the fourth column will be the \n"
   "second column divided by the number of lipids selected.\n"
   "\n"
        ;
    return(s);
    }



void Usage()
    {
    cerr << "Usage: contacts model trajectory selection1 selection_2 max" 
         << endl;
    }

int main (int argc, char *argv[])
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
        exit(-1);
        }
    }
if (argc < 6)
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

