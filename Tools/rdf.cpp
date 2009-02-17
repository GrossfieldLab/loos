/*
  Compute 3d radial distribution function for two selections
  This version works on groups of atoms, not individual atoms, meaning that
  the selections on the command line are divided up by molecule, and the 
  per-molecule center of mass is used.  If that's not what you want, you may 
  want to take a look at atomic-rdf.

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
    cerr << "Usage: rdf SystemFile Trajectory selection1 selection2 "
         << "min max num_bins skip" 
         << endl;
    cerr << endl;
    cerr << "SystemFile can be a CHARMM/XPLOR PSF, an AMBER parmtop, a PDB"
         << " file, or any other format supported by LOOS"
         << endl;
    cerr << "Trajectory can be a CHARMM/NAMD DCD or AMBER trajectory," 
         << " or any other format supported by LOOS"
         << endl;
    }

int main (int argc, char *argv[])
{
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 9)
   )
    {
    Usage();
    exit(-1);
    }

// Print the command line arguments
cout << "# " << invocationHeader(argc, argv) << endl;

// Create the system and the trajectory file
// Note: The pTraj type is a Boost shared pointer, so we'll need
//       to use pointer semantics to access it
AtomicGroup system = createSystem(argv[1]);
pTraj traj = createTrajectory(argv[2], system);

// Get the selection patterns to identify the atoms to be used
char *selection1 = argv[3];  // String describing the first selection
char *selection2 = argv[4];  // String describing the second selection

// Get the histogram parameters
double hist_min = atof(argv[5]); // Lower edge of the histogram
double hist_max = atof(argv[6]); // Upper edge of the histogram
int num_bins = atoi(argv[7]); // Number of bins in the histogram

// Get the number of frames to discard as equilibration
int skip = atoi(argv[8]);  

double bin_width = (hist_max - hist_min)/num_bins;

// The groups may describe a set of individual atoms (eg water oxygens) 
// or chunks (eg lipid headgroups).  What we're doing here is dividing up
// each group by molecule, so that (for example) if selection1 was all of the 
// water, we'd break group1 up into pieces with each water being it's own piece.
// The RDF would then use the center of mass of each water molecule.
// The problem is that you can only do split by molecule when all of the atoms
// are present (so you can walk the connectivity tree correctly), what we'll
// do is first split the atoms by molecules, then rejoin the ones which match
// each of the patterns.
vector<AtomicGroup> molecules = system.splitByMolecule();


// Set up the selector to define group1 atoms
Parser parser1(selection1);
KernelSelector parsed_sel1(parser1.kernel());

// Set up the selector to define group2 atoms
Parser parser2(selection2);
KernelSelector parsed_sel2(parser2.kernel());

// Loop over the molecules and add them to group1 or group2
vector<AtomicGroup> g1_mols, g2_mols;
vector<AtomicGroup>::iterator m;
for (m=molecules.begin(); m!=molecules.end(); m++)
    {
    AtomicGroup tmp = m->select(parsed_sel1);
    if (tmp.size() > 0)
        {
        g1_mols.push_back(tmp);
        }

    AtomicGroup tmp2 = m->select(parsed_sel2);
    if (tmp2.size() > 0)
        {
        g2_mols.push_back(tmp2);
        }
    }


// Skip the initial frames as equilibration
traj->readFrame(skip); 

// read the initial coordinates into the system
traj->updateGroupCoords(system);

// Create the histogram and zero it out
vector<double> hist;
hist.reserve(num_bins);
hist.insert(hist.begin(), num_bins, 0.0);

double min2 = hist_min*hist_min;
double max2 = hist_max*hist_max;

// loop over the frames of the trajectory
int frame = 0;
double volume = 0.0;
int unique_pairs=0;
while (traj->readFrame())
    {
    cerr << "Processing frame " << frame << endl;
    // update coordinates and periodic box
    traj->updateGroupCoords(system);
    GCoord box = system.periodicBox(); 
    volume += box.x() * box.y() * box.z();

    unique_pairs = 0;

    // compute the distribution of g2 around g1 
    for (unsigned int j = 0; j < g1_mols.size(); j++)
        {
        GCoord p1 = g1_mols[j].centerOfMass();
        for (unsigned int k = 0; k < g2_mols.size(); k++)
            {
            // skip "self" pairs -- in case selection1 and selection2 overlap
            if (g1_mols[j] == g2_mols[k])
                {
                continue;
                }
            unique_pairs++;
            GCoord p2 = g2_mols[k].centerOfMass();
            
            // Compute the distance squared, taking periodicity into account
            double d2 = p1.distance2(p2, box);
            if ( (d2 < max2) && (d2 > min2) )
                {
                double d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist[bin]++;
                }
            }
        }
    frame++;
    }

volume /= frame;



double expected = frame * unique_pairs / volume;
double cum1 = 0.0;
double cum2 = 0.0;

// Output the results
cout << "# Dist\tRDF\tCumAround1\tCumAround2" << endl;
for (int i = 0; i < num_bins; i++)
    {
    double d = bin_width*(i + 0.5);

    double d_inner = bin_width*i;
    double d_outer = d_inner + bin_width;
    double norm = 4.0/3.0 * M_PI*(d_outer*d_outer*d_outer 
                                - d_inner*d_inner*d_inner);

    double total = hist[i]/ (norm*expected);
    cum1 += hist[i] / (frame*g1_mols.size());
    cum2 += hist[i] / (frame*g2_mols.size());

    cout << d << "\t" << total << "\t" 
         << cum1 << "\t" << cum2 << endl;

    }
}

