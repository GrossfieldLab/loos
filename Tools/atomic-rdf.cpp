/*
  Compute 3d radial distribution function for two selections.
  This version treats each atom in each selection independently.  If you
  want to look at the distribution of centers of mass, look at rdf instead.

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



void Usage()
    {
    cerr << "Usage: atomic-rdf system trajectory selection1 selection2 "
         << "min max num_bins skip" 
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
cout << "# " << loos::invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names
// Create the system and read the trajectory file
AtomicGroup system = loos::createSystem(argv[1]);
pTraj traj = loos::createTrajectory(argv[2], system);

char *selection1 = argv[3];  // String describing the first selection
char *selection2 = argv[4];  // String describing the second selection
double hist_min = atof(argv[5]); // Lower edge of the histogram
double hist_max = atof(argv[6]); // Upper edge of the histogram
int num_bins = atoi(argv[7]); // Number of bins in the histogram
int skip = atoi(argv[8]);  // Number of frames to skip as equilibration

double bin_width = (hist_max - hist_min)/num_bins;


// Set up the selector to define group1 atoms
AtomicGroup group1 = loos::selectAtoms(system, selection1);

// Set up the selector to define group2 atoms
AtomicGroup group2 = loos::selectAtoms(system, selection2);

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
    // update coordinates and periodic box
    traj->updateGroupCoords(system);
    GCoord box = system.periodicBox(); 
    volume += box.x() * box.y() * box.z();

    unique_pairs = 0;
    // compute the distribution of g2 around g1 
    for (int j = 0; j < group1.size(); j++)
        {
        pAtom a1 = group1[j];
        GCoord p1 = a1->coords();
        for (int k = 0; k < group2.size(); k++)
            {
            pAtom a2 = group2[k];
            // skip "self" pairs 
            if (a1 == a2)
                {
                continue;
                }
            unique_pairs++;
            GCoord p2 = a2->coords();
            // Compute the distance squared, taking periodicity into account
            double d2 = p1.distance2(p2, box);
            if ( (d2 <= max2) && (d2 >= min2) )
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
    cum1 += hist[i] / (frame*group1.size());
    cum2 += hist[i] / (frame*group2.size());

    cout << d << "\t" << total << "\t" 
         << cum1 << "\t" << cum2 << endl;

    }
}

