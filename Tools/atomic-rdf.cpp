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
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;



void Usage()
    {
    cerr << "Usage: atomic-rdf system trajectory selection1 selection2 "
         << "min max num_bins skip" 
         << endl;
    }

int main (int argc, char *argv[])
{

// Build options
opts::BasicOptions* bopts = new opts::BasicOptions;
opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
opts::RequiredArguments* ropts = new opts::RequiredArguments;

// These are required command-line arguments (non-optional options)
ropts->addArgument("selection1", "selection1");
ropts->addArgument("selection2", "selection2");
ropts->addArgument("min", "min radius");
ropts->addArgument("max", "max radius");
ropts->addArgument("num_bins", "number of bins");

opts::AggregateOptions options;
options.add(bopts).add(tropts).add(ropts);
if (!options.parse(argc, argv))
  exit(-1);


// Print the command line arguments
cout << "# " << invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names
// Create the system and read the trajectory file
AtomicGroup system = tropts->model;
pTraj traj = tropts->trajectory;
if (!(system.isPeriodic() || traj->hasPeriodicBox()))
  {
  cerr << "Error- Either the model or the trajectory must have periodic box information.\n";
  exit(-1);
  }

// Extract our required command-line arguments
string selection1 = ropts->value("selection1");  // String describing the first selection
string selection2 = ropts->value("selection2");  // String describing the second selection
double hist_min = parseStringAs<double>(ropts->value("min")); // Lower edge of the histogram
double hist_max = parseStringAs<double>(ropts->value("max")); // Upper edge of the histogram
int num_bins = parseStringAs<double>(ropts->value("num_bins")); // Number of bins in the histogram


double bin_width = (hist_max - hist_min)/num_bins;


// Set up the selector to define group1 atoms
AtomicGroup group1 = selectAtoms(system, selection1);

// Set up the selector to define group2 atoms
AtomicGroup group2 = selectAtoms(system, selection2);

// Create the histogram and zero it out
vector<double> hist;
hist.reserve(num_bins);
hist.insert(hist.begin(), num_bins, 0.0);

double min2 = hist_min*hist_min;
double max2 = hist_max*hist_max;

// loop over the frames of the trajectory
int frame = 0;
double volume = 0.0;
unsigned long unique_pairs=0;
while (traj->readFrame())
    {

    // update coordinates and periodic box
    traj->updateGroupCoords(system);
    GCoord box = system.periodicBox(); 
    volume += box.x() * box.y() * box.z();

    unique_pairs = 0;
    // compute the distribution of g2 around g1 
    for (uint j = 0; j < group1.size(); j++)
        {
        pAtom a1 = group1[j];
        GCoord p1 = a1->coords();
        for (uint k = 0; k < group2.size(); k++)
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
    cum1 += hist[i] / (frame*group1.size());
    cum2 += hist[i] / (frame*group2.size());

    cout << d << "\t" << total << "\t" 
         << cum1 << "\t" << cum2 << endl;

    }
}

