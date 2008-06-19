/*
  Compute 2d radial distribution function for two selections

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
using namespace std;

#include <iostream>
#include <math.h>

#include <loos.hpp>
#include <psf.hpp>
#include <dcd.hpp>
#include <utils.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>

void Usage()
    {
    cerr << "Usage: xy_rdf PSF DCD selection1 selection2 "
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

cout << "# " << invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names
PSF psf(argv[1]);
DCD dcd(argv[2]);
char *selection1 = argv[3];
char *selection2 = argv[4];
double hist_min = atof(argv[5]);
double hist_max = atof(argv[6]);
int num_bins = atoi(argv[7]);
int skip = atoi(argv[8]);

double bin_width = (hist_max - hist_min)/num_bins;

Parser parser1(selection1);
KernelSelector parsed_sel1(parser1.kernel());
AtomicGroup group1 = psf.select(parsed_sel1);

Parser parser2(selection2);
KernelSelector parsed_sel2(parser2.kernel());
AtomicGroup group2 = psf.select(parsed_sel2);

// The groups are presumed to be units like lipid headgroups.  We want to 
// treat each headgroup as a single unit.  This could reasonably be done
// with splitByUniqueSegid() or splitByMolecule(), depending on the details.
// For lipids and lipopeptides, it's the same (although splitByMolecule will 
// have more overhead), but if you want water, then you need splitbyMolecule.
// If you get weird-looking behavior, you might want to change this to use
// splitByMolecule.
vector<AtomicGroup> g1_mols = group1.splitByUniqueSegid();
vector<AtomicGroup> g2_mols = group2.splitByUniqueSegid();

// Skip the initial frames
dcd.readFrame(skip); 

// read the initial coordinates into the psf
dcd.updateGroupCoords(psf);

// Now that we have some real coordinates, we need to subdivide the groups
// one more time, into upper and lower leaflets. This assumes that the 
// coordinates are properly centered and imaged.
vector<AtomicGroup> g1_upper, g1_lower;
vector<AtomicGroup> g2_upper, g2_lower;

for (int i = 0; i < g1_mols.size(); i++)
    {
    GCoord c = g1_mols[i].centerOfMass();
    if (c.z() >=0.0)
        {
        g1_upper.push_back(g1_mols[i]);
        }
    else
        {
        g1_lower.push_back(g1_mols[i]);
        }
    }

for (int i = 0; i < g2_mols.size(); i++)
    {
    GCoord c = g2_mols[i].centerOfMass();
    if (c.z() >=0.0)
        {
        g2_upper.push_back(g2_mols[i]);
        }
    else
        {
        g2_lower.push_back(g2_mols[i]);
        }
    }


// Create 2 histograms -- one for top, one for bottom
vector<double> hist_lower, hist_upper;
hist_lower.reserve(num_bins);
hist_upper.reserve(num_bins);
hist_lower.insert(hist_lower.begin(), num_bins, 0.0);
hist_upper.insert(hist_upper.begin(), num_bins, 0.0);

double min2 = hist_min*hist_min;
double max2 = hist_max*hist_max;

// loop over the frames of the dcd file
int frame = 0;
double area = 0.0;
while (dcd.readFrame())
    {
    // update coordinates and periodic box
    dcd.updateGroupCoords(psf);
    GCoord box = psf.periodicBox(); 
    area += box.x() * box.y();

#if 0
    if (frame % 10 == 0)
        {
        cout << "#Processing file " << file_list[i] << endl;
        }
#endif

    // compute the distribution of g2 around g1 for the lower leaflet
    for (int j = 0; j < g1_lower.size(); j++)
        {
        GCoord p1 = g1_lower[j].centerOfMass();
        for (int k = 0; k < g2_lower.size(); k++)
            {
            // skip "self" pairs
            if (g1_lower[j] == g2_lower[k])
                {
                continue;
                }
            GCoord p2 = g2_lower[k].centerOfMass();
            GCoord displ = (p2 - p1);
            displ.reimage(box);
            double d2 = (displ.x() * displ.x()) + (displ.y() * displ.y());
            if ( (d2 <= max2) && (d2 >= min2) )
                {
                double d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist_lower[bin]++;
                }
            }
        }

    // compute the distribution of g2 around g1 for the upper leaflet
    for (int j = 0; j < g1_upper.size(); j++)
        {
        GCoord p1 = g1_upper[j].centerOfMass();
        for (int k = 0; k < g2_upper.size(); k++)
            {
            // skip "self" pairs
            if (g1_upper[j] == g2_upper[k])
                {
                continue;
                }
            GCoord p2 = g2_upper[k].centerOfMass();
            GCoord displ = (p2 - p1);
            displ.reimage(box);
            double d2 = (displ.x() * displ.x()) + (displ.y() * displ.y());
            if ( (d2 <= max2) && (d2 >= min2) )
                {
                double d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist_upper[bin]++;
                }
            }
        }

    frame++;
    }

area /= frame;

int num_upper, num_lower;
// Normalization is slightly different is this is a self distribution
if (group1 == group2)
    {
    // the two selection groups are the same
    num_upper = g1_upper.size()*(g1_upper.size()-1);
    num_lower = g1_lower.size()*(g1_lower.size()-1);
    }
else
    {
    // the two selection groups are different
    num_upper = g1_upper.size() * g2_upper.size();
    num_lower = g1_lower.size() * g2_lower.size();
    }

double upper_expected = frame * num_upper / area;
double lower_expected = frame * num_lower / area;

// Output the results
cout << "# Dist\tTotal\tUpper\tLower" << endl;
for (int i = 0; i < num_bins; i++)
    {
    double d = bin_width*(i + 0.5);

    double d_inner = bin_width*i;
    double d_outer = d_inner + bin_width;
    double norm = M_PI*(d_outer*d_outer - d_inner*d_inner);

    double upper = hist_upper[i]/(norm*upper_expected);
    double lower = hist_lower[i]/(norm*lower_expected);
    double total = (hist_upper[i] + hist_lower[i])/
                        (norm*(upper_expected + lower_expected));

    cout << d << "\t";
    cout << total << "\t";
    cout << upper << "\t";
    cout << lower << endl;

    }
}

