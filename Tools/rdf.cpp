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

using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

string selection1, selection2;
string split_by;
double hist_min, hist_max;
int num_bins;
int skip;

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage
{
public:
  
  void addGeneric(po::options_description& o) 
  {
    o.add_options()
      ("split-mode",po::value<string>(&split_by), "how to split the selections");
  }

  void addHidden(po::options_description& o)
  {
    o.add_options()
      ("sel1", po::value<string>(&selection1), "first selection")
      ("sel2", po::value<string>(&selection2), "second selection")
      ("hist-min", po::value<double>(&hist_min), "Histogram minimum")
      ("hist-max", po::value<double>(&hist_max), "Histogram maximum")
      ("num-bins", po::value<int>(&num_bins), "Histogram bins"); 
  }

  bool check(po::variables_map& vm)
  {
    return(!(vm.count("sel1") 
             && vm.count("hist-min")
             && vm.count("hist-max")
             && vm.count("num-bins")));
  }

  bool postCondition(po::variables_map& vm)
  {
    if (vm.count("sel1") && !vm.count("sel2"))
      selection2 = selection1;
    return(true);
  }

  string help() const
  {
    return("first-selection second-selection histogrm-min histogram-min histogram-bins");
  }

  string print() const
  {
    ostringstream oss;
    oss << boost::format("split-mode='%s', sel1='%s', sel2='%s', hist-min=%f, hist-max=%f, num-bins=%f")
      % split_by
      % selection1
      % selection2
      % hist_min
      % hist_max
      % num_bins;
    return(oss.str());
  }

};
// @endcond


enum split_mode { BY_RESIDUE, BY_SEGMENT, BY_MOLECULE };

split_mode parseSplit(const string &split_by)
    {
    split_mode split;
    // figure out the splitting mode
    if (!split_by.compare("by-residue"))
        {
        split = BY_RESIDUE;
        }
    else if (!split_by.compare("by-segment"))
        {
        split = BY_SEGMENT;
        }
    else if (!split_by.compare("by-molecule"))
        {
        split = BY_MOLECULE;
        }
    else
        {
        cerr << "--split-mode must be: by-residue|by-segment|by-molecule"
             << endl;
        cerr << "If unset, defaults to by-molecule";
        cerr << endl;
        exit(-1);
        }
    return (split);
    }


int main (int argc, char *argv[])
{

// parse the command line options
string hdr = invocationHeader(argc, argv);
opts::BasicOptions* bopts = new opts::BasicOptions;
opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
ToolOptions* topts = new ToolOptions;

opts::AggregateOptions options;
options.add(bopts).add(tropts).add(topts);
if (!options.parse(argc, argv))
  exit(-1);

// parse the split mode, barf if you can't do it
split_mode split=parseSplit(split_by);

// Print the command line arguments
cout << "# " << hdr << endl;

// Create the system and the trajectory file
// Note: The pTraj type is a Boost shared pointer, so we'll need
//       to use pointer semantics to access it
AtomicGroup system = tropts->model;
pTraj traj = tropts->trajectory;

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
vector<AtomicGroup> grouping;
if (split == BY_MOLECULE)
    {
    grouping = system.splitByMolecule();
    }
else if (split == BY_RESIDUE)
    {
    grouping = system.splitByResidue();
    }
else if (split == BY_SEGMENT)
    {
    grouping = system.splitByUniqueSegid();
    }

// Set up the selector to define group1 atoms
Parser parser1(selection1);
KernelSelector parsed_sel1(parser1.kernel());

// Set up the selector to define group2 atoms
Parser parser2(selection2);
KernelSelector parsed_sel2(parser2.kernel());

// Loop over the molecules and add them to group1 or group2
vector<AtomicGroup> g1_mols, g2_mols;
vector<AtomicGroup>::iterator g;
for (g=grouping.begin(); g!=grouping.end(); g++)
    {
    AtomicGroup tmp = g->select(parsed_sel1);
    if (tmp.size() > 0)
        {
        g1_mols.push_back(tmp);
        }

    AtomicGroup tmp2 = g->select(parsed_sel2);
    if (tmp2.size() > 0)
        {
        g2_mols.push_back(tmp2);
        }
    }


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

