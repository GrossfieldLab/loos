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
string split_by, split_by2;
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
      ("split-mode",po::value<string>(&split_by)->default_value("by-molecule"), "how to split the selections (by-residue, molecule, segment, none)")
      ("split-mode2",po::value<string>(&split_by2)->default_value("by-molecule"), "how to split the second selection (by-residue, molecule, segment, none)")
      ;
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

  void addPositional(po::positional_options_description& p) {
    p.add("sel1", 1);
    p.add("sel2", 1);
    p.add("hist-min", 1);
    p.add("hist-max", 1);
    p.add("num-bins", 1) ;
  }

  bool check(po::variables_map& vm)
  {
    return(!(vm.count("sel1")
             && vm.count ("hist-min")
             && vm.count ("hist-max")
             && vm.count ("num-bins")));
  }

  string help() const
  {
    return("first-selection second-selection histogram-min histogram-max histogram-bins");
  }

  string print() const
  {
    ostringstream oss;
    oss << boost::format("split-mode='%s', sel1='%s', sel2='%s', hist-min=%f, hist-max=%f, num-bins=%f, split-mode2='%'")
      % split_by
      % selection1
      % selection2
      % hist_min
      % hist_max
      % num_bins
      % split_by2;
    return(oss.str());
  }
};
// @endcond


string fullHelpMessage(void)
{
string s =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Compute the radial distribution function for 2 selections, \n"
    "treating the selections as groups as opposed to individual atoms.\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "This tool computes the radial distribution function for 2 selections,\n"
    "treating the selections as groups.  There are 4 ways to group the atoms,\n"
    "controlled by the arguments to --split-mode and --split-mode2: \n"
    "    by-residue: the selection is split into unique residues\n"
    "    by-molecule: the selection is split into unique molecules (only available\n"
    "                if the system file contains connectivity information)\n"
    "    by-segment: the selection is split using the segid (this is present in \n"
    "                CHARMM/NAMD/XPLOR derived files, and some PDB files)\n"
    "    none: treat the entire selection as a single unit\n"
    "\n"
    "The default mode if --split-mode and --split-mode2 aren't set is \"by-molecule\".\n"
    "In all cases, the splitting is performed _before_ the selection is \n"
    "performed, because by-molecule requires the whole system to work correctly.  \n"
    "\n"
    "The distance is then computed between the centers of mass of the grouped \n"
    "objects.\n"
    "If you want to consider individual atoms instead of the centers of mass \n"
    "e.g. if you want to consider all of the individual atoms in a residue),\n"
    "use the tool atomic-rdf instead.\n"
    "\n"
    "histogram-min, histogram-max, and histogram-bins specify the range over \n"
    "which the radial distribution function is computed and the number of bins \n"
    "used.\n"
    "\n"
    "EXAMPLE\n"
    "\n"
    "If the selection string looked like \n"
    "    'resname == \"TRP\" and name =~\"^C\"'\n"
    "with \"by-residue\" splitting, then the full system would first be split \n"
    "into separate residues, and then the selection string would be applied to \n"
    "those individual residues, in this case returning the carbon atoms from \n"
    "the tryptophan residues.  The program would use the center of mass of the\n"
    "carbon atoms to as the point from which to compute the RDF.\n"
    "\n"
    "See also atomic-rdf and xy_rdf.\n"
    ;

    return(s);
    }


enum split_mode { BY_RESIDUE, BY_SEGMENT, BY_MOLECULE, NONE };

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
    else if (!split_by.compare("none"))
        {
        split = NONE;
        }
    else
        {
        cerr << "--split-mode(2) must be: by-residue|by-segment|by-molecule|none"
             << endl;
        exit(-1);
        }
    return (split);
    }

uint doSplit(const AtomicGroup &system, const string selection,
             const split_mode split, vector<AtomicGroup> &grouping)
    {

    // make sure the "result" vector<AG> is empty to start
    grouping.clear();

    vector<AtomicGroup> tmp;
    if (split == BY_MOLECULE)
        {
        tmp = system.splitByMolecule();
        }
    if (split == BY_RESIDUE)
        {
        tmp = system.splitByResidue();
        }
    else if (split == BY_SEGMENT)
        {
        tmp = system.splitByUniqueSegid();
        }
    else if (split == NONE)
        {
        tmp.push_back(system);
        }

    Parser parser(selection);
    KernelSelector parsed_sel(parser.kernel());

    vector<AtomicGroup>::iterator t;
    for (t=tmp.begin(); t!=tmp.end(); ++t)
        {
        AtomicGroup newgroup = t->select(parsed_sel);
        if (newgroup.size() > 0)
            {
            grouping.push_back(newgroup);
            }
        }
    return grouping.size();
    }


int main (int argc, char *argv[])
{

// parse the command line options
string hdr = invocationHeader(argc, argv);
opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
opts::WeightsOptions* wopts = new opts::WeightsOptions;
ToolOptions* topts = new ToolOptions;

opts::AggregateOptions options;
options.add(bopts).add(tropts).add(wopts).add(topts);
if (!options.parse(argc, argv))
  exit(-1);

// parse the split mode, barf if you can't do it
split_mode split=parseSplit(split_by);
split_mode split2=parseSplit(split_by2);

// Print the command line arguments
cout << "# " << hdr << endl;

// Create the system and the trajectory file
// Note: The pTraj type is a Boost shared pointer, so we'll need
//       to use pointer semantics to access it
AtomicGroup system = tropts->model;
pTraj traj = tropts->trajectory;
if (!(system.isPeriodic() || traj->hasPeriodicBox()))
  {
  cerr << "Error- Either the model or the trajectory must have periodic box information.\n";
  exit(-1);
  }

// Attach trajectory to weights
if (wopts->has_weights)
    {
    wopts->weights->add_traj(traj);
    }


double bin_width = (hist_max - hist_min)/num_bins;

// Select the 2 groups, then split them appropriately
vector<AtomicGroup> g1_mols, g2_mols;
uint numgroups = doSplit(system, selection1, split, g1_mols);
if (numgroups == 0)
    {
    cerr << "No groups created by selection1" << endl;
    exit(1);
    }
numgroups = doSplit(system, selection2, split2, g2_mols);
if (numgroups == 0)
    {
    cerr << "No groups created by selection2" << endl;
    exit(1);
    }



// read the initial coordinates into the system
vector<uint> framelist = tropts->frameList();
traj->readFrame(framelist[0]);
traj->updateGroupCoords(system);

// Create the histogram and zero it out
vector<double> hist;
hist.reserve(num_bins);
hist.insert(hist.begin(), num_bins, 0.0);

double min2 = hist_min*hist_min;
double max2 = hist_max*hist_max;

// Precompute the overlap between the two groups (this can be an
// expensive operation, so it's better to have it outside the
// while-loop)

unsigned long unique_pairs = 0;
Math::Matrix<int, Math::RowMajor> group_overlap(g1_mols.size(), g2_mols.size());

for (uint j=0; j<g1_mols.size(); ++j)
    {
    for (uint i=0; i<g2_mols.size(); ++i)
      {
      bool b = (g1_mols[j] == g2_mols[i]);
      group_overlap(j, i) = b;
      if ( !b )
        {
          ++unique_pairs;
        }
      }
    }


// loop over the frames of the trajectory
uint framecount = framelist.size();
double volume = 0.0;
for (uint index = 0; index<framecount; ++index)
    {
    traj->readFrame(framelist[index]);
    // update coordinates and periodic box
    traj->updateGroupCoords(system);
    // if no frame weights file provided, defaults to 1.0 
    const double weight = wopts->weights->get();
    wopts->weights->accumulate();


    GCoord box = system.periodicBox();
    volume += weight*(box.x() * box.y() * box.z());

    // compute the distribution of g2 around g1
    for (unsigned int j = 0; j < g1_mols.size(); j++)
        {
        GCoord p1 = g1_mols[j].centerOfMass();
        for (unsigned int k = 0; k < g2_mols.size(); k++)
            {
            // skip "self" pairs -- in case selection1 and selection2 overlap
            if (group_overlap(j, k))
                {
                continue;
                }
            GCoord p2 = g2_mols[k].centerOfMass();

            // Compute the distance squared, taking periodicity into account
            double d2 = p1.distance2(p2, box);
            if ( (d2 < max2) && (d2 > min2) )
                {
                double d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist[bin]+=weight;
                }
            }
        }
    }
// totalWeight() defaults to frameCount() if no weights file provided
volume /= wopts->weights->totalWeight();
double expected = unique_pairs / volume;
expected *= wopts->weights->totalWeight();

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
    cum1 += hist[i] / (wopts->weights->totalWeight()*g1_mols.size());
    cum2 += hist[i] / (wopts->weights->totalWeight()*g2_mols.size());

    cout << d << "\t" << total << "\t"
         << cum1 << "\t" << cum2 << endl;

    }
}
