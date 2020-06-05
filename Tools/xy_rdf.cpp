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

#include <loos.hpp>


using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

string selection1, selection2;
string split_by;
double hist_min, hist_max;
int num_bins;
int timeseries_interval;
string output_directory;
bool sel1_spans, sel2_spans;
bool reselect_leaflet = false;


// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage
{
public:

  void addGeneric(po::options_description& o)
  {
    o.add_options()
      ("split-mode",po::value<string>(&split_by)->default_value("by-molecule"), "how to split the selections (by-residue, molecule, segment)")
      ("timeseries", po::value<int>(&timeseries_interval)->default_value(0), "Interval to write out timeseries, 0 means never")
      ("timeseries-directory", po::value<string>(&output_directory)->default_value(string("output")))
      ("sel1-spans", "Selection 1 appears in both leaflets")
      ("sel2-spans", "Selection 2 appears in both leaflets")
      ("reselect", "Recompute leaflet location for each frame")
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
    p.add("num-bins", 1);
  }

  bool check(po::variables_map& vm)
  {
    return(!(vm.count("sel1")
             && vm.count("hist-min")
             && vm.count("hist-max")
             && vm.count("num-bins")));
  }

  bool postConditions(po::variables_map& vm)
  {
    if (vm.count("sel1") && !vm.count("sel2"))
      selection2 = selection1;
    sel1_spans = false;
    sel2_spans = false;
    if (vm.count("sel1-spans"))
        {
        sel1_spans = true;
        }
    if (vm.count("sel2-spans"))
        {
        sel2_spans = true;
        }

    if (vm.count("reselect"))
        {
        reselect_leaflet = true;
        }
    if ( (timeseries_interval > 0) && vm.count("weights"))
        {
        cerr << "Cannot specify reweighting and time series at the same time"
             << endl;
        return(false);
        }
    return(true);
  }

  string help() const
  {
    return("first-selection second-selection histogram-min histogram-max histogram-bins");
  }

  string print() const
  {
    ostringstream oss;
    oss << boost::format("split-mode='%s', sel1='%s', sel2='%s', hist-min=%f, hist-max=%f, num-bins=%f, timeseries=%d, timeseries-directory='%s', sel1-spans=%d, sel2-spans=%d reselect=%d")
      % split_by
      % selection1
      % selection2
      % hist_min
      % hist_max
      % num_bins
      % timeseries_interval
      % output_directory
      % sel1_spans
      % sel2_spans
      % reselect_leaflet;
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

string fullHelpMessage(void)
{
string s =
    "\n"
    "SYNOPSIS\n"
    "\n"
    "Compute a radial distribution function in the x-y plane\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "This tool is intended primarily for analyzing lateral structure of\n"
    "membrane systems. As with rdf, it operates primarily on groups of atoms.\n"
    "There are 3 ways to group the atoms, controlled by the\n"
    "arguments to --split-mode: \n"
    "\n"
    "  by-residue:  the selection is split into unique residues\n"
    "  by-molecule: the selection is split into unique molecules (only \n"
    "               available if the system file contains connectivity \n"
    "               information)\n"
    "  by-segment:  the selection is split using the segid (this is present \n"
    "               in CHARMM/NAMD/XPLOR derived files, and some PDB files)\n"
    "\n"
    "The default mode if --split-mode isn't set is \"by-molecule\".\n"
    "In each case, the splitting is performed _before_ the selection is \n"
    "performed.  \n"
    "\n"
    "The distance is then computed between the centers of mass of the grouped\n"
    "objects, only considering the x and y coordinates.  The program treats\n"
    "the two leaflets of the separately, based on the sign of the z-coordinate\n"
    "of the center of mass of the selection in the first frame; this can \n"
    "cause problems if the membrane has not already been centered at the \n"
    "origin (the merge-traj tool can do this for you).\n"
    "\n"
    "If one or both of the selections should be included in both leaflets\n"
    "(e.g. a transmembrane helix or protein), the user should specify \n"
    "the --sel1-spans or --sel2-spans flags (applying to the first and \n"
    "second selections respectively).  Otherwise, the selection will \n"
    "be included in only one leaflet, which could lead to non-sensical \n"
    "results.  Note: this flag is applied to _all_ of the components \n"
    "in the selection --- they are each assumed to span the membrane.\n"
    "\n"
    "The --timeseries flag lets you track the evolution of the RDF over time.\n"
    "\n"
    "The --reselect flag causes the program to recompute which leaflet\n"
    "each molecule is in at each time step. This will impose a small \n"
    "overhead, but is necessary if you're dealing with molecules that \n"
    "can flip from one leaflet to the other.\n"
    "\n"
    "EXAMPLE\n"
    "\n"
    "To look at the distribution of PE lipid headgroups in a lipid\n"
    "bilayer, you might use a command line like\n"
    "\n"
    "xy_rdf model-file traj-file 'resname == \"PEGL\"' 'resname == \"PEGL\"' 0 40 40\n"
    "    --split-mode=by-molecule\n"
    "\n"
    "Assuming the CHARMM27-style lipid naming, the headgroup would be its own\n"
    "residue with name \"PEGL\", and the result would be the lateral RDF for \n"
    "the headgroup centers of mass.  \n"
    "\n"
    "As with the other rdf tools (rdf, atomic-rdf), histogram-min,\n"
    "histogram-max, and histogram-bins specify the range over which the\n"
    "radial distribution function is computed and the number of bins used.  \n"
    "\n"
    "The --timeseries flag lets you track the evolution of the rdf over time, by\n"
    "writing out a windowed average as it is accumulated.  So, the adding the \n"
    "flags\n"
    "\n"
    "      --timeseries 100 --timeseries-directory \"foo\"\n"
    "\n"
    "to the end of the command above would cause the program to write out a \n"
    "new average every 100 frames considering only the frames in that interval.\n"
    "The files will appear in the directory \"foo\", with names rdf_0.dat, \n"
    "rdf_1.dat, etc.  The program does not attempt to create \"foo\" if it \n"
    "doesn't exist, and instead will simply exit.\n"
    "\n"
    "Note: the 5th column (\"Cum\") is not a density like the other values, \n"
    "but rather the absolute number of molecules of the second selection \n"
    "found around the first selection.\n"
    ;

    return (s);
    }

void assign_leaflet(vector<AtomicGroup> &molecules,
                    vector<AtomicGroup> &upper,
                    vector<AtomicGroup> &lower,
                    const bool spans
                    )
    {
    upper.clear();
    lower.clear();
    if (spans)
        {
        upper = molecules;
        lower = molecules;
        return;
        }

    for (unsigned int i = 0; i < molecules.size(); i++)
        {
        GCoord c = molecules[i].centerOfMass();
        if (c.z() >=0.0)
            {
            upper.push_back(molecules[i]);
            }
        else
            {
            lower.push_back(molecules[i]);
            }
        }
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
options.add(bopts).add(tropts).add(topts).add(wopts);
if (!options.parse(argc, argv))
  exit(-1);



// parse the split mode, barf if you can't do it
split_mode split=parseSplit(split_by);

cout << "# " << invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names
AtomicGroup system = tropts->model;
pTraj traj = tropts->trajectory;;
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

AtomicGroup group1 = selectAtoms(system, selection1);
if (group1.empty())
    {
    cerr << "Error- no atoms selected by '" << selection1 << "'\n";
    exit(-1);
    }
group1.pruneBonds();

AtomicGroup group2 = selectAtoms(system, selection2);
if (group2.empty())
    {
    cerr << "Error- no atoms selected by '" << selection2 << "'\n";
    exit(-1);
    }
group2.pruneBonds();

// Split the groups into chunks, depending on how the user asked
// us to.
vector<AtomicGroup> g1_mols, g2_mols;
if (split == BY_MOLECULE)
    {
    g1_mols = group1.splitByMolecule();
    g2_mols = group2.splitByMolecule();
    }
else if (split == BY_RESIDUE)
    {
    g1_mols = group1.splitByResidue();
    g2_mols = group2.splitByResidue();
    }
else if (split == BY_SEGMENT)
    {
    g1_mols = group1.splitByUniqueSegid();
    g2_mols = group2.splitByUniqueSegid();
    }

// read the initial coordinates into the system
traj->updateGroupCoords(system);

// Now that we have some real coordinates, we need to subdivide the groups
// one more time, into upper and lower leaflets. This assumes that the
// coordinates are properly centered and imaged.
vector<AtomicGroup> g1_upper, g1_lower;
vector<AtomicGroup> g2_upper, g2_lower;

assign_leaflet(g1_mols, g1_upper, g1_lower, sel1_spans);
assign_leaflet(g2_mols, g2_upper, g2_lower, sel2_spans);


// Create 2 histograms -- one for top, one for bottom
// Also create 2 histograms to store the total
vector<double> hist_lower, hist_upper;
vector<double> hist_lower_total, hist_upper_total;
hist_lower.reserve(num_bins);
hist_upper.reserve(num_bins);
hist_lower_total.reserve(num_bins);
hist_upper_total.reserve(num_bins);
hist_lower.insert(hist_lower.begin(), num_bins, 0.0);
hist_upper.insert(hist_upper.begin(), num_bins, 0.0);
hist_lower_total.insert(hist_lower_total.begin(), num_bins, 0.0);
hist_upper_total.insert(hist_upper_total.begin(), num_bins, 0.0);


double min2 = hist_min*hist_min;
double max2 = hist_max*hist_max;

// loop over the frames of the traj file
double area = 0.0;
double interval_area = 0.0;
double cum_upper_pairs = 0.0;
double cum_lower_pairs = 0.0;
double interval_upper_pairs = 0.0;
double interval_lower_pairs = 0.0;

vector<uint> framelist = tropts->frameList();
uint framecnt = framelist.size();


for (uint index = 0; index<framecnt; ++index)
    {
    // update coordinates and periodic box
    traj->readFrame(framelist[index]);
    traj->updateGroupCoords(system);

    // If this is a reweighting calculation, get the weight.
    // Otherwise, use 1.0
    double weight = 1.0;
    if (wopts->has_weights)
        {
        weight = wopts->weights->get();
        wopts->weights->accumulate();
        }

    GCoord box = system.periodicBox();
    area += weight*(box.x() * box.y());
    interval_area += weight*(box.x() * box.y());

    if (reselect_leaflet)
        {
        assign_leaflet(g1_mols, g1_upper, g1_lower, sel1_spans);
        assign_leaflet(g2_mols, g2_upper, g2_lower, sel2_spans);
        }

    // compute the distribution of g2 around g1 for the lower leaflet
    for (unsigned int j = 0; j < g1_lower.size(); j++)
        {
        GCoord p1 = g1_lower[j].centerOfMass();
        for (unsigned int k = 0; k < g2_lower.size(); k++)
            {
            // skip "self" pairs
            if (g1_lower[j] == g2_lower[k])
                {
                continue;
                }

            cum_lower_pairs+=weight;
            interval_lower_pairs+=weight;

            GCoord p2 = g2_lower[k].centerOfMass();
            GCoord displ = (p2 - p1);
            displ.reimage(box);
            double d2 = (displ.x() * displ.x()) + (displ.y() * displ.y());
            if ( (d2 < max2) && (d2 > min2) )
                {
                double d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist_lower[bin]+=weight;
                }
            }
        }

    // compute the distribution of g2 around g1 for the upper leaflet
    for (unsigned int j = 0; j < g1_upper.size(); j++)
        {
        GCoord p1 = g1_upper[j].centerOfMass();
        for (unsigned int k = 0; k < g2_upper.size(); k++)
            {
            // skip "self" pairs
            if (g1_upper[j] == g2_upper[k])
                {
                continue;
                }

            cum_upper_pairs+=weight;
            interval_lower_pairs+=weight;

            GCoord p2 = g2_upper[k].centerOfMass();
            GCoord displ = (p2 - p1);
            displ.reimage(box);
            double d2 = (displ.x() * displ.x()) + (displ.y() * displ.y());
            if ( (d2 < max2) && (d2 > min2) )
                {
                double d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist_upper[bin]+=weight;
                }
            }
        }


    // if requested, write out timeseries as well
    if (timeseries_interval && (index % timeseries_interval == 0))
        {
        interval_area /= timeseries_interval;
        double upper_expected = interval_upper_pairs / interval_area;
        double lower_expected = interval_lower_pairs / interval_area;


        // create the output file
        ostringstream outfilename;
        outfilename << output_directory << "/" << "rdf_" << index << ".dat";
        ofstream out(outfilename.str().c_str());
        if (out.fail())
            {
            cerr << "couldn't open " << outfilename.str() << " ... exiting" << endl;
            exit(-1);
            }
        out << "# Dist\tTotal\tUpper\tLower\tCum" << endl;
        double cum = 0.0;
        for (int m = 0; m < num_bins; m++)
            {
            double d = bin_width*(m + 0.5);

            double d_inner = bin_width*m;
            double d_outer = d_inner + bin_width;
            double norm = M_PI*(d_outer*d_outer - d_inner*d_inner);

            double upper = 0.0;
            if (interval_upper_pairs > 0)
            {
              upper = hist_upper[m]/(norm*upper_expected);
            }

            double lower = 0.0;
            if (interval_lower_pairs > 0)
            {
              lower = hist_lower[m]/(norm*lower_expected);
            }

            double total = (hist_upper[m] + hist_lower[m])/
                                (norm*(upper_expected + lower_expected) );
            cum += (hist_upper[m] + hist_lower[m])/(group1.size()*timeseries_interval);

            out << d << "\t"
                << total << "\t"
                << upper << "\t"
                << lower << "\t"
                << cum   << endl;

            }

        out << endl; // blank line for gnuplot
        out.close();

        // sum up the total histograms
        for (int i=0; i<num_bins; ++i)
            {
            hist_upper_total[i] += hist_upper[i];
            hist_lower_total[i] += hist_lower[i];
            }

        // rezero the histograms
        hist_upper.assign(num_bins, 0.0);
        hist_lower.assign(num_bins, 0.0);

        // zero out the area
        interval_area = 0.0;
        interval_upper_pairs = 0;
        interval_lower_pairs = 0;
        }

    }

// normalize the area
// Not necessary if we're doing reweighting, since they're already
// normalized
if (!wopts->has_weights) area /= framecnt;

// If we didn't write timeseries, then we need to copy the interval histograms
// to the total ones.  If we did, we need to add in the additional data points
// since the last time we wrote out a times series file
if (!timeseries_interval)
    {
    hist_lower_total = hist_lower;
    hist_upper_total = hist_upper;
    }
else if (framecnt % timeseries_interval != 0)
    {
    for (int i=0; i< num_bins; i++)
        {
        hist_lower_total[i] += hist_lower[i];
        hist_upper_total[i] += hist_upper[i];
        }
    }

double upper_expected = cum_upper_pairs / area;
double lower_expected = cum_lower_pairs / area;

// Output the results
cout << "# Dist\tTotal\tUpper\tLower\tCum" << endl;
double cum = 0.0;
for (int i = 0; i < num_bins; i++)
    {
    double d = bin_width*(i + 0.5);

    double d_inner = bin_width*i;
    double d_outer = d_inner + bin_width;
    double norm = M_PI*(d_outer*d_outer - d_inner*d_inner);

    double upper = 0.0;
    if (cum_upper_pairs > 0)
    {
      upper = hist_upper_total[i]/(norm*upper_expected);
    }

    double lower = 0.0;
    if (cum_lower_pairs > 0)
    {
      lower = hist_lower_total[i]/(norm*lower_expected);
    }

    double total = (hist_upper_total[i] + hist_lower_total[i])/
                        (norm*(upper_expected + lower_expected));
    if (wopts->has_weights)
        {
        upper /= wopts->weights->totalWeight();
        lower /= wopts->weights->totalWeight();
        total /= wopts->weights->totalWeight();
        }

    double cum_increment = (hist_upper_total[i] + hist_lower_total[i]) /
                                    group1.size();
    if (wopts->has_weights)
        {
        cum_increment *= framecnt / wopts->weights->totalWeight();
        }

    cum += cum_increment;

    cout << d     << "\t"
         << total << "\t"
         << upper << "\t"
         << lower << "\t"
         << cum   << endl;

    }
}
