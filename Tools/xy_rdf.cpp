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

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage
{
public:
  
  void addGeneric(po::options_description& o) 
  {
    o.add_options()
      ("split-mode",po::value<string>(&split_by)->default_value("by-molecule"), "how to split the selections (by-residue, molecule, segment)")
      ("timeseries", po::value<int>(&timeseries_interval)->default_value(0), "Interval to write out timeseries, 0 means never")
      ("timeseries-directory", po::value<string>(&output_directory)->default_value(string("output")));

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

  bool postCondition(po::variables_map& vm)
  {
    if (vm.count("sel1") && !vm.count("sel2"))
      selection2 = selection1;
    return(true);
  }

  string help() const
  {
    return("first-selection second-selection histogram-min histogram-min histogram-bins");
  }

  string print() const
  {
    ostringstream oss;
    oss << boost::format("split-mode='%s', sel1='%s', sel2='%s', hist-min=%f, hist-max=%f, num-bins=%f, timeseries=%d, timeseries-directory='%s'")
      % split_by
      % selection1
      % selection2
      % hist_min
      % hist_max
      % num_bins
      % timeseries_interval
      % output_directory;
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

cout << "# " << invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names
AtomicGroup system = tropts->model;
pTraj traj = tropts->trajectory;

double bin_width = (hist_max - hist_min)/num_bins;

AtomicGroup group1 = selectAtoms(system, selection1);
group1.pruneBonds();

AtomicGroup group2 = selectAtoms(system, selection2);
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

for (unsigned int i = 0; i < g1_mols.size(); i++)
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

for (unsigned int i = 0; i < g2_mols.size(); i++)
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


// Figure out the number of unique pairs -- if the groups are the same,
// we skip self pairs, so the normalization is different
uint num_upper, num_lower;
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


double min2 = hist_min*hist_min;
double max2 = hist_max*hist_max;

// loop over the frames of the traj file
int frame = 0;
double area = 0.0;
double interval_area = 0.0;


while (traj->readFrame())
    {
    // update coordinates and periodic box
    traj->updateGroupCoords(system);
    GCoord box = system.periodicBox(); 
    area += box.x() * box.y();
    interval_area += box.x() * box.y();

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
            GCoord p2 = g2_lower[k].centerOfMass();
            GCoord displ = (p2 - p1);
            displ.reimage(box);
            double d2 = (displ.x() * displ.x()) + (displ.y() * displ.y());
            if ( (d2 < max2) && (d2 > min2) )
                {
                double d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist_lower[bin]++;
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
            GCoord p2 = g2_upper[k].centerOfMass();
            GCoord displ = (p2 - p1);
            displ.reimage(box);
            double d2 = (displ.x() * displ.x()) + (displ.y() * displ.y());
            if ( (d2 < max2) && (d2 > min2) )
                {
                double d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist_upper[bin]++;
                }
            }
        }

    frame++;

    // if requested, write out timeseries as well
    if (timeseries_interval && (frame % timeseries_interval == 0))
        {
        interval_area /= timeseries_interval;
        double upper_expected = timeseries_interval * num_upper / interval_area;
        double lower_expected = timeseries_interval * num_lower / interval_area;


        // create the output file
        ostringstream outfilename;
        outfilename << output_directory << "/" << "rdf_" << frame << ".dat";
        ofstream out(outfilename.str().c_str());
        if (out.fail())
            {
            cerr << "couldn't open " << outfilename << " ... exiting" << endl;
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

            double upper = hist_upper[m]/(norm*upper_expected);
            double lower = hist_lower[m]/(norm*lower_expected);
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
        }

    }

// normalize the area
area /=frame;

// If we didn't write timeseries, then we need to copy the interval histograms
// to the total ones.  If we did, we need to add in the additional data points
// since the last time we wrote out a times series file
if (!timeseries_interval)
    {
    hist_lower_total = hist_lower;
    hist_upper_total = hist_upper;
    }
else if (frame % timeseries_interval != 0)
    {
    for (int i=0; i< num_bins; i++)
        {
        hist_lower_total[i] += hist_lower[i];
        hist_upper_total[i] += hist_upper[i];
        }
    }

double upper_expected = frame * num_upper / area;
double lower_expected = frame * num_lower / area;

// Output the results
cout << "# Dist\tTotal\tUpper\tLower\tCum" << endl;
double cum = 0.0;
for (int i = 0; i < num_bins; i++)
    {
    double d = bin_width*(i + 0.5);

    double d_inner = bin_width*i;
    double d_outer = d_inner + bin_width;
    double norm = M_PI*(d_outer*d_outer - d_inner*d_inner);

    double upper = hist_upper_total[i]/(norm*upper_expected);
    double lower = hist_lower_total[i]/(norm*lower_expected);
    double total = (hist_upper_total[i] + hist_lower_total[i])/
                        (norm*(upper_expected + lower_expected));
    cum += (hist_upper_total[i] + hist_lower_total[i])
                    /(group1.size()*frame);

    cout << d     << "\t"
         << total << "\t"
         << upper << "\t"
         << lower << "\t"
         << cum   << endl;

    }
}

