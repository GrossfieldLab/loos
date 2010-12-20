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
#include <boost/program_options.hpp>

using namespace std;
using namespace loos;

namespace po = boost::program_options;

#if 0
void Usage()
    {
    cerr << "Usage: xy_rdf system traj selection1 selection2 "
         << "min max num_bins skip" 
         << endl;
    }
#endif

string system_filename;
string traj_filename;
string selection1, selection2;
string split_by;
double hist_min, hist_max;
int num_bins;
int skip;
int timeseries_interval;
string output_directory;

void parseOptions(int argc, char *argv[])
    {
    try 
        {
        po::options_description generic("Allowed options");
        generic.add_options()
            ("help,h", "Produce this help message")
            //("fullhelp", "Even more help")
            ("sel1", po::value<string>(&selection1), "first selection")
            ("sel2", po::value<string>(&selection2), "second selection")
            ("split-mode",po::value<string>(&split_by), "how to split the selections")
            ("skip", po::value<int>(&skip)->default_value(0), "frames to skip")
            ("timeseries", po::value<int>(&timeseries_interval)->default_value(0), "Interval to write out timeseries, 0 means never")
            ("timeseries-directory", po::value<string>(&output_directory)->default_value(string("foo")));
        
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("model", po::value<string>(&system_filename), "Model filename")
            ("traj", po::value<string>(&traj_filename), "Trajectory filename")
            ("hist-min", po::value<double>(&hist_min), "Histogram minimum")
            ("hist-max", po::value<double>(&hist_max), "Histogram maximum")
            ("num-bins", po::value<int>(&num_bins), "Histogram bins"); 

        po::options_description command_line;
        command_line.add(generic).add(hidden);

        po::positional_options_description p;
        p.add("model", 1);
        p.add("traj", 1);
        p.add("hist-min", 1);
        p.add("hist-max", 1);
        p.add("num-bins", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help") || 
            !vm.count("model") || !vm.count("traj") || 
            !vm.count("hist-min") || !vm.count("hist-max") ||
            !vm.count("num-bins") ||
            !vm.count("sel1") 
           )
            {
            cerr << "Usage: " << argv[0] << " "
                 << "model-name trajectory-name hist-min hist-max num-bins "
                 << "--split-mode by-residue|by-segment|by-molecule"
                 << "--sel1 SELECTION [--sel2 SELECTION] "
                 << "--timeseries VALUE --timeseries-directory DIR"
                 << endl;
            cerr << generic;
            exit(-1);
            }
        
        // if there's only 1 selection, duplicate it
        if (!vm.count("sel2"))
            {
            selection2 = selection1;
            }
        }
    catch(exception& e) 
        {
        cerr << "Error - " << e.what() << endl;
        exit(-1);
        }
    }

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

#if 0
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 9)
   )
    {
    Usage();
    exit(-1);
    }
#endif

// parse the command line options
parseOptions(argc, argv);

// parse the split mode, barf if you can't do it
split_mode split=parseSplit(split_by);

cout << "# " << invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names
AtomicGroup system = createSystem(system_filename);
pTraj traj = createTrajectory(traj_filename, system);

double bin_width = (hist_max - hist_min)/num_bins;

AtomicGroup group1 = selectAtoms(system, selection1);
AtomicGroup group2 = selectAtoms(system, selection2);

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

// Skip the initial frames
traj->readFrame(skip); 

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
int num_upper, num_lower;
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
double total_area = 0.0;


while (traj->readFrame())
    {
    // update coordinates and periodic box
    traj->updateGroupCoords(system);
    GCoord box = system.periodicBox(); 
    area += box.x() * box.y();

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
        total_area += area;
        area /= timeseries_interval;
        double upper_expected = timeseries_interval * num_upper / total_area;
        double lower_expected = timeseries_interval * num_lower / total_area;


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
            double norm = M_PI*(d_outer*d_outer - d_inner*d_inner)/area;

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
        area = 0;
        
        }

    }

// normalize the area
if (timeseries_interval)
    {
    total_area /= frame;
    }
else
    {
    total_area = area/frame;
    }

// if we didn't write timeseries, then we need to copy the interval histograms
// to the total ones
if (timeseries_interval)
    {
    hist_lower_total = hist_lower;
    hist_upper_total = hist_upper;
    }

double upper_expected = frame * num_upper / total_area;
double lower_expected = frame * num_lower / total_area;

// Output the results
cout << "# Dist\tTotal\tUpper\tLower\tCum" << endl;
double cum = 0.0;
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
    cum += (hist_upper[i] + hist_lower[i])/(group1.size()*timeseries_interval);

    cout << d     << "\t"
         << total << "\t"
         << upper << "\t"
         << lower << "\t"
         << cum   << endl;

    }
}

