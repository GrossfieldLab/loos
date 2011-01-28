/*
 *     Compute lipid order params
 *     Assumes the selections are single-bond carbons, and the following
 *     two atoms are hydrogens
 *
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
#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;

string system_filename;
string timeseries_filename;
int skip;
string selection;
int first_carbon, last_carbon;
int axis_index;
bool one_res_lipid = false;
bool three_res_lipid = false;

bool dump_timeseries = false;
string traj_filename;

bool block_average = false;
string block_filename;
int ba_first, ba_last;

void parseOptions(int argc, char *argv[])
    {
    try
        {
        po::options_description generic("Allowed options");
        generic.add_options()
            ("help,h", "Produce this help message")
            ("1", "Use 1 residue lipids")
            ("3", "Use 3 residue lipids")
            ("y_axis,y", "Use y axis as magnetic field")
            ("x_axis,x", "Use x axis as magnetic field")
            ("timeseries,t", po::value<string>(&timeseries_filename), "File name for outputing timeseries")
            ("block_average,b", po::value<string>(&block_filename),"File name for block averaging data")
            ("ba_first,f", po::value<int>(&ba_first), "Lower range of blocks to average over to calculate uncertainty")
            ("ba_last,l", po::value<int>(&ba_last), "Upper range of blocks to average over to calculate uncertainty");

        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("model", po::value<string>(&system_filename), "Model filename")
            ("traj", po::value<string>(&traj_filename), "Trajectory filename")
            ("sel", po::value<string>(&selection), "Selection string for carbons")
            ("skip", po::value<int>(&skip), "Frames to skip")
            ("first_carbon", po::value<int>(&first_carbon), "Number of first carbon")
            ("last_carbon", po::value<int>(&last_carbon), "Number of last carbon");
            
        po::options_description command_line;
        command_line.add(generic).add(hidden);

        po::positional_options_description p;
        p.add("model", 1);
        p.add("traj", 1);
        p.add("skip", 1);
        p.add("sel", 1);
        p.add("first_carbon", 1);
        p.add("last_carbon", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc,argv).
                options(command_line).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")              ||
                !vm.count("model")        || !vm.count("traj")         ||
                !vm.count("skip")         || !vm.count("sel")          ||
                !vm.count("first_carbon") || !vm.count("last_carbon")
           )
            {
            cerr << "Usage: " << argv[0] << " "
                 << "model-name trajectory-name skip-frames selection-string "
                 << "first-carbon last-carbon "
                 << "[--1|--3]"
                 << endl;
            cerr << generic;
            exit(-1);
            }

        // Verify sanity of the options for number of residues per lipid
        if (vm.count("1") && vm.count("3"))
            {
            cerr << "Can't select \"--1\" and \"--3\" at the same time" 
                 << endl
                 << "Your lipids either have 1 residue or 3, not both" 
                 << endl;
            exit(-1);
            }
        else if (vm.count("1"))
            {
            one_res_lipid = true;
            }
        else if (vm.count("3"))
            {
            three_res_lipid = true;
            }

        // Verify sanity of options for which axis to use as magnetic field direction
        if (vm.count("y_axis") && vm.count("x_axis"))
            {
            cerr << "Can't specify \"--y_axis\" and \"--x_axis\" at the same time"
                 << endl
                 << "You can only compute the order parameters for 1 magnetic field "
                 << endl
                 << "at a time.  Exiting..."
                 << endl;
            exit(-1);
            }
        if (vm.count("y_axis"))
            {
            axis_index = 1;
            }
        else if (vm.count("x_axis"))
            {
            axis_index = 0;
            }
        else  // default to using the z-axis
            {
            axis_index = 2;
            }

        // Did the user request time series output?
        if (vm.count("timeseries"))
            {
            dump_timeseries = true;
            }
        else
            {
            dump_timeseries = false;
            }

        // Did the user request block averaging?
        // If so, set the boolean flag
        if (vm.count("block_average"))
            {
            block_average = true;
            if (!vm.count("ba_first")) ba_first = 2;
            if (!vm.count("ba_last")) ba_first = 5;
            }

        }
    catch(exception& e)
        {
        cerr << "Error - " << e.what() << endl;
        exit(-1);
        }

    }


int main (int argc, char *argv[])
{

// parse the command line options
parseOptions(argc, argv);

cout << "# " << invocationHeader(argc, argv) << endl;

// Create the data structures for the system and trajectory
AtomicGroup system = createSystem(system_filename);
pTraj traj = createTrajectory(traj_filename, system);

// TODO: add a check for connectivity

// NOTE: We assume the selection is a list of all of the relevant carbon atoms. 
//       We'll split it into the individual carbon positions ourselves, then
//       figure out the relevant hydrogens from connectivity.
AtomicGroup main_selection = selectAtoms(system, selection);

// Do we need to figure out how many residues per lipid?
// If both are false, then neither option was set on the command line
if ( !one_res_lipid && !three_res_lipid)
    {
    // Number of residues per lipid wasn't set, try to figure it out.
    // If we find things that look like "C3" (one number after the carbon),
    // it has to be 3 residues per lipid.  If we find things that look like
    // "C213" (3 digits after the C), it has to be 1.  If we find neither, 
    // it could be because the user just selected a subset of carbon positions
    // to look at; we'll guess 3-residues per and inform the user.  
    // If we find both, we're screwed, and will die screaming.
    // Note: this gets CHARMM (pre- and post-C36 parameters) right, but hasn't been 
    // tested on anything else.  GROMACS doesn't have hydrogens, so this whole code
    // won't work, and I don't have access to anything else (eg AMBER GAFF).
    string sel1 = string("name =~ \"^C[1-9]$\"");
    string sel3 = string("name =~ \"^C[1-9][0-9][0-9]$\"");
    AtomicGroup a1, a3;
    try 
        {
        a1 = selectAtoms(main_selection, sel1.c_str());
        }
    catch (loos::NullResult e)
        {
        // discard this error, it just means nothing matched
        }
    try 
        {
        a3 = selectAtoms(main_selection, sel3.c_str());
        }
    catch (loos::NullResult e)
        {
        // discard this error, it just means nothing matched
        }

    if ( (a1.size() > 0) && (a3.size() == 0) )
        {
        // found only 1 digit numbers and no 3 digit numbers,
        // so it must be 3 residues per lipid
        three_res_lipid = true;
        cout << "# guessing there are 3 residues per lipid"
             << endl;
        }
    else if ( (a3.size() > 0) && (a1.size() == 0) )
        {
        // found only 3 digit numbers and no 1 digit ones, so
        // it must be 1 residue per lipid
        one_res_lipid = true;
        cout << "# guessing there is 1 residue per lipid"
             << endl;
        }
    else if ( (a1.size() > 0) && (a3.size() > 0) )
        {
        // found both 1 and 3 digit numbers, which is very weird.
        // Notify the user and tell them to pick manually.
        cerr << "Couldn't figure out whether you have 1 or 3 residues "
             << "per lipid molecules."
             << endl;
        cerr << "You'll need to specify this manually.  Exiting...."
             << endl;
        exit(-1);
        }
    else if ( (a1.size() == 0) && (a3.size() == 0) )
        {
        // we didn't find 1 or 3 digit numbers.  This can happen if the
        // user only looks at a subset of carbons -- eg carbons 5-7 in one 
        // residue format, or carbons 12-15 in 3 residue format.
        cerr << "Can't unambiguously tell whether you've got 1 or 3 "
             << "residues per lipid.  I'm guessing 3, "
             << endl;
        cerr << "but if this guess is wrong, you'll need to rerun with "
             << "the correct value specified."
             << endl;
        three_res_lipid = true;
        }
    }

// At this point, we must have either three_res_lipid or one_res_lipid set.
// This means we only have to test one to decide what to do.
assert(one_res_lipid != three_res_lipid);

// Now break into individual carbons
vector<AtomicGroup> selections;
for (int i =first_carbon; i<=last_carbon; i++)
    {
    char carbon_name[4];
    string sel_string = selection;
    sprintf(carbon_name, "%d", i);
    string name;
    if (three_res_lipid)
        {
        name = string(" && name == \"C");
        }
    else  // implies one_res_lipid is true
        {
        name = string(" && name =~ \"C[23]");
        }
    name += string(carbon_name) + string("\"");
    sel_string.insert(sel_string.size(), name);
    selections.push_back(selectAtoms(main_selection, sel_string.c_str()));
    }

// Now, figure out which hydrogens go with each carbon selected
// hydrogen_list mimics the structure of selections, so the hydrogens
// bound to the jth carbon of the ith selection will be found at
// hydrogen_list[i][j].
vector< vector<AtomicGroup> > hydrogen_list(selections.size());
HydrogenSelector hyd_sel;
for (unsigned int i=0; i<selections.size(); i++) 
    {
    AtomicGroup *s = &(selections[i]);
    AtomicGroup::Iterator iter(*s);
    pAtom p;
    while (p = iter()) 
        {
        vector<int> atom_ids = p->getBonds();
        AtomicGroup bonded = system.groupFromID(atom_ids);
        AtomicGroup bonded_hydrogens = bonded.select(hyd_sel);
        hydrogen_list[i].push_back(bonded_hydrogens);
        }
    }

#ifdef DEBUG
// Check to see if the correct hydrogens were found
for (unsigned int i=0; i<selections.size(); i++)
    {
    AtomicGroup *g = &(selections[i]);
    cerr << "total atoms in sel " << i << "= " << g->size() << endl;
    for (int j=0; j<g->size(); j++)
        {
        pAtom carbon = g->getAtom(j);
        // get the relevant hydrogens
        AtomicGroup *hyds = &(hydrogen_list[i][j]);
        cerr << *carbon << endl;
        cerr << *hyds << endl;
        }
    }
#endif



// skip the equilibration frames
traj->readFrame(skip);

const int num_frames = traj->nframes() - skip;

// We're going to accumulate the time series of average values for each
// carbon position, and turn this into an average at the end.
// This will let us do better uncertainty analysis.
vector<vector<float> > values;
vector<int> counts;
values.resize(selections.size());
for (uint i=0; i<values.size(); i++)
    {
    values[i].insert(values[i].begin(), num_frames, 0.0);
    }
counts.insert(counts.begin(), selections.size(), 0);

ofstream timeseries_outfile;
if (dump_timeseries)
    {
    timeseries_outfile.open(timeseries_filename.c_str());
    if (!timeseries_outfile.good())
        {
        cerr << "Error opening time series output file " 
             << timeseries_filename
             << endl
             << "Proceeding without outputing time series"
             << endl;
        dump_timeseries = false;
        }
    else
        {
        // write the header
        timeseries_outfile <<"# Timestep\t";
        for (int i=first_carbon; i<=last_carbon; i++)
            {
            timeseries_outfile << i << "\t";
            }
        timeseries_outfile << endl;
        }
        
    }

// loop over frames in the trajectory
int frame_index = 0;
while (traj->readFrame())
    {
    traj->updateGroupCoords(system);

    
    // loop over sets of selected carbons
    for (unsigned int i=0; i<selections.size(); i++)
        {
        AtomicGroup *g = &(selections[i]);
        for (int j=0; j<g->size(); j++)
            {
            // get the carbon
            pAtom carbon = g->getAtom(j);
            // get the relevant hydrogens
            AtomicGroup *hyds = &(hydrogen_list[i][j]);
            
            AtomicGroup::Iterator iter(*hyds);
            pAtom h;
            while (h = iter() )
                {
                GCoord v = carbon->coords() - h->coords();
                double length = v.length();
                double cos_val =  v[axis_index]/length;
                double order = 0.5 - 1.5*cos_val*cos_val;
                values[i][frame_index] += order;
                counts[i]++;
                }
            }
        }
    // turn the sum into an average
    if (dump_timeseries)
        {
        timeseries_outfile << frame_index << "\t";
        }
    for (unsigned int i=0; i<selections.size(); i++)
        {
        values[i][frame_index] /= counts[i];
        if (dump_timeseries)
            {
            timeseries_outfile << boost::format("%8.3f") % values[i][frame_index];
            }
        counts[i] = 0;
        }
    if (dump_timeseries)
        {
        timeseries_outfile << endl;
        }
    frame_index++;
    }

// Print header
if (!block_average)
    {
    cout << "# Carbon  S_cd   +/-" << endl;
    }
else
    {
    cout << "# Carbon  S_cd   +/-     BSE" << endl;
    }

ofstream ba_outfile;
int ba_maxblocks = -1;
if (block_average)
    {
    ba_outfile.open(block_filename.c_str());
    if (!ba_outfile.good())
        {
        cerr << "Failed opening block averaging output file " 
             << block_filename 
             << endl
             << "Turning off block averaging"
             << endl;
        block_average = false;
        }
    else
        {
        // Write a header
        ba_outfile << "# Carb\tBlock\tBlockSize\tStdErr" << endl;

        // Figure out the maximum number of blocks to try
        ba_maxblocks = (int)(frame_index / 10.0);
        }
    }


for (unsigned int i = 0; i < selections.size(); i++)
    {
    TimeSeries<float> t(values[i]);
    double ave = t.average();
    double dev = t.stdev();
    
    // get carbon number
    pAtom pa = selections[i].getAtom(0);
    string name = pa->name();
    name.erase(0,1); // delete the C
    int index = atoi(name.c_str());

    float sum = 0.0;
    cout << boost::format("%d\t%8.5f%8.5f") %
                      index % 
                      ave % dev;
    if (block_average)
        {
        for (int j=2; j<ba_maxblocks; j++)
            {
            float variance = t.block_var(j);
            float std_err = sqrt(variance/j);

            // The "plateau" region of many block averaging plots
            // is noisy, so we average over a range of block sizes
            if ( (j >= ba_first) && (j <= ba_last) )
                {
                sum += std_err;
                }
            float block_size = frame_index / j;
            ba_outfile << boost::format("%d\t%d\t%8.3f\t%8.5f") %
                                        index % j % block_size % std_err;
            ba_outfile << endl;
            }
        ba_outfile << endl;
        float bse = sum / (ba_last - ba_first + 1);
        cout << boost::format("%8.5f") % bse;
        }
    cout << endl;

    }


}

