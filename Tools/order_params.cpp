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


using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


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


// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage
{
public:

  void addGeneric(po::options_description& o)
  {
    o.add_options()
      ("1", "Use 1 residue lipids")
      ("3", "Use 3 residue lipids")
      ("y_axis,Y", "Use y axis as magnetic field")
      ("x_axis,X", "Use x axis as magnetic field")
      ("timeseries,T", po::value<string>(&timeseries_filename), "File name for outputing timeseries")
      ("block_average", po::value<string>(&block_filename),"File name for block averaging data")
      ("ba_first", po::value<int>(&ba_first), "Lower range of blocks to average over to calculate uncertainty")
      ("ba_last", po::value<int>(&ba_last), "Upper range of blocks to average over to calculate uncertainty");
  }

  void addHidden(po::options_description& o)
  {
    o.add_options()
      ("sel", po::value<string>(&selection), "Selection string for carbons")
      ("first_carbon", po::value<int>(&first_carbon), "Number of first carbon")
      ("last_carbon", po::value<int>(&last_carbon), "Number of last carbon");
  }

  void addPositional(po::positional_options_description& p) {
    p.add("sel", 1);
    p.add("first_carbon", 1);
    p.add("last_carbon", 1);
  }

  bool check(po::variables_map& vm)
  {
    if (!(vm.count("sel") && vm.count("first_carbon") && vm.count("last_carbon")))
      return(true);
    if (vm.count("1") && vm.count("3"))
      {
        cerr << "Can't select \"--1\" and \"--3\" at the same time"
             << endl
             << "Your lipids either have 1 residue or 3, not both"
             << endl;
        return(true);
      }

    if (vm.count("y_axis") && vm.count("x_axis"))
      {
        cerr << "Can't specify \"--y_axis\" and \"--x_axis\" at the same time"
             << endl
             << "You can only compute the order parameters for 1 magnetic field "
             << endl
             << "at a time."
             << endl;
        return(true);
      }

    return(false);
  }

  bool postConditions(po::variables_map& vm) {
    if (vm.count("1"))
      {
        one_res_lipid = true;
      }
    else if (vm.count("3"))
      {
        three_res_lipid = true;
      }

    if (vm.count("y_axis"))
      {
        axis_index = 1;
      }
    else if (vm.count("x_axis"))
      {
        axis_index = 0;
      }
    else // default to z-axis
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
        if (!vm.count("ba_last")) ba_last = 5;
      }


    return(true);
  }

  string help() const {
    return("selection first_carbon_number last_carbon_number");
  }

  string print() const {
    ostringstream oss;

    oss << boost::format("axis_index=%d, one_res_lipid=%d, three_res_lipid=%d, dump_timeseries=%d, block_average=%d, ba_first=%d, ba_last=%d")
      % axis_index
      % one_res_lipid
      % three_res_lipid
      % dump_timeseries
      % block_average
      % ba_first
      % ba_last;
    return(oss.str());
  }



};

// @endcond

string fullHelpMessage(void)
    {
    string s =
"\n"
"    SYNOPSIS\n"
"\n"
"    Compute lipid order parameters from C-H bond vectors\n"
"\n"
"    DESCRIPTION\n"
"\n"
"    This program computes order parameters analogous to those measured \n"
"    via deuterium quadrupolar splitting on perdeuterated lipid chains.\n"
"    It can correctly handle lipids specified using either 1 residue per lipid\n"
"    (e.g. POPC) or 3 residues per lipid (e.g. PALM-PCGL-OLEO), although it \n"
"    requires the presence of explicit hydrogens on the chains.  In addition, \n"
"    it can calculate the order parameters assuming the magnetic field is \n"
"    oriented along the x- or y-axis, in addition to the more conventional z. \n"
"    Finally, it has the ability to correctly estimate statistical errors in \n"
"    the order parameters using block averaging applied to the full system.\n"
"\n"
"    Formats for lipids\n"
"\n"
"    The CHARMM27 force field typically represented lipids using 3 residues\n"
"    (e.g. POPC was PALM-PCGL-OLEO), while the newer CHARMM36 force field\n"
"    represents each lipid as a single residue.  The code handles both\n"
"    formats; the user has the option to explicitly tell it which format to \n"
"    expect using the flags --1 and --3 (for 1 and 3 residues per lipid, \n"
"    respectively).  If neither flag is invoked, the code attempts to figure \n"
"    it out based on the names of the carbons in the selection; if the \n"
"    selection contains names like C3, its presumed to be 3/lipid format, \n"
"    while the presence of names like C213 suggests 1/lipid format.  The \n"
"    code outputs which case it thinks it has found, so that the user can \n"
"    override it if necessary.\n"
"\n"
"\n"
"    Specifying which carbons to examine\n"
"\n"
"    This is done in 2 parts: First, you give a selection on the command line\n"
"    that says which chains you want to examine.  Second, you say which carbon\n"
"    positions you want to look at within that range (e.g. 2 - 16 for a \n"
"    palmitoyl chain).  The program figures out which atoms match, and \n"
"    identifies the relevant hydrogen atoms.  See the example section for more \n"
"    details.\n"
"\n"
"    Note: the file used to specify the system _must_ have connectivity \n"
"    information (e.g. a psf file) for this program to work correctly.  \n"
"    If you're using gmxdump2pdb.pl to make a PSF/PDB from gromacs, you'll\n"
"    need to run it with the --hydrogens option.\n"
"\n"
"    Orientation of the magnetic field\n"
"\n"
"    As a rule, order parameters are measured with the magnetic field \n"
"    presumed to be oriented along the membrane normal.  However, with \n"
"    oriented samples one has the option to tilt the sample with respect \n"
"    to the field.  The present version of the code supports a subset of \n"
"    this capability, allowing you to compute the order parameters as if\n"
"    the field were oriented along the x or y axis, using the flags --x_axis\n"
"    and --y_axis.\n"
"\n"
"    Note: these results are indistinguishable experimentally, and should \n"
"    be identical in the simulation up to statistical noise.  If you're\n"
"    interested in these results, it's probably best to compute both x- \n"
"    and y-axis fields in order to verify that nothing funny is going on.\n"
"\n"
"    Estimating statistical error using bootstrapping\n"
"\n"
"    --block_average filename        turn on block averaging, output \n"
"                                    additional info to filename (see\n"
"                                    below for explanation of the file\n"
"                                    format)\n"
"    --ba_first  integer             lower number of blocks to consider\n"
"                                    as part of the plateau (col 2 in the\n"
"                                    block_average file).\n"
"    --ba_last  integer              upper number of blocks to consider\n"
"                                    as part of the plateau (col 2 in the\n"
"                                    block_average file).\n"
"    --timeseries  filename          output a timeseries for the average \n"
"                                    order parameter for each carbon position\n"
"\n"
"    Estimating the statistical error in the order parameters computed from\n"
"    molecular dynamics simulations is tricky, but it's not clear what \n"
"    constitutes a single independent sample.  Clearly, the orientation of \n"
"    a given C-H bond is correlated in successive frames.  However, there is\n"
"    also the correlation between the 2 or 3 hydrogens bound to the same \n"
"    carbon.  More problematically, _all_ of the bonds are correlated to the\n"
"    area per lipid for the system as a whole, which fluctuates slowly, \n"
"    often on the timescale of tens of nanoseconds.  So, while there may\n"
"    be a total of 100 lipids in a system, that doesn't mean each frame has \n"
"    anything like 200 (or even 100) independent measurements for the order\n"
"    parameter at a particular carbon position.  \n"
"\n"
"    The approach taken here is use apply block averaging (Flyvbjerg and \n"
"    Petersen, J Chem Phys, 1989, 91, 461-467), applied to the instantenous\n"
"    value of the order parameter for a given carbon position, averaged over \n"
"    the full system.  In this method, the standard error in the mean from\n"
"    a time series is estimated by computing the variance in the average \n"
"    values computed over consecutive subsets of the time series; if the \n"
"    value plateaus as the blocks get long, then that value is the standard\n"
"    error.  However, since the plateau is often noisy (since you don't have\n"
"    as many of the long blocks), the code offers the option to average \n"
"    over a range of block sizes.  \n"
"    \n"
"    Since I can't know in advance what range of block sizes will produce a \n"
"    plateau for any individual series, I let the user control this range.  \n"
"    I suggest running the program twice: run it with a guess for the \n"
"    averaging range (or use the defaults). Then, plot the data in the file \n"
"    containing the block average output and look at what range of block sizes\n"
"    (if any) contains the plateau.  This file is formatted as follows: the \n"
"    first column specifies the carbon number, the second the block number,\n"
"    the third the block size in frames, and the 4th the standard error \n"
"    estimate; plot the 4th vs the 3rd to get the variation in the standard\n"
"    error with block size (plot \"file\" using 3:4 w lp in gnuplot, for \n"
"    example), and then read off from the file which blocks are in that range.\n"
"\n"
"    If you wish to do your own analysis of the variation of the average \n"
"    order parameter, you can specify the --timeseries flag, which will dump\n"
"    out the average value for each carbon position at each time point, with\n"
"    each carbon in its own column.\n"
"    \n"
"    \n"
"    Notes for use with Gromacs:\n"
"\n"
"    If you want to use order_params to analyze data generated with gromacs, \n"
"    there are two caveats you must be aware of.  First, this program assumes\n"
"    you have explicit hydrogens, and will not work if you use a lipid model\n"
"    without them; if you used implicit hydrogens, you'll have to use the \n"
"    builtin gromacs tool instead.  Second, as mentioned above you'll need\n"
"    to use a model file that has explicit bonds in it, INCLUDING BONDS TO\n"
"    HYDROGEN.  The most common usage of gromacs will have this these bonds\n"
"    replaced by constraints, so when you use gmxdump2pdb.pl to make your\n"
"    fake psf file, be sure to give the \"--hydrogens\" option so that the\n"
"    bonds to hydrogen are generated correctly.\n"
"\n"
"    EXAMPLE\n"
"\n"
"    order_params namd.psf merged_1ns.dcd '(segid =~ \"L([0-9]+)\") && resname == \"PALM\"' 2 16 --block_average foo.dat --ba_first 2 --ba_last 5\n"
"\n"
"    This is a moderately complicated example: namd.psf specifies the system,\n"
"    while merged_1ns.dcd is the trajectory.  The next argument is a selection\n"
"    that specifies the palmitate chains for lipids with segment names L and a\n"
"    number; this assumes that the lipids are in the old C27 \n"
"    3-residues-per-lipid format).  The values 2 and 16 specify that we want\n"
"    to look at carbons 2 to 16 within the chain (all of them, in this case).\n"
"    The rest of the command line turns on block averaging (sending the \n"
"    output to foo.dat), and specifies that the plateau region should be\n"
"    estimated as the range between the 2 block data point and the 5 block\n"
"    data point.  \n"
"\n"
"    If you're using 1-residue-per-lipid format (as in the new C36 force \n"
"    field), the selection string would look a bit more complicated.  For\n"
"    example, to select the docosahexaenoyl chain of an SDPE lipid, the string would\n"
"    be \n"
"    'resname == \"SDPE\" && name =~ \"^C2[1-9][0-9]?\"'\n"
"    This selects from residues named \"SDPE\" those atoms with names beginning\n"
"    with C2, followed by a number from 1-9, and perhaps a second integer.\n"
"\n"
"    The resulting output will look something like:\n"
"    # order_params 'namd.psf' 'npgt_bgw.2.dcd' 'resname == \"SDPE\" && name =~ \"^C2[1-9][0-9]?\"' '2' '22' '--block_average' 'foo.dat' - alan (Thu Mar 15 09:23:44 2012) {/home/nleioatts/Public/c36} [2.0.0 120313]\n"
"    # guessing there is 1 residue per lipid\n"
"    # Carbon  S_cd   +/-     BSE\n"
"    2	 0.01811 0.01311 0.00164\n"
"    3	 0.05558 0.02481 0.00787\n"
"    4	 0.08293 0.03523 0.01002\n"
"    5	 0.08138 0.03229 0.00512\n"
"    6	 0.01762 0.02208 0.00303\n"
"    (with more lines)\n"
"\n"
"    The first column is the carbon number, the second is the order parameter,\n"
"    the third is the stdev of the value, and the last column is the \n"
"    estimate of the statistical error from block averaging.  If the command\n"
"    line hadn't included --block_average, the last column would be absent.\n"
        ;
    return(s);
    }


int main (int argc, char *argv[])
{

// parse the command line options
string hdr = invocationHeader(argc, argv);
opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
ToolOptions* topts = new ToolOptions;

opts::AggregateOptions options;
options.add(bopts).add(tropts).add(topts);
if (!options.parse(argc, argv))
  exit(-1);


// *** Copy options back into globals (that aren't already there...) ***
skip = tropts->skip;
system_filename = tropts->model_name;
traj_filename = tropts->traj_name;
vector<uint> framelist = tropts->frameList();

// Create the data structures for the system and trajectory
// (just copied from the BasicTrajectory object)
AtomicGroup system = tropts->model;
pTraj traj = tropts->trajectory;



cout << "# " << invocationHeader(argc, argv) << endl;
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
    AtomicGroup a1 = selectAtoms(main_selection, sel1.c_str());
    AtomicGroup a3 = selectAtoms(main_selection, sel3.c_str());

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
    snprintf(carbon_name, 4, "%d", i);
    string name;
    if (three_res_lipid)
        {
        name = string(" && name == \"C");
        }
    else  // implies one_res_lipid is true
        {
        name = string(" && name =~ \"C[123]");
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
    while ((p = iter()))
        {
        vector<int> atom_ids = p->getBonds();
        AtomicGroup bonded = system.groupFromID(atom_ids);
        AtomicGroup bonded_hydrogens = bonded.select(hyd_sel);
        if (bonded_hydrogens.size() == 0)
            {
            cerr << "No hydrogens found bound to carbon atom: "
                 << p->name()
                 << endl;
       cerr << "This could happen if your lipid doesn't have explicit\n "
            << "hydrogens, if you gave an incorrect selection string,\n "
            << "or if the model file you supplied doesn't have "
            << "connectivity information.";
       cerr << endl;
       cerr << "If you're using a PSF generated by gmxdump2pdb.pl, you\n"
            << "should rerun it using the --hydrogens option so that \n"
            << "bonds to hydrogen are correctly represented.\n";
       cerr << endl;
            exit(-1);
            }

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



const int num_frames = framelist.size();

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

for (uint i=0; i<framelist.size(); ++i)
    {
    traj->readFrame(framelist[i]);
    traj->updateGroupCoords(system);


    // loop over sets of selected carbons
    for (unsigned int i=0; i<selections.size(); i++)
        {
        AtomicGroup *g = &(selections[i]);
        for (uint j=0; j<g->size(); j++)
            {
            // get the carbon
            pAtom carbon = g->getAtom(j);
            // get the relevant hydrogens
            AtomicGroup *hyds = &(hydrogen_list[i][j]);

            AtomicGroup::Iterator iter(*hyds);
            pAtom h;
            while ( (h = iter()) )
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
            timeseries_outfile << boost::format("%8.3f") % fabs(values[i][frame_index]);
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
    double ave = fabs(t.average());
    double dev = t.stdev();

    // get carbon number
    pAtom pa = selections[i].getAtom(0);
    string name = pa->name();
    name.erase(0,1); // delete the C
    if (one_res_lipid)
        {
        name.erase(0,1); // erase the leading 2 or 3
        }
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
