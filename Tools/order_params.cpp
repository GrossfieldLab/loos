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

#if 0
void Usage()
    {
    cerr << "Usage: order_params system traj skip selection "
         << "first_carbon last_carbon [1|3]"
         << endl
         << endl;
    cerr << "The code will attempt to deduce whether you're using "
         << "one or three residues per lipid molecule.  To force it, "
         << "give the value 1 or 3 as the last argument."
         << endl;
    }
#endif

string system_filename;
string traj_filename;
int skip;
string selection;
int first_carbon, last_carbon;
bool one_res_lipid = false;
bool three_res_lipid = false;

void parseOptions(int argc, char *argv[])
    {
    try
        {
        po::options_description generic("Allowed options");
        generic.add_options()
            ("help,h", "Produce this help message")
            ("1", "Use 1 residue lipids")
            ("3", "Use 3 residue lipids");

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

        // Verify sanity
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

        }
    catch(exception& e)
        {
        cerr << "Error - " << e.what() << endl;
        exit(-1);
        }

    }


int main (int argc, char *argv[])
{

#if 0
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 7)
   )
    {
    Usage();
    exit(-1);
    }
#endif

// parse the command line options
parseOptions(argc, argv);

cout << "# " << invocationHeader(argc, argv) << endl;

#if 0
// copy the command line variables to real variable names
char *system_filename = argv[1];
char *traj_filename = argv[2];
int skip = atoi(argv[3]);
char *sel= argv[4];
int first_carbon = atoi(argv[5]);
int last_carbon = atoi(argv[6]);
#endif

#if 0
int one_or_three = 0;
if (argc > 7)
    {
    one_or_three = atoi(argv[7]);
    if ( (one_or_three != 1) && (one_or_three != 3) )
        {
        cerr << "If set, the last argument must be \"1\" or \"3\": "
             << one_or_three
             << endl;
        Usage();
        exit(-1);
        }
    }
#endif 

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


// set up storage to accumulate the averages
// we're going to dump all of the data from a given selection into one
// big lump, so the dimension of sums and counts should match the number of 
// selections specified 
vector<float> sums_x, sums_y, sums_z;
vector<float> sums2_x, sums2_y, sums2_z;
vector<int> counts;
sums_x.insert(sums_x.begin(), selections.size(), 0.0);
sums2_x.insert(sums2_x.begin(), selections.size(), 0.0);

sums_y.insert(sums_y.begin(), selections.size(), 0.0);
sums2_y.insert(sums2_y.begin(), selections.size(), 0.0);

sums_z.insert(sums_z.begin(), selections.size(), 0.0);
sums2_z.insert(sums2_z.begin(), selections.size(), 0.0);

counts.insert(counts.begin(), selections.size(), 0);

// loop over pdb files
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
                double cos_val =  v.z()/length;
                double order = 0.5 - 1.5*cos_val*cos_val;
                sums_z[i] += order;
                sums2_z[i] += order*order;
            
                cos_val =  v.x()/length;
                order = 0.5 - 1.5*cos_val*cos_val;
                sums_x[i] += order;
                sums2_x[i] += order*order;
            
                cos_val =  v.y()/length;
                order = 0.5 - 1.5*cos_val*cos_val;
                sums_y[i] += order;
                sums2_y[i] += order*order;

                counts[i]++;
                }
            }
        }
    }

// Print header
cout << "# Carbon  S_cd(z)   +/-   S_cd(x)   +/-   S_cd(y)   +/-" << endl;

for (unsigned int i = 0; i < selections.size(); i++)
    {
    double ave_x = sums_x[i] / counts[i];
    double ave2_x = sums2_x[i] / counts[i];
    double dev_x = sqrt(ave2_x - ave_x*ave_x);

    double ave_y = sums_y[i] / counts[i];
    double ave2_y = sums2_y[i] / counts[i];
    double dev_y = sqrt(ave2_y - ave_y*ave_y);

    double ave_z = sums_z[i] / counts[i];
    double ave2_z = sums2_z[i] / counts[i];
    double dev_z = sqrt(ave2_z - ave_z*ave_z);

    // get carbon number
    pAtom pa = selections[i].getAtom(0);
    string name = pa->name();
    name.erase(0,1); // delete the C
    int index = atoi(name.c_str());

#if 0
    cout << index 
         << "\t" << ave_z << "\t" << dev_z 
         << "\t" << ave_x << "\t" << dev_x 
         << "\t" << ave_y << "\t" << dev_y 
         << endl;
#endif
    cout << boost::format("%d\t%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f") %
                          index % 
                          ave_z % dev_z %
                          ave_x % dev_x %
                          ave_y % dev_y;
    cout << endl;

    }

}

