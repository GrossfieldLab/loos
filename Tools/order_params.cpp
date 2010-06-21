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


void Usage()
    {
    cerr << "Usage: order_params system traj skip selection "
         << "first_carbon last_carbon [1|3]"
         << endl;
    cerr << "The code will attempt to deduce whether you're using "
         << "one or three residues per lipid molecule.  To force it, "
         << "give the value 1 or 3 as the last argument."
         << endl;
    }

int main (int argc, char *argv[])
{
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 7)
   )
    {
    Usage();
    exit(-1);
    }

cout << "# " << invocationHeader(argc, argv) << endl;

// copy the command line variables to real variable names

char *system_filename = argv[1];
char *traj_filename = argv[2];
int skip = atoi(argv[3]);
char *sel= argv[4];
int first_carbon = atoi(argv[5]);
int last_carbon = atoi(argv[6]);

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

// Create the data structures for the system and trajectory
AtomicGroup system = createSystem(system_filename);
pTraj traj = createTrajectory(traj_filename, system);

// TODO: add a check for connectivity

// NOTE: We assume the selection is a list of all of the relevant carbon atoms. 
//       We'll split it into the individual carbon positions ourselves, then
//       figure out the relevant hydrogens from connectivity.
AtomicGroup main_selection = selectAtoms(system, sel);

// Do we need to figure out how many residues per lipid?
if ( (one_or_three != 1) && (one_or_three != 3) )
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
        one_or_three = 3;
        cout << "# guessing there are 3 residues per lipid"
             << endl;
        }
    else if ( (a3.size() > 0) && (a1.size() == 0) )
        {
        // found only 3 digit numbers and no 1 digit ones, so
        // it must be 1 residue per lipid
        one_or_three = 1;
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
        one_or_three = 3;
        }
    }

// Now break into individual carbons
vector<AtomicGroup> selections;
for (int i =first_carbon; i<=last_carbon; i++)
    {
    string sel_string = string(sel);
    char carbon_name[4];
    sprintf(carbon_name, "%d", i);
    string name;
    if (one_or_three == 3)
        {
        name = string(" && name == \"C");
        }
    else if (one_or_three == 1)
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
vector<float> sums;
vector<float> sums2;
vector<int> counts;
sums.insert(sums.begin(), selections.size(), 0.0);
sums2.insert(sums2.begin(), selections.size(), 0.0);
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
                double cos_val =  v.z() / v.length();
                double order = 0.5 - 1.5*cos_val*cos_val;
                sums[i] += order;
                sums2[i] += order*order;
                counts[i]++;
                }
            }
        }
    }

for (unsigned int i = 0; i < selections.size(); i++)
    {
    double ave = sums[i] / counts[i];
    double ave2 = sums2[i] / counts[i];
    double dev = sqrt(ave2 - ave*ave);

    // get carbon number
    pAtom pa = selections[i].getAtom(0);
    string name = pa->name();
    name.erase(0,1); // delete the C
    int index = atoi(name.c_str());

    cout << index << "\t" << ave << "\t" << dev << endl;

    }

}

