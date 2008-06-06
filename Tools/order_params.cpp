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



using namespace std;

#include <iostream>
#include <math.h>

#include <loos.hpp>
#include <psf.hpp>
#include <dcd.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>

void Usage()
    {
    cerr << "Usage: order_params psf dcd skip selection "
         << "first_carbon last_carbon"
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

char *psf_filename = argv[1];
char *dcd_filename = argv[2];
int skip = atoi(argv[3]);
char *sel= argv[4];
int first_carbon = atoi(argv[5]);
int last_carbon = atoi(argv[6]);

// Create the data structures for the system and dcd
PSF psf(psf_filename);
DCD dcd(dcd_filename);

// NOTE: We assume the selection is a list of all of the relevant carbon atoms. 
//       We'll break it into invidual carbons ourselves (assuming the normal
//       convention of C2, C3, ... 
//       Afterward, we'll figure out the relevant hydrogens ourselves
Parser main_p(sel);
KernelSelector main_parsed(main_p.kernel());
AtomicGroup main_selection = psf.select(main_parsed);

// Now break into individual carbons
vector<AtomicGroup> selections;
for (int i =first_carbon; i<=last_carbon; i++)
    {
    string sel_string = string(sel);
    char carbon_name[4];
    sprintf(carbon_name, "%d", i);
    string name = string(" && name == \"C") + string(carbon_name)
                  + string("\"");
    sel_string.insert(sel_string.size(), name);
    Parser p(sel_string.c_str());
    KernelSelector parsed(p.kernel());
    selections.push_back(main_selection.select(parsed));
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
        AtomicGroup bonded = psf.groupFromID(atom_ids);
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
dcd.readFrame(skip);


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
while (dcd.readFrame())
    {
    dcd.updateGroupCoords(psf);

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

