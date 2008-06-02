/*
 *     Compute lipid order params
 *     Assumes the selections are single-bond carbons, and the following
 *     two atoms are hydrogens
 *
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
    cerr << "Usage: order_params psf dcd skip "
         << "selection1 [selection2 ...]"
         << endl;
    }

int main (int argc, char *argv[])
{
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 5)
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

// Create the data structures for the system and dcd
PSF psf(psf_filename);
DCD dcd(dcd_filename);

// NOTE: We assume each selection is a list of carbon atoms, and we'll figure
//       out the identities of the bound hydrogens ourselves
vector<AtomicGroup> selections;
for (int i =4; i<argc; i++)
    {
    Parser p(argv[i]);
    KernelSelector parsed(p.kernel());
    selections.push_back(psf.select(parsed));
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

