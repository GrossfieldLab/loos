/*
 *     Compute lipid order params
 *     Assumes the selections are single-bond carbons, and the following
 *     two atoms are hydrogens
 *
 */
using namespace std;

#include <iostream>
#include "PDB.hpp"
#include <math.h>

void Usage()
    {
    cerr << "Usage: order_params pdb_filelist skip "
         << "selection_file1 [selection_file2...] "
         << endl;
    }

int main (int argc, char *argv[])
{
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 4)
   )
    {
    Usage();
    exit(-1);
    }

// Echo the command line to stdout
cout << "# ";
for (int i =0; i<argc; i++)
    {
    cout << argv[i] << " ";
    }
cout << endl;

// copy the command line variables to real variable names

char *pdb_filelist = argv[1];
int skip = atoi(argv[2]);

vector<Selection> selections;
for (int i =3; i<argc; i++)
    {
    selections.push_back(Selection(argv[i]));
    }

// Read the list of pdb files
vector<string> file_list;
ifstream pdb_files(pdb_filelist);
if (pdb_files.fail())
    {
    cerr << "couldn't read " << pdb_filelist << ": exiting..." << endl;
    exit(-1);
    }
string name;
while ( getline(pdb_files, name) )
    {
    file_list.push_back(name);
    }

// remove the number of frames specified by skip
file_list.erase(file_list.begin(), file_list.begin()+skip);

int num_pdbfiles = file_list.size();

// read the first file
PDBFile file(file_list[0].c_str());

// make the selections
vector<Group *> carbons;
for (vector<Selection>::iterator s =selections.begin();
                                 s!=selections.end();
                                 s++)
    {
    Group *g = new Group;
    file.select(*s, *g);
    carbons.push_back(g);
    }

// set up storage to accumulate the averages
vector<float> sums;
vector<int> counts;
sums.insert(sums.begin(), carbons.size(), 0.0);
counts.insert(counts.begin(), carbons.size(), 0);

// loop over pdb files
for (int i = 0; i < num_pdbfiles; i++)
    {
    // get the new coordinates
    file.update_coor(file_list[i].c_str());

    // loop over carbon groups
    for (int j=0; j<carbons.size(); j++)
        {
        Group *g = carbons[j];
        for (vector<Atom*>::iterator c =g->atoms.begin();
                                     c!=g->atoms.end();
                                     c++)
            {
            // identify the two hydrogens attached
            // a->index = the index of a in file + 1
            // so the index of the 2 hydrogens will by a->index and a->index+1
            int i = (*c)->index;
            Coor *h1 = &(file.atoms[i].coor);
            Coor *h2 = &(file.atoms[i+1].coor);

            Coor v = (*c)->coor - h1;
            float cos_val = v.z / v.length();
            sums[j] += 0.5 - 1.5*cos_val*cos_val;

            v = (*c)->coor - h2;
            cos_val = v.z / v.length();
            sums[j] += 0.5 - 1.5*cos_val*cos_val;

            counts[j] += 2;
            }
        }
    }

for (int i = 0; i < carbons.size(); i++)
    {
    float val = sums[i] / counts[i];
    // get carbon number
    string name = carbons[i]->atoms[0]->name;
    name.erase(0,1); // delete the C
    int index = atoi(name.c_str());

    //cout << i << "\t" 
    //     << val << "\t"
    //     << carbons[i]->atoms[0]->print() << endl;
    cout << index << "\t" << val << endl;

    }

}

