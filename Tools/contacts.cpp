/*
 *     Compute 2d radial distribution function for two selections
 */
using namespace std;

#include <iostream>
#include "PDB.hpp"
#include <math.h>

void Usage()
    {
    cerr << "Usage: contacts pdb_filelist selection_file1 selection_file2 "
         << "max" 
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

// Echo the command line to stdout
cout << "# ";
for (int i =0; i<argc; i++)
    {
    cout << argv[i] << " ";
    }
cout << endl;

// copy the command line variables to real variable names
char *pdb_filelist = argv[1];
char *selection_file1 = argv[2];
char *selection_file2 = argv[3];
float max= atof(argv[4]);
float max2 = max*max;


Selection s1 = Selection(selection_file1);
Selection s2 = Selection(selection_file2);

// Read the list of pdb files
vector<string> file_list;
ifstream pdb_files(pdb_filelist);
if (pdb_files.bad())
    {
    cerr << "couldn't read " << pdb_filelist << ": exiting..." << endl;
    exit(-1);
    }
string name;
while ( getline(pdb_files, name) )
    {
    file_list.push_back(name);
    }

int num_pdbfiles = file_list.size();


// read the first file
PDBFile file(file_list[0].c_str());

// make the selections
Group group1, group2;
file.select(s1, group1);
file.select(s2, group2);

if (group1.num_atoms == 0)
    {
    cerr << "No atoms in group1" << endl;
    exit(0);
    }

if (group2.num_atoms == 0)
    {
    cerr << "No atoms in group2" << endl;
    exit(0);
    }

cout << "#Frame\tPairs\tPerGroup1\tPerGroup2" << endl;

// loop over the pdb files
for (int i = 0; i < num_pdbfiles; i++)
    {
    // get the new coordinates
    file.update_coor(file_list[i].c_str());
    int count = 0;

    // compute the number of contacts between group1 atoms and group2 atoms
    for (int j = 0; j < group1.num_atoms; j++)
        {
        Atom *a1 = group1.atoms[j];
        for (int k = 0; k < group2.num_atoms; k++)
            {
            Atom *a2 = group2.atoms[k];
            if (a1 == a2)
                {
                continue;
                }
            float d2 = a1->dist_squ(*a2, file.box);
            if (d2 <= max2)
                {
                count++;
                }
            }
        }
    
    // Output the results
    float per_g1_atom = (float)count / group1.num_atoms;
    float per_g2_atom = (float)count / group2.num_atoms;
    cout << i << "\t" 
         << count << "\t"
         << per_g1_atom << "\t"
         << per_g2_atom << endl;
    }
}

