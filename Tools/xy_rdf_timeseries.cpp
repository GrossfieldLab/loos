/*
 *     Compute 2d radial distribution function for two selections
 */
using namespace std;

#include <iostream>
#include "PDB.hpp"
#include <math.h>

void Usage()
    {
    cerr << "Usage: xy_rdf pdb_filelist selection_file1 selection_file2 "
         << "min max num_bins skip interval out_dir" 
         << endl;
    }

int main (int argc, char *argv[])
{
if ( (argc <= 1) || 
     ( (argc >= 2) && (strncmp(argv[1], "-h", 2) == 0) ) ||
     (argc < 10)
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
float hist_min = atof(argv[4]);
float hist_max = atof(argv[5]);
int num_bins = atoi(argv[6]);
int skip = atoi(argv[7]);
int interval = atoi(argv[8]);
char *dir_name = argv[9];

float bin_width = (hist_max - hist_min)/num_bins;


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

// remove the number of frames specified by skip
file_list.erase(file_list.begin(), file_list.begin()+skip);

int num_pdbfiles = file_list.size();


// read the first file
PDBFile file(file_list[0].c_str());

// make the selections
Group group1, group2;
file.select(s1, group1);
file.select(s2, group2);


// split the selections into upper and lower leaflets
// Assumes proper centering and imaging
// Also assumes that the relevant pieces of the groups are individual atoms
Group g1_upper, g1_lower;
Group g2_upper, g2_lower;

for (int i = 0; i < group1.num_atoms; i++)
    {
    Atom *a = group1.atoms[i];
    if (a->coor.z >=0.0)
        {
        g1_upper.add_atom(a);
        }
    else
        {
        g1_lower.add_atom(a);
        }
    }

for (int i = 0; i < group2.num_atoms; i++)
    {
    Atom *a = group2.atoms[i];
    if (a->coor.z >= 0.0)
        {
        g2_upper.add_atom(a);
        }
    else
        {
        g2_lower.add_atom(a);
        }
    }


// Create 2 histograms -- one for top, one for bottom
vector<float> hist_lower, hist_upper;
hist_lower.reserve(num_bins);
hist_upper.reserve(num_bins);
hist_lower.insert(hist_lower.begin(), num_bins, 0.0);
hist_upper.insert(hist_upper.begin(), num_bins, 0.0);

float min2 = hist_min*hist_min;
float max2 = hist_max*hist_max;

int num_upper, num_lower;
if (string(selection_file1).compare(string(selection_file2)) == 0)
    {
    // the two selection groups are the same
    num_upper = g1_upper.num_atoms*(g1_upper.num_atoms-1)/2;
    num_lower = g1_lower.num_atoms*(g1_lower.num_atoms-1)/2;
    }
else
    {
    // the two selection groups are different
    num_upper = g1_upper.num_atoms * g2_upper.num_atoms/2;
    num_lower = g1_lower.num_atoms * g2_lower.num_atoms/2;
    }

float area = file.box.x * file.box.y;  // ASSUMES FIXED AREA
float upper_expected = interval * num_upper / area;
float lower_expected = interval * num_lower / area;

#if 0
cout << "Upper: " << upper_expected << endl;
cout << "Lower: " << lower_expected << endl;
cout << "g1: " << g1_upper.num_atoms << " " << g1_lower.num_atoms << endl;
cout << "g2: " << g2_upper.num_atoms << " " << g2_lower.num_atoms << endl;
#endif


// loop over the pdb files
for (int i = 0; i < num_pdbfiles; i++)
    {
    // get the new coordinates
    file.update_coor(file_list[i].c_str());

    // compute the distribution of g2 around g1 for the lower leaflet
    for (int j = 0; j < g1_lower.num_atoms; j++)
        {
        Atom *a1 = g1_lower.atoms[j];
        for (int k = 0; k < g2_lower.num_atoms; k++)
            {
            Atom *a2 = g2_lower.atoms[k];
            if (a1 == a2)
                {
                continue;
                }
            Coor displ = a1->displacement(*a2, file.box);
            float d2 = (displ.x * displ.x) + (displ.y * displ.y);
            if ( (d2 <= max2) && (d2 >= min2) )
                {
                float d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist_lower[bin]++;
                //cerr << d << "\t" << bin << "\t" << hist_lower[bin] << endl;
                }
            }
        }


    // compute the distribution of g2 around g1 for the upper leaflet
    for (int j = 0; j < g1_upper.num_atoms; j++)
        {
        Atom *a1 = g1_upper.atoms[j];
        for (int k = 0; k < g2_upper.num_atoms; k++)
            {
            Atom *a2 = g2_upper.atoms[k];
            if (a1 == a2)
                {
                continue;
                }
            Coor displ = a1->displacement(*a2, file.box);
            float d2 = (displ.x * displ.x) + (displ.y * displ.y);
            if ( (d2 <= max2) && (d2 >= min2) )
                {
                float d = sqrt(d2);
                int bin = int((d-hist_min)/bin_width);
                hist_upper[bin]++;
                }
            }
        }

    // periodically dump out rdf
    if (i % interval == (interval-1) )
        {
        // create the output file
        ostringstream outfilename; 
        outfilename << dir_name << "/" << "rdf_" << i+1 << ".dat";
        ofstream out(outfilename.str().c_str());
        if (out.fail())
            {
            cerr << "couldn't open " << outfilename << " ... exiting" << endl;
            exit(-1);
            }
        out << "# Dist\tTotal\tUpper\tLower" << endl;
        float cum = 0.0;
        for (int m = 0; m < num_bins; m++)
            {
            float d = bin_width*(m + 0.5);

            float d_inner = bin_width*m;
            float d_outer = d_inner + bin_width;
            float norm = M_PI*(d_outer*d_outer - d_inner*d_inner);

            float upper = hist_upper[m]/(norm*upper_expected);
            float lower = hist_lower[m]/(norm*lower_expected);
            // The 2.0 is because there are 2 leaflets, so twice as much area
            float total = (hist_upper[m] + hist_lower[m])/
                                (2.0*norm*(upper_expected + lower_expected) );
            cum += (hist_upper[m] + hist_lower[m])/(group1.num_atoms*interval);

            out << d << "\t"
                << total << "\t"
                << upper << "\t"
                << lower << "\t"
                << cum   << endl;

            }

        out << endl; // blank line for gnuplot
        out.close();
        // rezero the histograms
        hist_upper.assign(num_bins, 0.0);
        hist_lower.assign(num_bins, 0.0);
        }
    }
}

