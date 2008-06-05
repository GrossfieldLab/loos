/*
  density-dist-windowed.cpp
  (c) 2008 Tod D. Romo and Alan Grossfield

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Compute the charge/mass/electron density along the z dimension of a system

*/

#include <ios>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctype.h>

#include <loos.hpp>
#include <psf.hpp>
#include <dcd.hpp>
#include <utils.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>

void Usage()
    {
    cerr << "Usage: density-dist-windowed "
         << " PSF DCD E|C|M num_frames_to_skip min_z max_z num_bins window_size"
         << " filename_prototype [extra_selection_1 ...]"
         << endl;
    }

int main(int argc, char *argv[]) {

    if ( (argc <= 1) || 
         (strncmp(argv[1], "-h", 2) == 0) ||
         (argc < 10) ) {
        Usage();
        exit(-1);
    }

    cout << "# " << invocationHeader(argc, argv) << endl;


    PSF psf(argv[1]);
    DCD dcd(argv[2]);
    char calc_type = toupper(*argv[3]); // bad programmer, no cookie for you
    int num_skip = atoi(argv[4]);
    double min_z = atof(argv[5]);
    double max_z = atof(argv[6]);
    int num_bins = atoi(argv[7]);
    int window = atoi(argv[8]);
    string filename_proto = string(argv[9]);

    bool do_charge = false;
    bool do_mass = false;
    bool do_elec = false;

    // Figure out if we're doing charge, electron density, or mass
    if (calc_type == 'C')
        do_charge = true;
    else if (calc_type == 'E')
        do_elec = true;
    else if (calc_type == 'M')
        do_mass = true;
    else
        throw(runtime_error("calc type must be C, E, or M"));


    // Add an arbitrary number of selections, and track the 
    // density from each selection
    vector<AtomicGroup> subsets;
    subsets.push_back(psf);
    for (int i=10; i<argc; i++) {
        Parser parser(argv[i]);
        KernelSelector parsed_sel(parser.kernel());
        AtomicGroup g = psf.select(parsed_sel);
        subsets.push_back(g);
    }
  
    double bin_width = (max_z - min_z) / num_bins;
  
    // Create and zero out the total charge_distribution
    vector< vector<double> > dists;
    dists.reserve(subsets.size());
    for (unsigned int i=0; i<=subsets.size(); i++) {
        vector<double> *v = new vector<double>;
        v->insert(v->begin(), num_bins, 0.0);
        dists.push_back(*v);
    }
    // skip the equilibration frames
    dcd.readFrame(num_skip);
  
    // loop over the remaining frames
    int frame = 0;
    while (dcd.readFrame()) {
        // update coordinates
        dcd.updateGroupCoords(psf);

        // Loop over the subsets and compute charge distribution.
        // (first set is all atoms)
        for (unsigned int i=0; i<subsets.size(); i++) {
            pAtom pa;
            AtomicGroup::Iterator iter(subsets[i]);
            double weight;
            while ( pa = iter() ) {
                if (do_charge)
                    weight = pa->charge();
                else if (do_mass)
                    weight = pa->mass();
                else if (do_elec)
                    weight = pa->atomic_number() - pa->charge();
                double z = pa->coords().z();

                if ( (z > min_z) && (z < max_z) ) {
                    int bin = (int)( (z - min_z) / bin_width );
                    dists[i][bin] += weight;
                }
            }
        }
        frame++;
        // output the results for this averaging window
        if (frame % window == 0) {
            int output_val = frame/window;

            // set up the output file
            char piece[128];
            snprintf(piece, 128, "_%d.dat", output_val);
            string filename = filename_proto + string(piece);
            ofstream outfile(filename.c_str());
            if (!outfile.is_open() ) {
                throw(runtime_error("couldn't open output file"));
                }
            outfile << "# Z\tAllAtoms";
            for (unsigned int i=1; i<subsets.size(); i++) {
                outfile << " Set(" << i <<") "; 
                }
            outfile << endl;

            // Normalize by the number of frames and 
            // output the average density
            for (int i=0; i<num_bins; i++) {
                double z = (i+0.5)*bin_width + min_z;
                outfile << z << "\t"; 
                for (unsigned int j=0; j<subsets.size(); j++) {
                    double c = dists[j][i] / window;
                    outfile << c << "\t";
                }
                outfile << endl;
            }

            // zero out the distributions
            vector< vector<double> >::iterator d;
            for (d=dists.begin(); d!=dists.end(); d++) {
                d->assign(d->size(), 0.0);
            }
            
        }

    }
  
}
