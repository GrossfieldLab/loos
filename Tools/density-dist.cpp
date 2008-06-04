/*
  charge_density.cpp
  (c) 2008 Tod D. Romo and Alan Grossfield

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Compute the charge/mass/electron density along the z dimension of a system

*/

#include <ios>
#include <iostream>
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
    cerr << "Usage: density-dist "
         << " PSF DCD E|C|M num_frames_to_skip min_z max_z num_bins"
         << " [extra_selection_1 ...]"
         << endl;
    }

int main(int argc, char *argv[]) {

    if ( (argc <= 1) || 
         (strncmp(argv[1], "-h", 2) == 0) ||
         (argc < 8) ) {
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


    // TODO: add an arbitrary number of selections, and track the charge 
    // density from each selection
    vector<AtomicGroup> subsets;
    subsets.push_back(psf);
    for (int i=8; i<argc; i++) {
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
    }

    // Normalize by the number of frames and output the average charge density
    cout << "# Z\tAllAtoms";
    for (unsigned int i=1; i<subsets.size(); i++) {
        cout << " Set(" << i <<") "; 
    }
    cout << endl;

    for (int i=0; i<num_bins; i++) {
        double z = (i+0.5)*bin_width + min_z;
        cout << z << "\t"; 
        for (unsigned int j=0; j<subsets.size(); j++) {
            double c = dists[j][i] / frame;
            cout << c << "\t";
        }
        cout << endl;
    }


  
}
