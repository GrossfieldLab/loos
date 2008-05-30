/*
  charge_density.cpp
  (c) 2008 Tod D. Romo and Alan Grossfield

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Compute the charge density along the z dimension of a system

*/

#include <ios>
#include <iostream>
#include <iomanip>

#include <loos.hpp>
#include <psf.hpp>
#include <dcd.hpp>
#include <utils.hpp>
#include <Parser.hpp>
#include <Selectors.hpp>



int main(int argc, char *argv[]) {

    cout << "# " << invocationHeader(argc, argv) << endl;

    PSF psf(argv[1]);
    DCD dcd(argv[2]);
    int num_skip = atoi(argv[3]);
    double min_z = atof(argv[4]);
    double max_z = atof(argv[5]);
    int num_bins = atoi(argv[6]);

    // TODO: add an arbitrary number of selections, and track the charge 
    // density from each selection
    vector<AtomicGroup> subsets;
    subsets.push_back(psf);
    for (int i=7; i<argc; i++) {
        Parser parser(argv[i]);
        KernelSelector parsed_sel(parser.kernel());
        AtomicGroup g = psf.select(parsed_sel);
        subsets.push_back(g);
    }
  
    double bin_width = (max_z - min_z) / num_bins;
  
    // Create and zero out the total charge_distribution
    vector< vector<double> > charge_dists;
    charge_dists.reserve(subsets.size());
    for (unsigned int i=0; i<=subsets.size(); i++) {
        vector<double> *v = new vector<double>;
        v->insert(v->begin(), num_bins, 0.0);
        charge_dists.push_back(*v);
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
            while ( pa = iter() ) {
                double c = pa->charge();
                double z = pa->coords().z();

                if ( (z > min_z) && (z < max_z) ) {
                    int bin = (int)( (z - min_z) / bin_width );
                    charge_dists[i][bin] += c;
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
            double c = charge_dists[j][i] / frame;
            cout << c << "\t";
        }
        cout << endl;
    }


  
}
