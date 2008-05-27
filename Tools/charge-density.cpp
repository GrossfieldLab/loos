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
#include <Atom.hpp>
#include <psf.hpp>
#include <dcd.hpp>


int main(int argc, char *argv[]) {

    // First, read in the PSF...
    PSF psf(argv[1]);
    DCD dcd(argv[2]);
    int num_skip = atoi(argv[3]);
    double min_z = atof(argv[4]);
    double max_z = atof(argv[5]);
    int num_bins = atoi(argv[6]);
  
    double bin_width = (max_z - min_z) / num_bins;
  
    // Create and zero out the total charge_distribution
    vector<double> charge_dist;
    charge_dist.insert(charge_dist.begin(), num_bins, 0.0);
  
    // skip the equilibration frames
    dcd.readFrame(num_skip);
  
    // loop over the remaining frames
    int frame = 0;
    while (dcd.readFrame()) {
        dcd.updateGroupCoords(psf);
        AtomicGroup::Iterator iter(psf);
        pAtom pa;
        while ( pa = iter() ) {
            //cout << *pa << endl;
            double c = pa->charge();
            double z = pa->coords().z();

            if ( (z > min_z) && (z < max_z) ) {
                int bin = (int)( (z - min_z) / bin_width );
                charge_dist[bin] += c;
#if 0
                cout << z << "  "
                     << bin << "  "
                     << charge_dist[bin]
                     << endl;
#endif
            }
        }
        frame++;
        //if (frame % 100 == 0) cerr << "finished frame " << frame << endl;
    }

    // Normalize by the number of frames and output the average charge density
    cout << "# Z\tCharge(elec)" << endl;
    for (int i=0; i<num_bins; i++) {
        charge_dist[i] /= frame;
        double z = (i+0.5)*bin_width + min_z;
        cout << z << "\t" << charge_dist[i] << endl;
    }


  
}
