/*
  density-dist.cpp

  Compute the charge/mass/electron density along the z dimension of a system

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


#include <ctype.h>

#include <loos.hpp>

using namespace std;

void Usage()
    {
    cerr << "Usage: density-dist "
         << " system traj E|C|M num_frames_to_skip min_z max_z num_bins"
         << " [extra_selection_1 ...]"
         << endl;
    cerr << "Note: the system file must specify the mass and charge" << endl;
    }

int main(int argc, char *argv[]) {

    if ( (argc <= 1) || 
         (strncmp(argv[1], "-h", 2) == 0) ||
         (argc < 8) ) {
        Usage();
        exit(-1);
    }

    cout << "# " << loos::invocationHeader(argc, argv) << endl;


    AtomicGroup system = loos::createSystem(argv[1]);
    pTraj traj = loos::createTrajectory(argv[2], system);


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
    subsets.push_back(system);
    for (int i=8; i<argc; i++) {
        //Parser parser(argv[i]);
        //KernelSelector parsed_sel(parser.kernel());
        //AtomicGroup g = system.select(parsed_sel);
        AtomicGroup g = loos::selectAtoms(system, argv[i]);
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
    traj->readFrame(num_skip);
  
    // loop over the remaining frames
    int frame = 0;
    while (traj->readFrame()) {
        // update coordinates
        traj->updateGroupCoords(system);

        // Compute the bin volume for normalization purposes
        GCoord box = system.periodicBox();
        double bin_volume = bin_width * box.x() * box.y();

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
                else
                    throw(runtime_error(
                                    "must choose either charge, mass or elec"));
                double z = pa->coords().z();

                if ( (z > min_z) && (z < max_z) ) {
                    int bin = (int)( (z - min_z) / bin_width );
                    dists[i][bin] += weight/bin_volume;
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
