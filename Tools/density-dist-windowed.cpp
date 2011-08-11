/*
  density-dist-windowed.cpp

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
using namespace loos;


void Usage()
    {
    cerr << "Usage: density-dist-windowed "
         << " system traj E|C|M num_frames_to_skip min_z max_z num_bins window_size"
         << " filename_prototype [extra_selection_1 ...]"
         << endl;
    cerr << "Note: the system file must specify the mass and charge" << endl;
    }

int main(int argc, char *argv[]) {

    cout << "# This program is now deprecated: we suggest you use density-dist"
         << " with the --window option instead."
         << endl;


    if ( (argc <= 1) || 
         (strncmp(argv[1], "-h", 2) == 0) ||
         (argc < 10) ) {
        Usage();
        exit(-1);
    }

    cout << "# " << invocationHeader(argc, argv) << endl;


    AtomicGroup system = createSystem(argv[1]); 
    pTraj traj = createTrajectory(argv[2], system);
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
    subsets.push_back(system);
    for (int i=10; i<argc; i++) {
        //Parser parser(argv[i]);
        //KernelSelector parsed_sel(parser.kernel());
        //AtomicGroup g = system.select(parsed_sel);
        AtomicGroup g = selectAtoms(system, argv[i]);
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
    if (num_skip > 0)
      traj->readFrame(num_skip-1);
  
    // loop over the remaining frames
    int frame = 0;
    while (traj->readFrame()) {
        // update coordinates
        traj->updateGroupCoords(system);

        // Loop over the subsets and compute charge distribution.
        // (first set is all atoms)
        for (unsigned int i=0; i<subsets.size(); i++) {
            pAtom pa;
            AtomicGroup::Iterator iter(subsets[i]);
            double weight=0.0;
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
