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
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;

// Globals

enum CalculationType { ELECTRON, CHARGE, MASS };

CalculationType calc_type;
bool symmetrize;
bool auto_all;
uint skip, nbins;
uint window=0;
double min_z, max_z;
string model_name, traj_name, file_name_proto;
vector<string> selections;

void parseOptions(int argc, char *argv[]) {
  string calc_type_desc;

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("auto", po::value<bool>(&auto_all)->default_value(true), "Automatically compute system density")
      ("skip", po::value<uint>(&skip)->default_value(0), "Number of starting frames to skip")
      ("zsymmetry", po::value<bool>(&symmetrize)->default_value(false), "Symmetric with respect to Z")
      ("type", po::value<string>(&calc_type_desc)->default_value("electron"), "Calculation type (mass, charge, electron)")
      ("window", po::value<uint>(&window), "Window size (in frames) for time series")
      ("filename", po::value<string>(&file_name_proto), "Root of file name for windowed time series")
      ;
      

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "model")
      ("traj", po::value<string>(&traj_name), "trajectory")
      ("nbins", po::value<uint>(&nbins), "nbins")
      ("minz", po::value<double>(&min_z), "minz")
      ("maxz", po::value<double>(&max_z), "maxz")
      ("selection", po::value< vector<string> >(&selections), "selections");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("minz", 1);
    p.add("maxz", 1);
    p.add("nbins", 1);
    p.add("selection", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("traj") && vm.count("nbins") && vm.count("minz") && vm.count("maxz"))) {
      cerr << "Usage- " << argv[0] << " [options] model trajectory min-z max-z nbins [selection ...] >output\n";
      cerr << generic;
      exit(-1);
    }
    
    if (toupper(calc_type_desc[0]) == 'C')
      calc_type = CHARGE;
    else if (toupper(calc_type_desc[0]) == 'E')
      calc_type = ELECTRON;
    else if (toupper(calc_type_desc[0]) == 'M')
      calc_type = MASS;
    else {
      cerr << "Error- unknown calculation type '" << calc_type_desc << "' (should be either charge, mass or electron)\n";
      exit(-1);
    }

    if (vm.count("filename") != vm.count("window")) {
        cerr << "If you want windowed time series, you must specify both a window size and a filename root"
             << endl;
        exit(-1);
    }


  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);


  cout << "# " << hdr << endl;


  AtomicGroup system = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, system);

  // density from each selection
  vector<AtomicGroup> subsets;
  if (auto_all)
    subsets.push_back(system);
  for (vector<string>::iterator i = selections.begin(); i != selections.end(); ++i)
    subsets.push_back(selectAtoms(system, *i));

  // Verify properties...
  for (vector<AtomicGroup>::iterator i = subsets.begin(); i != subsets.end(); ++i) {
    bool ok;
    switch(calc_type) {
    case CHARGE: ok = i->allHaveProperty(Atom::chargebit); break;
    case ELECTRON: ok = (i->allHaveProperty(Atom::anumbit) && i->allHaveProperty(Atom::chargebit)); break;
    case MASS: ok = i->allHaveProperty(Atom::massbit); break;
    default: cerr << "Internal Error- unknown calculation type\n"; exit(-1);
    }
    if (!ok) {
      cerr << "WARNING- system is missing properties required for calculation type.  Default values will be used where possible.\n";
      break;
    }
  }
  
  double bin_width = (max_z - min_z) / nbins;

  // Create and zero out the total charge_distribution
  vector< vector<double> > dists(subsets.size(), vector<double>(nbins, 0.0));
  // If we do windowed time series too, we'll need to zero out as we go, which 
  // means we'll need separate storage for the sum
  vector< vector<double> > cum_dists(subsets.size(), vector<double>(nbins, 0.0));

  // skip the equilibration frames
  traj->readFrame(skip);
  
  // loop over the remaining frames
  uint frame = 0;
  while (traj->readFrame()) {
    // update coordinates
    traj->updateGroupCoords(system);

    // Compute the bin volume for normalization purposes
    GCoord box = system.periodicBox();
    double bin_volume = bin_width * box.x() * box.y();

    // Loop over the subsets and compute charge distribution.
    // (first set is all atoms)
    for (uint i=0; i<subsets.size(); i++) {
      double weight;
      for (AtomicGroup::iterator atom = subsets[i].begin(); 
                                 atom != subsets[i].end(); 
                                 ++atom) {
        switch(calc_type) {
        case CHARGE: weight = (*atom)->charge(); break;
        case MASS: weight = (*atom)->mass(); break;
        case ELECTRON: weight = (*atom)->atomic_number() - (*atom)->charge(); break;
        default:
          cerr << "ERROR- unknown calculation type\n";
          exit(-1);
        }
        double z = (*atom)->coords().z();
        if (symmetrize)
          z = abs(z);
        
        if ( (z > min_z) && (z < max_z) ) {
          uint bin = (int)( (z - min_z) / bin_width );
          dists[i][bin] += weight/bin_volume;
        }
      }
    }
    ++frame;

    // If windowed time series were requested, output them here
    if (window && (frame % window == 0)) {
        uint output_val = frame / window;

        // build the output file name
        char piece[128];
        snprintf(piece, 128, "_%d.dat", output_val);
        string file_name = file_name_proto + string(piece);
        ofstream outfile(file_name.c_str());
        if (!outfile.is_open() ) {
            throw(runtime_error("couldn't open output file"));
            }

        // write the header
        outfile << "# Z\tAllAtoms";
        for (uint i=1; i<subsets.size(); i++) {
          outfile << " Set(" << i <<") "; 
        }
        outfile << endl;

        // Normalize by the number of frames and
        // output the average density
        for (uint i=0; i<nbins; i++) {
            double z = (i+0.5)*bin_width + min_z;
            outfile << z << "\t";
            for (unsigned int j=0; j<subsets.size(); j++) {
                double c = dists[j][i] / window;
                outfile << c << "\t";
            }
            outfile << endl;
        }

        // Add the windowed density to the cumulative density and
        // zero out the distributions
        
        for (uint i=0; i<dists.size(); i++) {
            for (uint j=0; j<nbins; j++) {
                cum_dists[i][j] += dists[i][j];
            }

            // zero out the windowed dist
            dists[i].assign(dists[i].size(), 0.0);
        }
        
    }


  }

  // Normalize by the number of frames and output the average charge density
  cout << "# Z\t" << (auto_all ? "AllAtoms" : "");
  for (uint i=1; i<subsets.size(); i++) {
    cout << " Set(" << i <<") "; 
  }
  cout << endl;

  // If we didn't output windowed time series, we never used cum_dists, so we need 
  // to copy it over here
  // If we did, we may have some data in dists that hasn't been added to cum_dists yet
  if (!window) {
      cum_dists = dists;
  }
  else if (frame % window != 0){
     for (uint i=0; i<dists.size(); i++) {
         for (uint j=0; j<nbins; j++) {
             cum_dists[i][j] += dists[i][j];
         }
     }
  }

  for (uint i=0; i<nbins; i++) {
    double z = (i+0.5)*bin_width + min_z;
    cout << z << "\t"; 
    for (uint j=0; j<subsets.size(); j++) {
      double c = cum_dists[j][i] / frame;
      
      cout << c << "\t";
    }
    cout << endl;
  }


  
}
