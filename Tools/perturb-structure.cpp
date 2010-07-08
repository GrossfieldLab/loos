/*
  perturb-structure.cpp

  Apply a random perturbation to a structure (random directions, fixed magnitude)
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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


#include <loos.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;


uint seed;
string selection;
string model_name;
double magnitude;



void parseArgs(int argc, char *argv[]) {
  
  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection)->default_value("all"), "Selection to perturb")
      ("seed,S", po::value<uint>(&seed)->default_value(0l), "Random number seed (0 = use current time)");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("magnitude", po::value<double>(&magnitude), "Magnitude");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("magnitude", 1);
    p.add("model", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << "Usage- " << argv[0] << " [options] magnitude model >output.pdb\n";
      cout << generic;
      exit(0);
    }

    if (seed == 0)
      randomSeedRNG();
    else {
      base_generator_type& rng = rng_singleton();
      rng.seed(seed);
    }
      

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}





int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseArgs(argc, argv);

  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection);

  subset.perturbCoords(magnitude);
  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(hdr);

  cout << pdb;

}
