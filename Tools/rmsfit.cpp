/*
  rmsfit.cpp

  Superimposes one structure upon another
*/



/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008-2009 Tod D. Romo
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

namespace po = boost::program_options;

using namespace std;
using namespace loos;

// Globals...blech...

string source_name, target_name;
string source_selection, target_selection, apply_selection;
bool quiet = false;


void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options", 100);
    generic.add_options()
      ("help,h", "Produce this help message")
      ("apply,a", po::value<string>(&apply_selection)->default_value("all"), "Subset of source model to apply transformation to")
      ("source,s", po::value<string>(&source_selection)->default_value("name == 'CA'"), "Subset of the source model to align with")
      ("target,t", po::value<string>(&target_selection)->default_value("name == 'CA'"), "Subset of the target model to align with")
      ("quiet,Q", "Suppress writing out details of superposition");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("source_name", po::value<string>(&source_name), "Source model filename")
      ("target_name", po::value<string>(&target_name), "Target model filename");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("source_name", 1);
    p.add("target_name", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("source_name") && vm.count("target_name"))) {
      cerr << "Usage- rmsfit [options] source-model target-model >superimposed.pdb\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("quiet"))
      quiet = true;

  }

  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup source_model = createSystem(source_name);
  AtomicGroup target_model = createSystem(target_name);

  AtomicGroup source_subset = selectAtoms(source_model, source_selection);
  AtomicGroup apply_subset = selectAtoms(source_model, apply_selection);
  if (apply_selection != "all")
    apply_subset.clearBonds();

  AtomicGroup target_subset = selectAtoms(target_model, target_selection);

  if (source_subset.size() != target_subset.size()) {
    cerr << boost::format("ERROR - The source subset has %d atoms but the target subset has %d atoms.  The MUST be equal")
      % source_subset.size()
      % target_subset.size();
    exit(-10);
  }

  GMatrix M = source_subset.superposition(target_subset);
  XForm W(M);
  apply_subset.applyTransform(W);
  if (!quiet) {
    double d = source_subset.rmsd(target_subset);
    cerr << "Final RMSD = " << d << endl;
  }
  
  PDB pdb = PDB::fromAtomicGroup(apply_subset);
  pdb.remarks().add(hdr);
  cout << pdb;
}
