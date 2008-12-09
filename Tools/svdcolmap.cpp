/*
  svdcolmap.cpp

  Takes the magnitude of a left singular vector and maps this onto a
  PDB's B-values.
*/



/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
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
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;
using namespace loos;



typedef Math::Matrix<float, Math::ColMajor> Matrix;


struct Globals {
  string mapname;
  double scale;
  bool log;
  int index;
  string model_name;
  string prefix;
};

Globals globals;

void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("map,m", po::value<string>(&globals.mapname), "Use a map file to select atoms to color")
      ("scale,s", po::value<double>(&globals.scale)->default_value(1.0), "Scale magnitudes by this amount")
      ("log,l", po::value<bool>(&globals.log)->default_value(false),"Log-scale the output")
      ("index", po::value<int>(&globals.index)->default_value(0), "SVD Term index to use");


    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&globals.model_name), "Model filename")
      ("prefix", po::value<string>(&globals.prefix), "Prefix for SVD files");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("prefix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("prefix"))) {
      cerr << "Usage- svdcolmap [options] model-name svd-prefix\n";
      cerr << generic;
      exit(-1);
    }
  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


vector<int> readMap(const string& fname) {
  ifstream ifs(fname.c_str());
  if (!ifs) {
    cerr << "Error- cannot open " << fname << " for reading.\n";
    exit(-10);
  }

  vector<int> indices;
  int i, j, k;
  while (ifs >> i >> j >> k)
    indices.push_back(j);

  return(indices);
}


vector<pAtom> getAtoms(AtomicGroup& grp, const vector<int>& ids) {
  vector<pAtom> atoms;
  vector<int>::const_iterator i;

  for (i = ids.begin(); i != ids.end(); ++i) {
    pAtom pa(grp.findById(*i));
    if (!pa) {
      cerr << "Error- unable to find atom-id " << *i << endl;
      exit(-20);
    }
    atoms.push_back(pa);
  }

  return(atoms);
}



int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);
  AtomicGroup model = createSystem(globals.model_name);

  vector<pAtom> atoms;
  if (globals.mapname == "") {
    for (int i=0; i<model.size(); i++)
      atoms.push_back(model[i]);
  } else {
    vector<int> indices = readMap(globals.mapname);
    atoms = getAtoms(model, indices);
  }

  Matrix U;
  readAsciiMatrix(globals.prefix + "_U.asc", U);
  uint m = U.rows();
  uint n = U.cols();

  cerr << "Read in " << m << " x " << n << " matrix from " << globals.prefix + "_U.asc" << endl;

  if (m % 3 != 0) {
    cerr << "Error- dimensions of LSVs are bad.\n";
    exit(-11);
  }

  Matrix S;
  readAsciiMatrix(globals.prefix + "_s.asc", S);
  cerr << "Read in " << S.rows() << " singular values from " << globals.prefix + "_s.asc" << endl;

  pAtom pa;
  AtomicGroup::Iterator iter(model);
  while (pa = iter())
    pa->bfactor(0.0);

  double sval = S[globals.index];
  uint i;
  uint j;
  for (i=j=0; j<m; j += 3) {
    GCoord c(U(j, globals.index), U(j+1, globals.index), U(j+2, globals.index));
    c *= sval;
    double b = globals.scale * c.length();
    if (globals.log)
      b = log(b);
    atoms[i++]->bfactor(b);
  }

  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(header);
  cout << pdb;
}
