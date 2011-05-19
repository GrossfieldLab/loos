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


using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;


typedef Math::Matrix<float, Math::ColMajor> Matrix;

// @cond TOOL_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() :
    scale(1.0),
    log(false),
    index(0),
    mapname("")
  { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("map", po::value<string>(&mapname), "Use a map file to select atoms to color")
      ("scale", po::value<double>(&scale)->default_value(scale), "Scale magnitudes by this amount")
      ("log", po::value<bool>(&log)->default_value(log),"Log-scale the output")
      ("index", po::value<int>(&index)->default_value(index), "SVD Term index to use");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("map='%s',scale=%f,log=%d,index=%d")
      % mapname
      % scale
      % log
      % index;
    return(oss.str());
  }


  double scale;
  bool log;
  int index;
  string mapname;

};
// @endcond



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
  
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::OutputPrefix* popts = new opts::OutputPrefix;
  opts::ModelWithCoords* mopts = new opts::ModelWithCoords;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(popts).add(mopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);


  AtomicGroup model = mopts->model;

  vector<pAtom> atoms;
  if (topts->mapname == "") {
    for (int i=0; i<model.size(); i++)
      atoms.push_back(model[i]);
  } else {
    vector<int> indices = readMap(topts->mapname);
    atoms = getAtoms(model, indices);
  }

  Matrix U;
  readAsciiMatrix(popts->prefix + "_U.asc", U);
  uint m = U.rows();
  uint n = U.cols();

  cerr << "Read in " << m << " x " << n << " matrix from " << popts->prefix + "_U.asc" << endl;

  if (m % 3 != 0) {
    cerr << "Error- dimensions of LSVs are bad.\n";
    exit(-11);
  }

  Matrix S;
  readAsciiMatrix(popts->prefix + "_s.asc", S);
  cerr << "Read in " << S.rows() << " singular values from " << popts->prefix + "_s.asc" << endl;

  pAtom pa;
  AtomicGroup::Iterator iter(model);
  while (pa = iter())
    pa->bfactor(0.0);

  double sval = S[topts->index];
  uint i;
  uint j;
  bool warned = false;
  for (i=j=0; j<m; j += 3) {
    GCoord c(U(j, topts->index), U(j+1, topts->index), U(j+2, topts->index));
    c *= sval;
    double b = topts->scale * c.length();
    if (topts->log)
      b = log(b);

    if (b<0.0) {
      if (!warned) {
        cerr << "WARNING - There are negative B-values.  These will be reset to zero.\n";
        warned = true;
      }
      b = 0.0;
    }
    atoms[i++]->bfactor(b);
  }

  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(header);
  cout << pdb;
}
