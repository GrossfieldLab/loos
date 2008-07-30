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
#include <getopt.h>
#include <cmath>

string mapname;
double scale=1.0;
bool logscale = false;


static struct option long_options[] = {
  {"scale", required_argument, 0, 's'},
  {"log", no_argument, 0, 'l'},
  {"map", required_argument, 0, 'm'},
  {0,0,0,0}
};

static const char* short_options = "s:m:l";

void show_help(void) {
  cout << "Usage- svdcolmap [--map=fname] [--scale=float] [--log] index pdb svd_prefix >output.pdb\n";
}


void parseOptions(int argc, char *argv[]) {
  int opt, idx;

  while ((opt = getopt_long(argc, argv, short_options, long_options, &idx)) != -1)
    switch(opt) {
    case 'm': mapname = string(optarg); break;
    case 's': scale = strtod(optarg, 0); break;
    case 'l': logscale = true; break;
    case 0: break;
    default: cerr << "Unknown option '" << (char)opt << "' - ignored.\n";
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

  parseOptions(argc, argv);
  if (argc - optind != 3) {
    show_help();
    exit(-1);
  }

  int index = atoi(argv[optind++]);
  PDB pdb(argv[optind++]);
  string prefix(argv[optind++]);
  RawAsciiReader<double> reader;

  vector<pAtom> atoms;
  if (mapname == "") {
    for (int i=0; i<pdb.size(); i++)
      atoms.push_back(pdb[i]);
  } else {
    vector<int> indices = readMap(mapname);
    atoms = getAtoms(pdb, indices);
  }

  MatrixReader<double>::Result lsvs = reader.read(prefix + "U.asc");
  double *U = boost::get<0>(lsvs);
  uint m = boost::get<1>(lsvs);
  uint n = boost::get<2>(lsvs);

  cerr << "Read in " << m << " x " << n << " matrix from " << prefix + "U.asc" << endl;

  if (m % 3 != 0) {
    cerr << "Error- dimensions of LSVs are bad.\n";
    exit(-11);
  }

  MatrixReader<double>::Result svals = reader.read(prefix + "s.asc");
  double *s = boost::get<0>(svals);
  cerr << "Read in " << boost::get<1>(svals) << " x " << boost::get<2>(svals) << " matrix from " << prefix + "s.asc" << endl;

  pAtom pa;
  AtomicGroup::Iterator iter(pdb);
  while (pa = iter())
    pa->bfactor(0.0);

  double *lsv = U + index * m;
  double sval = s[index];
  uint i;
  uint j;
  for (i=j=0; j<m; j += 3) {
    GCoord c(lsv[j], lsv[j+1], lsv[j+2]);
    c *= sval;
    double b = scale * c.length();
    if (logscale)
      b = log(b);
    atoms[i++]->bfactor(b);
  }

  cout << pdb;
}
