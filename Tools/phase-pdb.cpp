/*
  phase_pdb

  Takes 3 colums from a matrix and sticks them into the coordinates for a fake PDB.
  This is an easy way to visualize PCA results in 3D using your
  favorite molecular visualization software...
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


vector<uint> columns;
vector<double> scales;
uint chunksize;
string segid_fmt;
string residue_name;
string atom_name;

string matrix_name;


void parseArgs(int argc, char *argv[]) {
  
  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("segid,S", po::value<string>(&segid_fmt)->default_value("P%03d"), "Segid printf format")
      ("atom,a", po::value<string>(&atom_name)->default_value("CA"), "Atom name to use")
      ("residue,r", po::value<string>(&residue_name)->default_value("SVD"), "Residue name to use")
      ("column,c", po::value< vector<uint> >(&columns), "Columns to use (default are first 3)")
      ("scales,s", po::value< vector<double> >(&scales), "Scale columns (default is 10,10,10)")
      ("chunk,C", po::value<uint>(&chunksize)->default_value(0), "Divide vector into chunks by these number of rows"); 

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("matrix", po::value<string>(&matrix_name), "Matrix file");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("matrix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("matrix"))) {
      cout << "Usage- " << argv[0] << " [options] matrix >output.pdb\n";
      cout << generic;
      exit(0);
    }

    if (columns.empty())
      for (uint i=0; i<3; ++i)
        columns.push_back(i);

    if (scales.empty())
      for (uint i=0; i<3; ++i)
        scales.push_back(10.0);
    else if (scales.size() == 1)
      for (uint i=0; i<2; ++i)
        scales.push_back(scales[0]);

    if (columns.size() != scales.size()) {
      cerr << "Error- number of columns selected does not equal number of scales\n";
      exit(-1);
    }

    if (columns.size() != 3) {
      cerr << "Error- must select 3 columns\n";
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
  parseArgs(argc, argv);

  RealMatrix A;
  readAsciiMatrix(matrix_name, A);
  uint m = A.rows();

  uint resid = 1;
  uint chunk = 1;
  AtomicGroup model;

  for (uint atomid = 0; atomid < m; ++atomid, ++resid) {
    if (chunksize && resid > chunksize) {
      for (uint i=atomid - chunksize + 1; i<atomid; ++i)
        model[i]->addBond(model[i-1]);

      resid = 1;
      ++chunk;
    }

    GCoord c(scales[0] * A(atomid, columns[0]),
             scales[1] * A(atomid, columns[1]),
             scales[2] * A(atomid, columns[2]));
    pAtom pa(new Atom(atomid+1, atom_name, c));
    pa->resid(resid);
    pa->resname(residue_name);

    ostringstream segstream;
    segstream << boost::format(segid_fmt) % chunk;
    pa->segid(segstream.str());

    model.append(pa);
  }

  if (chunksize && resid > chunksize)
    for (uint i=m - chunksize + 1; i<m; ++i)
      model[i]->addBond(model[i-1]);
  else
    for (uint i=1; i<m; ++i)
      model[i]->addBond(model[i-1]);
  
  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(hdr);
  cout << pdb;
}
