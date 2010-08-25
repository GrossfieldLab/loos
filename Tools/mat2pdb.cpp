/*
  mat2pdb.cpp
  
  Takes 3 columns from a matrix and turns that into coordinates for pseudoatoms in a PDB
*/





/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010 Tod D. Romo
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
#include <boost/format.hpp>
#include <limits>

using namespace std;
using namespace loos;
namespace po = boost::program_options;


string matrix_name;
vector<uint> rows;
vector<uint> cols;
bool connect;
double scale;

string residue_name("SVD");
string segment_name("");



void parseOptions(int argc, char *argv[]) {
  try {

    string col_desc;
    string row_desc;

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help,h", "Produce this help message")
      ("cols,c", po::value<string>(&col_desc)->default_value("0,1,2"), "Columns to use")
      ("rows,r", po::value<string>(&row_desc)->default_value("all"), "Rows to use")
      ("scale,s", po::value<double>(&scale)->default_value(100.0), "Scale coordinates")
      ("connect,C", po::value<bool>(&connect)->default_value(false), "Connect sequential atoms with bonds");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("matrix", po::value<string>(&matrix_name), "Matrix filename");
    
    po::options_description command_line;
    command_line.add(generic).add(hidden);
    
    po::positional_options_description p;
    p.add("matrix", 1);

    
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !vm.count("matrix")) {
      cerr << "Usage- " << argv[0] << " [options] matrix-name >output.pdb\n";
      cerr << generic;

      exit(-1);
    }

    cols = parseRangeList<uint>(col_desc);
    if (row_desc != "all")
      rows = parseRangeList<uint>(row_desc);
  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}



int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  DoubleMatrix V;
  readAsciiMatrix(matrix_name, V);

  if (rows.empty())
    for (uint j=0; j<V.rows(); ++j)
      rows.push_back(j);

  AtomicGroup model;
  for (uint j=0; j<rows.size(); ++j) {
    GCoord c(V(rows[j], cols[0]), V(rows[j], cols[1]), V(rows[j], cols[2]));
    c *= scale;
    
    pAtom p(new Atom(j+1, "CA", c));
    p->resid(j+1);
    p->resname(residue_name);
    p->segid(segment_name);
    p->bfactor( (100.0 * j) / rows.size() );

    model.append(p);
  }

  if (connect)
    for (uint j=0; j<rows.size()-1; ++j)
      model[j]->addBond(model[j+1]);

  PDB pdb = PDB::fromAtomicGroup(model);
  pdb.remarks().add(hdr);
  cout << pdb;
}
