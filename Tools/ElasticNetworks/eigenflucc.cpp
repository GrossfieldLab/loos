/*
  eigenflucc

  (c) 2010 Tod D. Romo, Grossfield Lab
           Department of Biochemistry
           University of Rochster School of Medicine and Dentistry


  Predict isotropic B-factors from a set of eigenpairs...
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
#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;

namespace po = boost::program_options;


// --------------- GLOBALS

bool verbose;
bool pca_input;
vector<uint> modes;
string eigvals_name, eigvecs_name, pdb_name;
string out_name;
double scale;
string selection;

const double kB = 6.950356e-9;  // \AA^{-1} K


void parseArgs(int argc, char *argv[]) {
  
  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("verbose,v", po::value<bool>(&verbose)->default_value(false), "Verbose output")
      ("selection,s", po::value<string>(&selection)->default_value("name == 'CA'"), "Selection used to make the ENM (only when altering a PDB)")
      ("pdb,p", po::value<string>(&pdb_name), "Alter the B-factors in a PDB")
      ("outpdb,o", po::value<string>(&out_name), "Filename to output PDB to")
      ("modes,m", po::value< vector<string> >(), "Modes to use (default is all)")
      ("scale,S", po::value<double>(&scale)->default_value(1.0), "Scaling factor to apply to eigenvalues")
      ("pca,P", po::value<bool>(&pca_input)->default_value(false), "Eigenpairs come from PCA, not ENM");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("eigvals", po::value<string>(&eigvals_name), "Eigenvalues filename")
      ("eigvecs", po::value<string>(&eigvecs_name), "Eigenvectors filename");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("eigvals", 1);
    p.add("eigvecs", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << "Usage- " << argv[0] << " [options] eigenvectors eigenvalues\n";
      cout << generic;
      exit(0);
    }

    if (vm.count("modes")) {
      vector<string> mode_list = vm["modes"].as< vector<string> >();
      modes = parseRangeList<uint>(mode_list);
    }

    if (!pdb_name.empty())
      if (out_name.empty())
        out_name = findBaseName(pdb_name) + "-ef.pdb";


  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

}



int main(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);
  parseArgs(argc, argv);

  DoubleMatrix eigvals;
  readAsciiMatrix(eigvals_name, eigvals);

  DoubleMatrix eigvecs;
  readAsciiMatrix(eigvecs_name, eigvecs);

  if (modes.empty())
    for (uint i=0; i<eigvals.rows(); ++i)
      modes.push_back(i);

  uint n = modes.size();
  uint m = eigvecs.rows();

  // V = m x n matrix of eigenvectors
  DoubleMatrix V(m, n);
  for (uint i=0; i<n; ++i)
    for (uint j=0; j<m; ++j)
      V(j, i) = eigvecs(j, modes[i]);

  // Make a copy and scale by the appropriate eigenvalue,
  // remembering that eigenvalues are inverted for ENM,
  // or squaring and not inverting in the case of PCA
  DoubleMatrix VS = V.copy();
  for (uint i=0; i<n; ++i) {
    double e = eigvals[modes[i]];
    double s = (pca_input) ? (scale * e * e): (scale / e);
    for (uint j=0; j<m; ++j)
      VS(j, i) *= s;
  }
  
  // U = Covariance matrix
  DoubleMatrix U = loos::Math::MMMultiply(VS,V,false,true);

  // B-factors come from trace of diagonal 3x3 superblocks
  vector<double> B;
  double prefactor = 8.0 *  M_PI * M_PI / 3.0;
  for (uint i=0; i<m; i += 3) {
    double b = prefactor * (U(i,i) + U(i+1, i+1) + U(i+2, i+2));
    B.push_back(b);
    cout << boost::format("%-8d %g\n") % (i/3) % b;
  }

  if (!pdb_name.empty()) {
    AtomicGroup model = createSystem(pdb_name);
    AtomicGroup subset = selectAtoms(model, selection);

    if ((unsigned int)subset.size() != B.size()) {
      cerr << boost::format("Error- selection has %d atoms, but %d were expected.\n") % subset.size() % B.size();
      exit(-10);
    }

    for (uint i=0; i<B.size(); ++i)
      subset[i]->bfactor(B[i]);

    PDB pdb = PDB::fromAtomicGroup(model);
    pdb.remarks().add(hdr);
    ofstream ofs(out_name.c_str());
    ofs << pdb;
  }

}

