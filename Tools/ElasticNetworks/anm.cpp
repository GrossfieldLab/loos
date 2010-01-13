/*
  anm

  (c) 2008 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

  Computes the anisotropic network model for a structure.  It does
  this by building a hessian for the structure, then computing the SVD
  of it and the corresponding pseudo-inverse (ignoring the 6 lowest
  modes).

  Usage:
    anm [selection string] radius model-name output-prefix

  Examples:
    anm 'resid >= 10 && resid <= 50 && name == "CA"' 15.0 foo.pdb foo

    This creates the following files:
      foo_H.asc      == The hessian
      foo_U.asc      == Left singular vectors
      foo_s.asc      == Singular values
      foo_V.asc      == Right singular vectors
      foo_Hi.asc     == Pseudo-inverse of H

  Notes:
    o The default selection (if none is specified) is to pick CA's
    o The output is in ASCII format suitable for use with Matlab/Octave/Gnuplot
      
  
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2009 Tod D. Romo
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
#include "hessian.hpp"

#include <limits>
#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace loos;
namespace po = boost::program_options;


typedef Math::Matrix<double, Math::ColMajor> Matrix;



// Globals...  Icky poo!

string selection;
string model_name;
string prefix;
double cutoff;

// Turns on parameter-free mode a la Yang et al, PNAS (2009) 106:12347
bool parameter_free;



void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection)->default_value("name == 'CA'"), "Which atoms to use for the network")
      ("pfree,P", po::value<bool>(&parameter_free)->default_value(false), "Use the parameter-free method rather than the cutoff")
      ("cutoff,c", po::value<double>(&cutoff)->default_value(15.0), "Cutoff distance for node contact");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("prefix", po::value<string>(&prefix), "Output prefix");

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
      cerr << "Usage- anm [options] model-name output-prefix\n";
      cerr << generic;
      exit(-1);
    }

    // Force the hessian to include all nodes...
    if (parameter_free)
      cutoff = numeric_limits<double>::max();

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection);
  cerr << boost::format("Selected %d atoms from %s\n") % subset.size() % model_name;

  Timer<WallTimer> timer;
  cerr << "Calculating hessian...";
  timer.start();
  Matrix H = hessian(subset, cutoff);
  timer.stop();
  cerr << boost::format(" done (%dx%d)\n") % H.rows() % H.cols();
  cerr << timer << endl;

  if (parameter_free)
    distanceWeight(H, subset);

  ScientificMatrixFormatter<double> sp(24, 18);

  writeAsciiMatrix(prefix + "_H.asc", H, header, false, sp);

  cerr << "Calculating SVD - ";
  timer.start();
  boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> result = svd(H);
  timer.stop();
  cerr << "done\n" << timer << endl;

  Matrix U = boost::get<0>(result);
  Matrix S = boost::get<1>(result);
  Matrix Vt = boost::get<2>(result);
  uint n = S.rows();

  reverseRows(S);
  reverseColumns(U);
  reverseRows(Vt);

  // Write out the LSVs (or eigenvectors)
  writeAsciiMatrix(prefix + "_U.asc", U, header, false, sp);
  writeAsciiMatrix(prefix + "_s.asc", S, header, false, sp);

  // Now go ahead and compute the pseudo-inverse...

  // Vt = Vt * diag(1./diag(S))
  // Remember, Vt is stored col-major but transposed, hence the
  // inverted indices...
  //
  // Note:  We have to toss the first 6 terms
  for (uint i=6; i<n; i++) {
    double s = 1.0/S[i];
    for (uint j=0; j<n; j++)
      Vt(i,j) *= s;
  }
  
  // Ki = Vt * U';
  // Again, Vt is internally transposed, so we have to specify
  // transposing it to sgemm in order to multiply the non-transposed
  // V...
  Matrix Hi = MMMultiply(Vt, U, true, true);
  writeAsciiMatrix(prefix + "_Hi.asc", Hi, header, false, sp);

}
