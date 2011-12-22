/*
  enmovie

  Elastic Network MOde VIsualizEr


  Usage:
    enmovie [options] model-name eigenvector-matrix output-name

  Notes:
    use the "--help" option for more information about how to run...

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008 Tod D. Romo
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
#include <iterator>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <cassert>

using namespace std;
using namespace loos;
namespace po = boost::program_options;


typedef Math::Matrix<double> Matrix;


// Globals, yuck!
string selection;
string supersel;
string eigvec_name, eigval_name;
string output_name;
string model_name;
bool invert_eigvals;
bool tag;
int steps;

vector<uint> modes;
vector<double> scales;


void parseOptions(int argc, char *argv[]) {
  string scale;
  string mode_string;

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection)->default_value("name == 'CA'"), "Selection used in computing the ANM")
      ("nsteps,n", po::value<int>(&steps)->default_value(100), "Number of steps to use in generating the trajectory")
      ("modes,m", po::value<string>(&mode_string)->default_value("0"), "List of modes to use")
      ("scales,S", po::value<string>(&scale)->default_value("1.0"), "Scale the displacements by this factor")
      ("invert,i", po::value<bool>(&invert_eigvals)->default_value(true), "Invert the eigenvalues for scaling")
      ("eigs,e", po::value<string>(&eigval_name), "Scale using corresponding eigenvalues from this file")
      ("superset,u", po::value<string>(&supersel)->default_value("all"), "Superset to use for frames in the output")
      ("tag,t", po::value<bool>(&tag)->default_value(false), "Tag ENM atoms with 'E' alt-loc");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("eigvec", po::value<string>(&eigvec_name), "Eigenvector name")
      ("outprefix", po::value<string>(&output_name), "Prefix for output filenames");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("eigvec", 1);
    p.add("outprefix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("eigvec") && vm.count("outprefix")) ) {
      cerr << "Usage- enmovie [options] model-name eigenvector-filename output-prefix\n";
      cerr << generic;
      exit(-1);
    }

    modes = parseRangeList<uint>(mode_string);
    if (modes.empty()) {
      cerr << "Error- invalid mode list.\n";
      exit(-1);
    }

    scales = parseRangeList<double>(scale);
    if (scales.size() == 1)
      for (vector<uint>::iterator i = modes.begin() + 1; i != modes.end(); ++i)
        scales.push_back(scales[0]);

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}




// Assume the scales.size() == modes.sizes()
void scaleByEigenvalues(vector<double>& scales, const vector<uint>& modes, const string& fname) {
  DoubleMatrix M;
  readAsciiMatrix(fname, M);

  for (uint i=0; i<modes.size(); ++i) {
    if (modes[i] >= M.rows()) {
      cerr << "ERROR- mode index " << modes[i] << " exceeds eigenvalue matrix dimensions.\n";
      exit(-10);
    }
    scales[i] *= (invert_eigvals ? 1./M[modes[i]] : M[modes[i]]);
  }
}



int main(int argc, char *argv[]) {
  // Get our invocation header for future metadata and parse
  // command-line args...
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  // Get the structure and select out the appropriate atoms...
  AtomicGroup model = createSystem(model_name);
  AtomicGroup superset = selectAtoms(model, supersel);
  AtomicGroup subset = selectAtoms(superset, selection);
  uint natoms = subset.size();

  Matrix U;
  readAsciiMatrix(eigvec_name, U);
  uint m = U.rows();
  if (m != subset.size() * 3) {
    cerr << boost::format("ERROR- subset size (%d) does not match eigenvector matrix size (%d)\n")
      % (subset.size() * 3)
      % m;
    exit(-10);
  }

  if (!eigval_name.empty())
    scaleByEigenvalues(scales, modes, eigval_name);
  

  // Instantiate a DCD writer...
  DCDWriter dcd(output_name + ".dcd");

  // We'll step along the Eigenvectors using a sin wave as a final scaling...
  double delta = 2*PI/steps;

  // Loop over requested number of frames...
  for (int frameno=0; frameno<steps; frameno++) {
    double k = sin(delta * frameno);

    // Have to make a copy of the atoms since we're computing a
    // displacement from the model structure...
    AtomicGroup frame = superset.copy();
    AtomicGroup frame_subset = selectAtoms(frame, selection);

    // Loop over all requested modes...
    for (uint j = 0; j<modes.size(); ++j) {
      // Loop over all atoms...
      for (uint i=0; i<natoms; i++) {
	GCoord c = frame_subset[i]->coords();

	// This gets the displacement vector for the ith atom of the
	// ci'th mode...
        GCoord d( U(i*3, modes[j]), U(i*3+1, modes[j]), U(i*3+2, modes[j]) );

	c += k * scales[j] * d;

	// Stuff those coords back into the Atom object...
	frame_subset[i]->coords(c);
        if (tag)
          frame_subset[i]->chainId("E");
      }
    }

    if (k == 0) {
      // Write out the selection, converting it to a PDB
      string outpdb(output_name + ".pdb");
      ofstream ofs(outpdb.c_str());
      PDB pdb;
      pdb = PDB::fromAtomicGroup(frame);
      pdb.remarks().add(header);
      pdb.renumber();
      pdb.clearBonds();
      ofs << pdb;
    }


    // Now that we've displaced the frame, add it to the growing DCD trajectory...
    dcd.writeFrame(frame);
  }
}
