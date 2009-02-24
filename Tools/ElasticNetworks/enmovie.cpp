/*
  enmovie

    (c) 2008,2009 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry

      Elastic Network MOde VIsualizEr


  Usage:
    enmovie [options] model-name prefix output-name

  Notes:
    use the "--help" option for more information about how to run...

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
string mode_string;
string input_prefix;
string output_name;
string model_name;
double scale;
int steps;
bool scale_by_sval(true);

vector<int> modes;


void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("selection,s", po::value<string>(&selection)->default_value("name == 'CA'"), "Selection used in computing the ANM")
      ("nsteps,n", po::value<int>(&steps)->default_value(100), "Number of steps to use in generating the trajectory")
      ("modes,m", po::value<string>(&mode_string)->default_value("0"), "List of modes to use")
      ("scale,S", po::value<double>(&scale)->default_value(1.0), "Scale the displacements by this factor")
      ("norelative", "Do NOT scale the displacements by the corresponding singular value");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("inprefix", po::value<string>(&input_prefix), "Prefix for input filenames")
      ("outprefix", po::value<string>(&output_name), "Prefix for output filenames");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model", 1);
    p.add("inprefix", 1);
    p.add("outprefix", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || !(vm.count("model") && vm.count("inprefix") && vm.count("outprefix")) ) {
      cerr << "Usage- enmovie [options] model-name input-prefix output-prefix\n";
      cerr << generic;
      exit(-1);
    }

    if (vm.count("norelative"))
      scale_by_sval = false;
    
    modes = parseRangeList(mode_string);
    if (modes.empty()) {
      cerr << "Error- invalid mode list.\n";
      exit(-1);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}



int main(int argc, char *argv[]) {
  // Get our invocation header for future metadata and parse
  // command-line args...
  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);

  // Get the structure and select out the appropriate atoms...
  AtomicGroup model = createSystem(model_name);
  AtomicGroup subset = selectAtoms(model, selection);
  int natoms = subset.size();

  cout << boost::format("Selected %u atoms out of %u from %s.\n") % subset.size() % model.size() % model_name;

  // Write out the selection, converting it to a PDB
  string outpdb(output_name + ".pdb");
  ofstream ofs(outpdb.c_str());
  PDB pdb;
  pdb = PDB::fromAtomicGroup(subset);
  pdb.remarks().add(header);
  ofs << pdb;

  Matrix U;
  readAsciiMatrix(input_prefix + "_U.asc", U);
  int m = U.rows();
  int n = U.cols();
  cout << boost::format("Read in a %dx%d matrix of left singular vectors.\n") % m % n;
  assert(m == 3*natoms && "Incorrect number of atoms in SVD versus selection");


  Matrix S;
  readAsciiMatrix(input_prefix + "_s.asc", S);
  int sm = S.rows();
  int sn = S.cols();
  cout << boost::format("Read in a %dx%d matrix of singular values.\n") % sm % sn;
  assert(sm == 3*natoms && "Incorrect number of atoms in SVD versus selection");

  // Since we've decomposed a hessian, the slowest modes have the 
  // smallest singular values, but we have to exclude the 6 lowest...
  vector<int>::iterator vi;
  for (vi = modes.begin(); vi != modes.end(); ++vi) {
    *vi = n - 7 - *vi;
    if (*vi < 0) {
      cerr << boost::format("ERROR - you requested an invalid mode (%d).\n") % (n - (*vi + 7));
      exit(-1);
    }
  }

  cout << "Transformed mode indices: ";
  copy(modes.begin(), modes.end(), ostream_iterator<int>(cout, ","));
  cout << endl;
  

  // Instantiate a DCD writer...
  DCDWriter dcd(output_name + ".dcd");

  // We'll step along the LSV's using a sin wave as a final scaling...
  double delta = 2*PI/steps;

  // Loop over requested number of frames...
  for (int frameno=0; frameno<steps; frameno++) {
    double k = sin(delta * frameno);

    // Have to make a copy of the atoms since we're computing a
    // displacement from the model structure...
    AtomicGroup frame = subset.copy();

    // Loop over all requested modes...
    vector<int>::const_iterator ci;
    for (ci = modes.begin(); ci != modes.end(); ++ci) {
      // Loop over all atoms...
      for (int i=0; i<natoms; i++) {
	GCoord c = frame[i]->coords();

	// This gets the displacement vector for the ith atom of the
	// ci'th mode...
        GCoord d( U(i*3, *ci), U(i*3+1, *ci), U(i*3+2, *ci) );

        if (scale_by_sval)
          d /= S[*ci];

	c += scale * k * d;

	// Stuff those coords back into the Atom object...
	frame[i]->coords(c);
      }
    }

    // Now that we've displaced the frame, add it to the growing DCD trajectory...
    dcd.writeFrame(frame);
  }
}
