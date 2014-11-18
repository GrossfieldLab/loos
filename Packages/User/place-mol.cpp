/*
  place-mol.cpp

  (c) 2011 Joshua N. Horn, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Given two selections, will place selection 1 the specified
  distance away from selection 2 in the z-dimension. Useful
  for creating structures for insertion.

  Usage:  place-mol.cpp model1 model2 sel1 sel2
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008,2009 Tod Romo
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



#include <cmath>
#include <loos.hpp>
#include <limits>
#include <boost/program_options.hpp>
#include <sys/time.h>

using namespace std;
using namespace loos;
namespace po = boost::program_options;

// Globals
string model1_name, model2_name;
string selection1, selection2;
double dist = 0;
double delete_dist = 2;
bool randomly_rotate;

void parseOptions(int argc, char *argv[]) {

  try {
    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("rotate,r", po::value<bool>(&randomly_rotate)->default_value(0), "Randomly rotate selection 1.")
      ("deletion_distance,d", po::value<double>(&delete_dist)->default_value(2), "Distance to delete overlap.");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model1", po::value<string>(&model1_name), "model1")
      ("model2", po::value<string>(&model2_name), "model2")
      ("selection1", po::value< string >(&selection1), "selection1")
      ("selection2", po::value< string >(&selection2), "selection2")
      ("distance", po::value< double >(&dist), "distance");

    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("model1", 1);
    p.add("model2", 1);
    p.add("selection1", 1);
    p.add("selection2", 1);
    p.add("distance", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") ) {
      cerr << "Usage- " << argv[0] << " [options] model trajectory sel-1 sel-2 [sel-3 ...] >output\n";
      cerr << generic;
      exit(-1);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}



int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);
  parseOptions(argc, argv);
  cout << "# " << header << endl;

  AtomicGroup model1 = createSystem(model1_name);
  AtomicGroup model2 = createSystem(model2_name);

  AtomicGroup sel1 = selectAtoms(model1, selection1);
  AtomicGroup sel2 = selectAtoms(model2, selection2);

  model1.periodicBox(model2.periodicBox());
  double z_translate = (sel2.centroid()).z() - (sel1.centroid()).z() + dist;
  double y_center = (sel1.centroid()).y() - (sel2.centroid()).y();
  double x_center = (sel1.centroid()).x() - (sel2.centroid()).x();
  GCoord translation((-x_center),(-y_center),z_translate);
  model1.translate(translation);
  
  if (randomly_rotate){
      // translate to origin
    GCoord temp = sel1.centerAtOrigin();

    // Alan's code for a random number generator
    boost::mt19937 rng;
    boost::uniform_real<double> anglemap(-180.0, 180.0);
    // seed the system RNG from the time, use that to seed the 
    // individual RNGs
    timeval tim;
    gettimeofday(&tim, NULL);
    unsigned int seed = static_cast<unsigned int>(tim.tv_usec);
    srand(seed);
    rng.seed( rand() );
    boost::variate_generator<loos::base_generator_type&, boost::uniform_real<> > angles(rng, anglemap);

    // apply random rotation
    XForm rand_rot;
    rand_rot.rotate('x', angles());
    rand_rot.rotate('y', angles());
    rand_rot.rotate('z', angles());
    sel1.applyTransform(rand_rot);

    // translate back
    sel1.translate(temp);
  }

  // now we have our model positioned correctly, next we will delete overlapping atoms from model 2,
  // and merge the two files.
  AtomicGroup overlap = model2.within(delete_dist, model1);
  model2.remove(overlap);
  AtomicGroup merged = model1.merge(model2);
  merged.renumber(1, 1);

  // output PDB structure of model 1
  PDB pdb = PDB::fromAtomicGroup(merged);
  pdb.remarks().add(invocationHeader(argc, argv));
  cout << pdb;
}
