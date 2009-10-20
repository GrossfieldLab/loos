/*
  alignto.cpp

  Aligns a trajectory to a reference model...
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo
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



using namespace loos;
using namespace std;
namespace po = boost::program_options;


string ref_name, ref_sel;
string model_name, model_sel, traj_name;
string transform_sel;
string out_name;

vector<uint> indices;





void fullHelp() {
  cout << "No extra help available at this time\n";
}



void parseOptions(int argc, char *argv[]) {
  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("reference,r", po::value<string>(&ref_sel)->default_value("name == 'CA'"), "Reference selection to align to")
      ("selection,s", po::value<string>(&model_sel)->default_value("name == 'CA'"), "Model selection to align with")
      ("transform,t", po::value<string>(&transform_sel)->default_value("all"), "Transform this selection")
      ("range,R", po::value< vector<string> >(), "Frames to align with");


    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("out", po::value<string>(&out_name), "Output filename")
      ("ref", po::value<string>(&ref_name), "Reference filename")
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filenames");


    po::options_description command_line;
    command_line.add(generic).add(hidden);

    po::positional_options_description p;
    p.add("out", 1);
    p.add("ref", 1);
    p.add("model", 1);
    p.add("traj", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || vm.count("fullhelp") || !(vm.count("model") && vm.count("traj") && vm.count("ref") && vm.count("out"))) {
      cerr << "Usage- " << argv[0] << " [options] output-name reference model trajectory\n";
      cerr << generic;
      if (vm.count("fullhelp"))
        fullHelp();
      exit(-1);
    }

    if (vm.count("range")) {
      vector<string> ranges = vm["range"].as< vector<string> >();
      indices = parseRangeList<uint>(ranges);
    }

  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  parseOptions(argc, argv);
  AtomicGroup reference = createSystem(ref_name);
  AtomicGroup ref_subset = selectAtoms(reference, ref_sel);
  
  AtomicGroup model = createSystem(model_name);
  AtomicGroup model_subset = selectAtoms(model, model_sel);
  AtomicGroup model_xform = selectAtoms(model, transform_sel);

  pTraj traj = createTrajectory(traj_name, model);

  DCDWriter dcdout(out_name);
  dcdout.setTitle(hdr);

  while (traj->readFrame()) {
    traj->updateGroupCoords(model);
    
    GMatrix M = model_subset.superposition(ref_subset);
    XForm W(M);
    model_xform.applyTransform(W);
    
    dcdout.writeFrame(model);
  }
}
