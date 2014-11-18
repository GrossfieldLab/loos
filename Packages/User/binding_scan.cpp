/*
  binding_scan
  
  (c) 2009 Tod D. Romo, Grossfield Lab
  Department of Biochemistry
  University of Rochster School of Medicine and Dentistry

  Loop over a selection and compute a distance based score that tells you
  how much contact is occuring between the probe and selection.

*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008-2009, Tod D. Romo
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
#include <numeric>

using namespace std;
using namespace loos;
namespace po = boost::program_options;

typedef vector<AtomicGroup> vGroup;
typedef vector<vector <AtomicGroup> > list_vGroup; 

double inner_cutoff, outer_cutoff;

string probe_selection;
string model_name, traj_name;
string target_selection;
int skip;

void fullHelp(void) {
    cout << "Sorry... can't help you";

}

void parseOptions(int argc, char *argv[]) {
  try {

    po::options_description generic("Allowed options");
    generic.add_options()
      ("help", "Produce this help message")
      ("fullhelp", "Even more help")
      ("probe,p", po::value<string>(&probe_selection)->default_value("segname =~ 'Rhod'"), "Main selection")
      ("skip,s", po::value<int>(&skip)->default_value(0), "Frames to skip");

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("model", po::value<string>(&model_name), "Model filename")
      ("traj", po::value<string>(&traj_name), "Trajectory filename")
      ("target", po::value<string>(&target_selection), "Target selection");
    
    po::options_description command_line;
    command_line.add(generic).add(hidden);
    
    po::positional_options_description p;
    p.add("model", 1);
    p.add("traj", 1);
    p.add("target", -1);

    
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || vm.count("fullhelp") || !(vm.count("model") && vm.count("traj") && !target_selection.empty())) {
      cerr << "Usage- " << argv[0] << " [options] model-name trajectory-name target\n";
      cerr << generic;
      if (vm.count("fullhelp"))
        fullHelp();
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
  parseOptions(argc, argv);

  AtomicGroup model = createSystem(model_name);
  pTraj traj = createTrajectory(traj_name, model);

  AtomicGroup probe = selectAtoms(model, probe_selection);
  vGroup probe_residues = probe.splitByResidue();

  cout << "# " << hdr << endl;
  // selections for targets
  AtomicGroup target = selectAtoms(model, target_selection);
  vGroup targets = target.splitByMolecule();

  // initialize scoring matrix
  vector<double> residue_score;
  for (unsigned int j = 0; j < probe_residues.size(); j++){

      residue_score.push_back(0);

    }


int frame_count = 0;

traj->readFrame(skip);

// Loop over trajectory
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);

    // Loop over residues
    for (unsigned int i = 0; i < residue_score.size(); i++){
        
        // Loop over molecules in target
        for (unsigned int j = 0; j < targets.size(); j++){
            double temp_score = 0.0;
            // Loop over atoms in target
            for (unsigned int k = 0; k < targets[j].size(); k++){
              
                // Loop over atoms in residue
                for (unsigned int l = 0; l < probe_residues[i].size(); l++){
                
                    // Calculate packing score for each atom and add to residue
                    double v = (targets[j][k]->coords()).distance(probe_residues[i][l]->coords());
                    temp_score = temp_score + (1 / (v*v*v*v*v*v));

                }

            }
            // Add score sum to running total for a given residue
            residue_score[i] = residue_score[i] + temp_score;
        }

    }

    frame_count++;
    // end while loop for frames
  }

  cout << "#Residue\tScore\tToAvg" << endl;

  double avg = 0;
  for (unsigned int a = 0; a < residue_score.size(); a++){

    avg = avg + (residue_score[a] / frame_count);

  }
  avg = avg / residue_score.size();
  cout << "#Avg " << avg << endl;
  // walk through matrix and normalize data by frame count & mol number and output
  for (unsigned int x = 0; x < residue_score.size(); x++){

      //double min = *min_element(residue_score.begin(), residue_score.end()) / frame_count;
      // normalized by frame and against average
      double temp_output = (residue_score[x]) / frame_count;
      double temp_normalized = (residue_score[x] / avg) / frame_count;

      // output current residue number
      cout << (probe_residues[x])[0]->resid() << "\t" 
           << boost::format("  %8.8f") % temp_output << "\t"
           << boost::format("  %8.8f") % temp_normalized << "\t"
           << endl;

    }

    //string sig;
    // find which score is x times greater than the others
    // unnecessary right now
    //for (unsigned int z = 0; z < score_matrix[x].size(); z++){

    //  for (unsigned int w = z+1; w < score_matrix[x].size(); w++){

    //    if (score_matrix[x][z] > threshold * score_matrix[x][w])
    //        sig = (targets[z])[0]->name();
    //    if (score_matrix[x][w] > threshold * score_matrix[x][z])
    //        sig = (targets[w])[0]->name();

    //  }

    //}

    //cout << "  " << sig;

    cout << endl;

  

// end main
} 

