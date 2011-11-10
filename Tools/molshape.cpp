/*
  molwatch

  Computes size/shape/positional information for a selection over time...
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


using namespace std;
using namespace loos;
namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

string model_name, traj_name;
string selection;
int verbosity = 0;
bool split_by_mol;
bool split_by_segid;
bool zabs;


// @cond TOOLS_INTERNAL

class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() : split_by_mol(false), split_by_segid(false), zabs(false) { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("molecule", po::value<bool>(&split_by_mol)->default_value(split_by_mol), "Split by molecule")
      ("segid", po::value<bool>(&split_by_segid)->default_value(split_by_segid), "Split by segid")
      ("abs", po::value<bool>(&zabs)->default_value(zabs), "Use absolute Z-value");
  }

  string print() const {
    ostringstream oss;
    oss << boost::format("molecule=%d, segid=%d, zabs=%d") % split_by_mol % split_by_segid % zabs;
    return(oss.str());
  }


  bool split_by_mol, split_by_segid, zabs;
};

// @endcond





string split(const loos::Coord<double>& g) {
  stringstream ss;

  ss << setw(10) << g[0] << " " << g[1] << " " << g[2];
  return(ss.str());
}


void modifyZ(AtomicGroup& grp) {

  for (AtomicGroup::iterator i = grp.begin(); i != grp.end(); ++i) {
    GCoord c = (*i)->coords();
    c.z() = abs(c.z());
    (*i)->coords(c);
  }
}


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);
  
  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection("!hydrogen || !(segid == 'BULK' || segid == 'SOLV')");
  opts::BasicTrajectory* tropts = new opts::BasicTrajectory;
  ToolOptions* topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  cout << "# " << hdr << endl;
  
  AtomicGroup model = tropts->model;
  AtomicGroup subset = selectAtoms(model, sopts->selection);
  pTraj traj = tropts->trajectory;

  vector<AtomicGroup> objects;
  if (split_by_mol)
    objects = subset.splitByMolecule();
  else if (split_by_segid)
    objects = subset.splitByUniqueSegid();
  else
    objects.push_back(subset);

  cout << "# 1     2  3  4  5   6    7    8    9    10      11  12  13  14:16 17:19 20:22\n";
  cout << "# frame cX cY cZ Vol BoxX BoxY BoxZ rgyr pA1/pA2 pA1 pA2 pA3 (pV1) (pV2) (pV3)\n";

  

  uint t=tropts->skip;
  while (traj->readFrame()) {
    traj->updateGroupCoords(subset);
    if (zabs)
      modifyZ(subset);
    for (uint i=0; i<objects.size(); ++i) {
      GCoord c = objects[i].centroid();
      vector<GCoord> bdd = objects[i].boundingBox();
      GCoord box = bdd[1] - bdd[0];
      double vol = box[0] * box[1] * box[2];
      vector<GCoord> paxes = objects[i].principalAxes();
      double ratio = paxes[3][0] / paxes[3][1];
      double rgyr = objects[i].radiusOfGyration();
      
      cout << setw(10) << t <<  " " << split(c) << " " << vol << " " << split(box) << " " << rgyr << " ";
      cout << ratio << " " << split(paxes[3]) << " " << split(paxes[0]) << " " << split(paxes[1]) << " " << split(paxes[2]) << endl;
    }
    
    ++t;
  }
}
