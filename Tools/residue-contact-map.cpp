/*
  contact-map
  
  Generates a heat map of contacts between selected residues for a trajectory...
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo
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


typedef vector<AtomicGroup>   vGroup;

// @cond TOOL_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() :
    use_centers(false)
  { }

  void addGeneric(po::options_description& o) {
    o.add_options()
      ("centers", po::value<bool>(&use_centers)->default_value(false), "Use center of mass of residues for distance");
  }

  string print() const {
    ostringstream oss;

    oss << "centers=" << use_centers;
    return(oss.str());
  }

  bool use_centers;
};
// @endcond




void accumulateFrameUsingCenters(DoubleMatrix& M, const vGroup& residues, const double threshold) {
  
  for (uint j=1; j<residues.size(); ++j) {
    GCoord v(residues[j].centerOfMass());

    for (uint i=0; i<j; ++i)
      if (v.distance2(residues[i].centerOfMass()) <= threshold) {
        M(j, i) += 1;
        M(i, j) += 1;
      }
  }

  for (uint i=0; i<residues.size(); ++i)
    M(i, i) += 1;
}


void accumulateFrameUsingAllAtoms(DoubleMatrix& M, const vGroup& residues, const double threshold) {
  
  for (uint j=1; j<residues.size(); ++j) {
    
    for (uint i=0; i<j; ++i) {
      bool flag = false;
      for (AtomicGroup::const_iterator a = residues[j].begin(); a != residues[j].end() && !flag; ++a)
        for (AtomicGroup::const_iterator b = residues[i].begin(); b != residues[i].end(); ++b)
          if ((*a)->coords().distance2((*b)->coords()) <= threshold) {
            flag = true;
            break;
          }
      if (flag) {
        M(j, i) += 1;
        M(i, j) += 1;
      }
    }
  }

  for (uint i=0; i<residues.size(); ++i)
    M(i, i) += 1;
}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* bopts = new opts::BasicOptions;
  opts::BasicSelection* sopts = new opts::BasicSelection;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  ToolOptions* topts = new ToolOptions;
  opts::RequiredArguments* ropts = new opts::RequiredArguments("threshold", "Distance threshold for contacts");

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(tropts).add(topts).add(ropts);
  if (!options.parse(argc, argv))
    exit(-1);

  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  vector<uint> indices = tropts->frameList();

  double thresh = parseStringAs<double>(ropts->value("threshold"));
  thresh *= thresh;

  AtomicGroup subset = selectAtoms(model, sopts->selection);
  vGroup residues = subset.splitByResidue();

  DoubleMatrix M(residues.size(), residues.size());
  for (vector<uint>::iterator i = indices.begin(); i != indices.end(); ++i) {
    traj->readFrame(*i);
    traj->updateGroupCoords(model);
    if (topts->use_centers)
      accumulateFrameUsingCenters(M, residues, thresh);
    else
      accumulateFrameUsingAllAtoms(M, residues, thresh);
  }

  for (ulong i=0; i<residues.size() * residues.size(); ++i)
    M[i] /= indices.size();

  writeAsciiMatrix(cout, M, hdr);
}
