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



string fullHelpMessage(void) {
  string msg =
    "\n"
    "SYNOPSIS\n"
    "\tCalculate a contact \"heat-map\" between residues in a simulation.\n"
    "\n"
    "DESCRIPTION\n"
    "\n"
    "\tThis tool will break apart the selection into residues.  At each time point\n"
    "in the trajectory, it will determine if any residues are in contact with each\n"
    "other.  This will be accumulated over the trajectory and a matrix representing\n"
    "the fractional contacts will be written out.  This matrix can be visualized as\n"
    "a \"heat-map\" using octave/matlab or gnuplot.\n"
    "\tA contact can be defined in two different ways.  It can be defined as occuring when\n"
    "the distance between any two atoms less than or equal to the\n"
    "threshold given on the command line.  Alternatively, it can be defined as occuring when\n"
    "the distance between the centers of mass of the two residues is less than or equal\n"
    "to the threshold.\n"
    "\n"
    "EXAMPLES\n"
    "\n"
    "\tresidue-contact-map --select 'segid == \"PROT\"' \\n"
    "\t  model.pdb simulation.dcd 4.0 >contacts.asc\n"
    "This example defines a contact when any pair of atoms between a given two residues is\n"
    "less than or equal to the 4.0 Angstroms.  Only residues with segid \"PROT\" are used.\n"
    "\n"
    "\tresidue-contact-map --select 'resid <= 100' --centers 1 \\\n"
    "\t  model.pdb simulation.dcd 6.5 >contacts.asc\n"
    "This example defines a contact when the centers of mass between two residues is less than\n"
    "or equal two 6.5 Angstroms.  Only the first 100 residues are used.\n"
    "\n"
    "SEE ALSO\n"
    "\trmsds\n";

  return(msg);
}



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
      bool flag = true;
      for (AtomicGroup::const_iterator a = residues[j].begin(); a != residues[j].end() && flag; ++a)
        for (AtomicGroup::const_iterator b = residues[i].begin(); b != residues[i].end(); ++b)
          if ((*a)->coords().distance2((*b)->coords()) <= threshold) {
            flag = false;
            break;
          }
      if (!flag) {
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

  opts::BasicOptions* bopts = new opts::BasicOptions(fullHelpMessage());
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
