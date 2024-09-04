/*
  Compute the hexagonal order parameter for a membrane


  Alan Grossfield
  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2024
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

#include <iostream>
#include <loos.hpp>

using namespace std;
using namespace loos;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMessage =
    // clang-format off
    "SYNOPSIS \n"
" \n"
"Compute the hexagonal order parameter for a membrane. \n"
" \n"
"DESCRIPTION \n"
"The hexagonal order parameter was proposed in \n"
"Nelson & Halperin, Phys. Rev. B 19, 2457–2484 (1979) \n"
" \n"
"POTENTIAL COMPLICATIONS \n"
"Splitting into leaflets assumes the membrane has already been centered at z=0"
" \n"
;
//clang-format on
class ToolOptions : public opts::OptionsPackage {
public:
  ToolOptions() {}
  // clang-format off
  void addGeneric(po::options_description& o) { 
    o.add_options()
      ("symmetry", po::int(&sym->default_value(6), "Symmetry of the lattice"))
      ("cutoff", po::double(&cutoff->default_value(10.0), "Cutoff distance for neighbors"))
    ;
  }
  /*
  // clang-format on
  string print() const {
    ostringstream oss;
    oss << boost::format("timeseries=%s,min_dist=%s,max_dist=%s,by_molecule=%b,"
                         "by_fragment=%b,num_bins=%d") %
               timeseries % min_dist % max_dist % by_molecule % by_fragment %
               num_bins;
    return (oss.str());
  }
  bool postConditions(po::variables_map &map) {
    if(by_molecule && by_fragment){
      cerr << "ERROR: --by-molecule and --by-fragment flags are mutually exclusive.\n";
      return(false);
    } else
      return(true);
    
  }
  */
  int sym;
  double cutoff;
};

int main(int argc, char *argv[]) {
  string header = invocationHeader(argc, argv);

  opts::BasicOptions *bopts = new opts::BasicOptions(fullHelpMessage);
  opts::BasicSelection *sopts = new opts::BasicSelection("all");
  opts::MultiTrajOptions *mtopts = new opts::MultiTrajOptions;
  ToolOptions *topts = new ToolOptions;

  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(topts);
  if (!options.parse(argc, argv))
    exit(-1);

  // Get the molecules to work with. First split by molecule, then apply the 
  // selection string 
  vector<AtomicGroup> molecules, lipids;
  molecules = mtopts->model.splitByMolecule();

  for (auto molecule : molecules) {
    AtomicGroup mol = selectAtoms(molecule, sopts->selection);
    if (mol.size() > 0) {
      lipids.push_back(mol);
    }
  }


  while (mtopts->trajectory->readFrame()) {
    mtopts->trajectory->updateGroupCoords(mtopts->model);

    vector<GCoord> upper, lower;
    // Split the lipids into upper and lower leaflets
    // Have to do this every time, because cholesterols can move between leaflets, 
    // especially in coarse-grained simulations
    for (auto lipid : lipids) {
      GCoord cen = lipid.centroid();
      // We zero out the z coordinate so we can use regular distance calculations
      if (cen.z > 0)
        cen.z()=0.0;
        upper.push_back(cen);
      else
        cen.z()=0.0;
        lower.push_back(cen);
    }
    box = mtopts->model.periodicBox();

    for (auto leaflet : {upper, lower}) {

      for (uint i = 0; i < leaflet.size()-1; i++) {
        for (uint j = i+1; j < leaflet.size(); j++) {
          double dist = leaflet[i].distance2(leaflet[j], box)
          x
          if (dist < topts->cutoff) {
            // Add the neighbor to the list of neighbors
          }
        }
        /*
        double sum = 0.0;
        for (auto neighbor : neighbors) {
          double dx = neighbor.x - lipid.x;
          double dy = neighbor.y - lipid.y;
          sum += cos(6.0 * atan2(dy, dx));
        }
        hex += 1.0 - (1.0 / 6.0) * sum;
        */
      }
    }

  }
}