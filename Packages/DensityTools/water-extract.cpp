/*
  water-extract.cpp

   usage:
     water-extract [options] model trajectory >output.pdb
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo, Alan Grossfield
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



#include <iterator>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <loos.hpp>
#include <DensityGrid.hpp>
#include <DensityTools.hpp>
#include <DensityOptions.hpp>


using namespace std;
using namespace loos;
using namespace loos::DensityTools;



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  opts::BasicOptions* basopts = new opts::BasicOptions;
  opts::TrajectoryWithFrameIndices* tropts = new opts::TrajectoryWithFrameIndices;
  opts::BasicWater* watopts = new opts::BasicWater;

  opts::AggregateOptions options;
  options.add(basopts).add(tropts).add(watopts);
  if (!options.parse(argc, argv))
    exit(-1);



  AtomicGroup model = tropts->model;
  pTraj traj = tropts->trajectory;
  vector<uint> frames = tropts->frameList();

  AtomicGroup subset = selectAtoms(model, watopts->prot_string);
  AtomicGroup waters = selectAtoms(model, watopts->water_string);

  AtomicGroup liquid;
  uint current_id = 1;

  for (vector<uint>::iterator t = frames.begin(); t != frames.end(); ++t) {

    traj->readFrame(*t);
    traj->updateGroupCoords(model);

    vector<int> mask = watopts->filter_func->filter(waters, subset);
    for (uint j=0; j<mask.size(); ++j)
      if (mask[j]) {
        pAtom atom(new Atom(*(waters[j])));
        atom->id(current_id);
        atom->resid(current_id++);
        atom->segid("WATER");
        liquid.append(atom);
      }

  }

  PDB pdb = PDB::fromAtomicGroup(liquid);
  pdb.remarks().add(hdr);
  cout << pdb;
}
