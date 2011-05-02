/*
  water-extract.cpp

  (c) 2011 Tod D. Romo, Grossfield Lab
      Department of Biochemistry
      University of Rochster School of Medicine and Dentistry


   usage:
     water-extract [options] model trajectory >output.pdb


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
  opts::BasicTrajectoryOptions* trajopts = new opts::BasicTrajectoryOptions;
  opts::BasicWaterOptions* watopts = new opts::BasicWaterOptions;

  opts::AggregateOptions options;
  options.add(basopts).add(trajopts).add(watopts);
  if (!options.parse(argc, argv))
    exit(-1);



  AtomicGroup model = createSystem(trajopts->model_name);
  pTraj traj = createTrajectory(trajopts->traj_name, model);
  vector<uint> frames = opts::assignFrameIndices(traj, trajopts->frame_index_spec, trajopts->skip);

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
