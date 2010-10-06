/*
  fidpick
  (c) 2010 Tod D. Romo, Grossfield Lab, URMC
  
  Picks fiducial structures for a structural histogram, a la Lyman &
  Zuckerman, Biophys J (2006) 91:164-172


  Usage- fidpick model trajectory selection output-name cutoff [seed]

*/



#include <loos.hpp>

using namespace std;
using namespace loos;

typedef vector<AtomicGroup>      vGroup;



vector<uint> findFreeFrames(const vector<int>& map) {
  vector<uint> indices;

  for (uint i=0; i<map.size(); ++i)
    if (map[i] < 0)
      indices.push_back(i);

  return(indices);
}



int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc < 7 || argc > 8) {
    cerr << "Usage - fidpick model trajectory range|all selection output-name cutoff [seed]\n";
    exit(-1);
  }

  int opti = 1;
  AtomicGroup model = createSystem(argv[opti++]);
  model.clearBonds();

  pTraj traj = createTrajectory(argv[opti++], model);
  string range(argv[opti++]);
  string selection(argv[opti++]);
  AtomicGroup subset = selectAtoms(model, selection);
  string outname(argv[opti++]);
  double cutoff = strtod(argv[opti++], 0);

  long seed;
  if (opti < argc) {
    seed = strtol(argv[opti++], 0, 10);
    rng_singleton().seed(static_cast<unsigned int>(seed));   // This addresses a bug in earlier BOOSTs
  } else
    randomSeedRNG();

  vector<uint> frames;
  if (range == "all")
    for (uint i=0; i<traj->nframes(); ++i)
      frames.push_back(i);
  else
    frames = parseRangeList<uint>(range);


  boost::uniform_real<> rmap;
  boost::variate_generator< base_generator_type&, boost::uniform_real<> > rng(rng_singleton(), rmap);

  vGroup fiducials;
  vector<int> assignments(frames.size(), -1);
  
  cout << "Frames picked:\n";

  vector<uint> possible_frames = findFreeFrames(assignments);
  while (! possible_frames.empty()) {
    uint pick = possible_frames[static_cast<uint>(floor(possible_frames.size() * rng()))];

    traj->readFrame(frames[pick]);
    traj->updateGroupCoords(model);

    AtomicGroup fiducial = subset.copy();
    fiducial.centerAtOrigin();
    uint myid = fiducials.size();
    if (assignments[pick] >= 0) {
      cerr << "INTERNAL ERROR - " << pick << " pick was already assigned to " << assignments[pick] << endl;
      exit(-99);
    }

    fiducials.push_back(fiducial);
    assignments[pick] = myid;
    
    uint cluster_size = 0;
    for (uint i = 0; i<assignments.size(); ++i) {
      if (assignments[i] >= 0 || i == pick)
        continue;
      traj->readFrame(frames[i]);
      traj->updateGroupCoords(model);
      subset.centerAtOrigin();
      subset.alignOnto(fiducial);
      double d = subset.rmsd(fiducial);
      if (d < cutoff) {
        assignments[i] = myid;
        ++cluster_size;
      }
    }

    cout << "\t" << frames[pick] << "\t" << cluster_size << endl;

    possible_frames = findFreeFrames(assignments);
  }

  cerr << "Done!\nWrote " << fiducials.size() << " fiducials to " << outname << endl;

  DCDWriter dcd(outname + ".dcd", fiducials, hdr);

  PDB pdb = PDB::fromAtomicGroup(fiducials[0]);
  pdb.renumber();
  pdb.remarks().add(hdr);
  string fname = outname + ".pdb";
  ofstream ofs(fname.c_str());
  ofs << pdb;

}
