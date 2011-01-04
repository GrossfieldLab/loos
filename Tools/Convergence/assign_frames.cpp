#include <loos.hpp>

#include "fid-lib.hpp"

using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  string hdr = invocationHeader(argc, argv);

  if (argc != 6) {
    cerr << "Usage - " << argv[0] << " model trajectory range selection fiducials\n";
    exit(-1);
  }
  
  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  string range(argv[k++]);
  string selection(argv[k++]);
  
  AtomicGroup subset = selectAtoms(model, selection);
  AtomicGroup ref_model = subset.copy();
  ref_model.renumber();
  pTraj fiducials = createTrajectory(argv[k++], ref_model);

  vecUint frames;
  if (range == "all")
    for (uint i=0; i<traj->nframes(); ++i)
      frames.push_back(i);
  else
    frames = parseRangeList<uint>(range);

  vecGroup refs;
  cerr << "Reading fiducials...\n";
  readTrajectory(refs, ref_model, fiducials);
  cerr << "Read in " << refs.size() << " fiducials.\n";
  cerr << "Assigning...\n";
  vecUint assigned = assignStructures(subset, traj, frames, refs);
  cout << "# " << hdr << endl;
  copy(assigned.begin(), assigned.end(), ostream_iterator<uint>(cout, "\n"));

}
