// (c) 2012 Tod D. Romo, Grossfield Lab, URMC

#include <loos.hpp>

using namespace std;
using namespace loos;


int main(int argc, char *argv[]) {
  
  AtomicGroup model = createSystem(argv[1]);
  pTraj traj = createTrajectory(argv[2], model);

  vector<AtomicGroup> frames;
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);
    AtomicGroup frame = model.copy();
    frames.push_back(frame);
  }

  bool b = traj->rewind();
  if (!b) {
    cerr << "Rewind reported error.\n";
    exit(-1);
  }
  uint k = 0;
  while (traj->readFrame()) {
    traj->updateGroupCoords(model);
    for (uint i=0; i<model.size(); ++i)
      if (model[i]->coords() != frames[k][i]->coords()) {
        cerr << "Error- mismatch for frame " << k << " at atom " << i << endl;
        cerr << "Expected:\n" << *(frames[k][i]) << endl;
        cerr << "Got:\n" << *(model[i]) << endl;
        exit(-1);
      }
    ++k;
  }
}
