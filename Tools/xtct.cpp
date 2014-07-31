#include <loos.hpp>

using namespace std;
using namespace loos;


#include <xtcwriter.hpp>

void copyTraj(AtomicGroup& model, pTraj& in, TrajectoryWriter& out) {
  while (in->readFrame()) {
    in->updateGroupCoords(model);
    out.writeFrame(model);
  }
}



int main(int argc, char* argv[]) {

  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  
  XTCWriter xtc("foo.xtc");
  
  copyTraj(model, traj, xtc);
}
