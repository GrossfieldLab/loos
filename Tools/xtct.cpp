#include <loos.hpp>

using namespace std;
using namespace loos;


#include <xtcwriter.hpp>

int main(int argc, char* argv[]) {

  int k = 1;
  AtomicGroup model = createSystem(argv[k++]);
  pTraj traj = createTrajectory(argv[k++], model);
  
  XTCWriter xtc("foo.xtc");
  
  for (uint i=0; i<10; ++i) {
    traj->readFrame();
    traj->updateGroupCoords(model);
    xtc.writeFrame(model);
  }

}
