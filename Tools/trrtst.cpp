#include <loos.hpp>
#include <xdr.hpp>

#include <boost/format.hpp>


using namespace std;
using namespace boost;
using namespace loos;


int main(int argc, char *argv[]) {

    
  // Now try it through trajectory interface...
  AtomicGroup model = createSystem("f.gro");
  pTraj traj = createTrajectory("f.trr", model);


  cout << "nframes = " << traj->nframes() << endl;
  cout << "natoms = " << traj->natoms() << endl;

  int n = 0;
  while (traj->readFrame()) {
    cout << format("Frame = %d\n") % n++;
    cout << format("\tBox = %s\n") % traj->periodicBox();
    traj->updateGroupCoords(model);
    for (uint i=0; i<5; ++i) {
      cout << *(model[i]) << endl;
    }
    cout << endl;
  }

  
}
