#include <loos.hpp>
#include <boost/format.hpp>

using namespace boost;

int main(int argc, char *argv[]) {
  TinkerArc tarc(argv[1]);
  
  cout << format("There are %u frames with %u atoms.\n") % tarc.nframes() % tarc.natoms();
  int i = 0;
  while (tarc.readFrame()) {
    cout << format("Reading frame %d...\n") % i++;
    TinkerXYZ bell = tarc.currentFrame();
    vector<GCoord> bdd = bell.boundingBox();
    cout << "\tCenter @ " << bell.centroid() << " with bdd " << bdd[0] << " x " << bdd[1] << endl;
  }

  cout << format("Read in a total of %d frames.\n") % i;

  tarc.readFrame(tarc.nframes() - 5);
  i = 0;
  while (tarc.readFrame()) {
    cout << format("Reading frame %d...\n") % i++;
    TinkerXYZ bell = tarc.currentFrame();
    vector<GCoord> bdd = bell.boundingBox();
    cout << "\tCenter @ " << bell.centroid() << " with bdd " << bdd[0] << " x " << bdd[1] << endl;
  }

  cout << format("Read in a total of %d frames.\n") % i;
  
}
