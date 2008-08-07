#include <loos.hpp>
#include <boost/format.hpp>

using namespace boost;


int main(int argc, char *argv[]) {
  uint start = atoi(argv[2]);
  uint stop = atoi(argv[3]);
  uint stride = atoi(argv[4]);

  PDBTraj pdbt(argv[1], start, stop, stride);
  
  cout << format("There are %u frames with %u atoms.\n") % pdbt.nframes() % pdbt.natoms();
  int i = 0;
  while (pdbt.readFrame()) {
    cout << format("Reading frame %d...\n") % i++;
    cout << "\tname = " << pdbt.currentName() << endl;
    PDB pdb = pdbt.currentFrame();
    vector<GCoord> bdd = pdb.boundingBox();
    cout << "\tCenter @ " << pdb.centroid() << " with bdd " << bdd[0] << " x " << bdd[1] << endl;
  }

  cout << format("Read in a total of %d frames.\n") % i;
}
