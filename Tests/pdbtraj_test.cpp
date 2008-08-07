#include <loos.hpp>
#include <boost/format.hpp>

using namespace boost;


int main(int argc, char *argv[]) {
  PDB pdb(argv[1]);
  uint start = atoi(argv[3]);
  uint stop = atoi(argv[4]);
  uint stride = atoi(argv[5]);

  PDBTraj pdbt(argv[2], start, stop, stride);
  
  cout << format("There are %u frames with %u atoms.\n") % pdbt.nframes() % pdbt.natoms();
  int i = 0;
  while (pdbt.readFrame()) {
    cout << format("Reading frame %d...\n") % i++;
    cout << "\tname = " << pdbt.currentName() << endl;
    pdbt.updateGroupCoords(pdb);
    vector<GCoord> bdd = pdb.boundingBox();
    cout << "\tCenter @ " << pdb.centroid() << " with bdd " << bdd[0] << " x " << bdd[1] << endl;
  }

  cout << format("Read in a total of %d frames.\n") % i;
}
