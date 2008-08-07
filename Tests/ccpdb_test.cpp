#include <loos.hpp>
#include <boost/format.hpp>

using namespace boost;

int main(int argc, char *argv[]) {
  PDB pdb(argv[1]);
  CCPDB ccpdb(argv[2]);
  
  cout << format("There are %u frames with %u atoms.\n") % ccpdb.nframes() % ccpdb.natoms();
  int i = 0;
  while (ccpdb.readFrame()) {
    cout << format("Reading frame %d...\n") % i++;
    ccpdb.updateGroupCoords(pdb);
    vector<GCoord> bdd = pdb.boundingBox();
    cout << "\tCenter @ " << pdb.centroid() << " with bdd " << bdd[0] << " x " << bdd[1] << endl;
  }

  cout << format("Read in a total of %d frames.\n") % i;
}
