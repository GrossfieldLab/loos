#include <loos.hpp>
#include <boost/format.hpp>

using namespace boost;

int main(int argc, char *argv[]) {
  CCPDB ccpdb(argv[1]);
  
  cout << format("There are %u frames with %u atoms.\n") % ccpdb.nframes() % ccpdb.natoms();
  int i = 0;
  while (ccpdb.readFrame()) {
    cout << format("Reading frame %d...\n") % i++;
    PDB pdb = ccpdb.currentFrame();
    vector<GCoord> bdd = pdb.boundingBox();
    cout << "\tCenter @ " << pdb.centroid() << " with bdd " << bdd[0] << " x " << bdd[1] << endl;
  }

  cout << format("Read in a total of %d frames.\n") % i;
  cout << "--------------------------------------\n";

  uint f = ccpdb.nframes() - 6;
  cout << "Reading frame " << f << endl;
  ccpdb.readFrame(f);
  i=0;
  while (ccpdb.readFrame()) {
    cout << format("Reading frame %d...\n") % i++;
    PDB pdb = ccpdb.currentFrame();
    vector<GCoord> bdd = pdb.boundingBox();
    cout << "\tCenter @ " << pdb.centroid() << " with bdd " << bdd[0] << " x " << bdd[1] << endl;
  }

  cout << format("Read in a total of %d frames.\n") % i;

}
