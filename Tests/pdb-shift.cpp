#include <pdb.hpp>

int main(int argc, char *argv[]) {
  PDB pdb(argv[1]);
  GCoord c = pdb.centroid();

  pdb.xform().translate(c);
  pdb.xform().translate(2.0, 3.0, 4.0);
  pdb.xform().rotate('x', 30.0);
  pdb.xform().rotate('y', 45.0);
  pdb.xform().rotate('z', 22.5);
  pdb.xform().translate(-c);
  pdb.applyTransform();

  cout << pdb;

}
