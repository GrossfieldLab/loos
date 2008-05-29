#include <pdb.hpp>
#include <Selectors.hpp>

int main(int argc, char *argv[]) {
  PDB one(argv[1]);
  PDB two(argv[2]);

  CAlphaSelector casel;
  AtomicGroup oneca = one.select(casel);
  AtomicGroup twoca = two.select(casel);

  GMatrix M = oneca.alignOnto(twoca);

  GMatrix W = oneca.xform().current();
  one.xform().load(W);
  one.applyTransform();

  cout << one;
}
