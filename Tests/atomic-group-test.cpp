#include <iostream>
#include <string>
#include <vector>
#include <tr1/memory>


using namespace std;
using namespace tr1;

#include <loos.hpp>
#include <Atom.hpp>
#include <AtomicGroup.hpp>


struct CASelector : public AtomSelector {
  bool operator()(const pAtom& atom) const {
    return(atom->name() == "CA");
  }
};
  


int main() {
  pAtom a(new Atom(1, "CA", GCoord(1,2,3)));
  a->resid(1);
  pAtom b(new Atom(2, "CB", GCoord(4,5,6)));
  b->resid(1);
  pAtom c(new Atom(3, "CG", GCoord(-1, -2, -3)));
  c->resid(1);

  pAtom d(new Atom(4, "CA", GCoord(7,3,4)));
  d->resid(2);
  pAtom e(new Atom(5, "CG", GCoord(1,7,2)));
  e->resid(2);
  pAtom f(new Atom(4, "C", GCoord(1,7,2)));
  f->resid(2);

  a->addBond(b);
  a->addBond(c);
  cout << "Atom a:\n" << *a << endl;

  AtomicGroup g1;


  g1.append(a);
  g1.append(b);
  g1.append(c);
  g1.append(d);
  g1.append(e);

  cout << g1 << endl;
  cout << "-------------------\n";
  cout << "Operator[] test:\n";
  pAtom tmpatm = g1[3];
  g1[3] = f;
  cout << g1 << endl;
  cout << "=====\n";
  g1[3] = tmpatm;
  cout << g1 << endl;

  cout << "-------------------\n";
  cout << "Operator+ tests:\n";
  AtomicGroup gg;
  gg.append(a);
  gg += b+c;
  cout << gg << endl;
  cout << "======\n";
  gg += g1;
  cout << gg << endl;
  cout << "=====\n";
  gg = a+b+c+d+e+f;
  cout << gg << endl;

  cout << "-------------------\n";
  cout << "CA selection:\n";
  CASelector sel;

  AtomicGroup s = g1.select(sel);
  cout << s << endl;



  cout << "-------------------\n";
  cout << "Iterator test:\n";
  AtomicGroup::Iterator ii(g1);
  pAtom pa;
  while (pa = ii())
    cout << *pa << endl;

  cout << "-------------------\n";
  cout << "Bounds test...\n";
  AtomicGroup::BoundingBox bb = g1.boundingBox();
  int k;
  for (k=0; k<3; k++)
    cout << "min[" << k << "] = " << bb.min[k] << ", max[" << k << "] = " << bb.max[k] << endl;

  cout << "-------------------\n";
  cout << "Radius = " << g1.radius() << endl;
  cout << "Rgyr = " << g1.radiusOfGyration() << endl;

  cout << "-------------------\n";
  cout << g1.numberOfResidues() << " Residues in group.\n";
  AtomicGroup g2 = g1.getResidue(b);
  cout << "Size = " << g2.size() << endl;
  cout << g2 << endl;

  cout << "-------------------\n";
  AtomicGroup g4 = *(g1.clone());
  cout << "Clone & sharing test:\n";
  g4[0]->resid(999);
  cout << g4 << endl;
  cout << "-" << endl;
  cout << g1 << endl;

  cout << "---\n";
  AtomicGroup g5 = g4;
  g5[0]->resid(111);
  cout << g4 << endl;

  cout << "-------------------\n";
  cout << "subset(1,2):\n";
  AtomicGroup g3 = g1.subset(1,2);
  cout << g3 << endl;
  cout << "\nsubset(-2):\n";
  g3 = g1.subset(-2);
  cout << g3 << endl;
  cout << "\nsubset(-2,4):\n";
  g3 = g1.subset(-2,2);
  cout << g3 << endl;

  cout << "-------------------\n";
  cout << "Excise(1,2):\n";
  g3 = g1.excise(1,2);
  cout << g3 << endl << g1 << endl;

}
