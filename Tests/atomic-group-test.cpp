/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <string>
#include <vector>
#include <tr1/memory>


using namespace std;
using namespace tr1;

#include <loos.hpp>
#include <Atom.hpp>
#include <AtomicGroup.hpp>


struct SelectorCA : public AtomSelector {
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
  cout << "Operator= test:\n";
  g1 = g1;
  cout << "operator== test:\n";
  AtomicGroup ggg = g1;
  cout << (ggg == g1) << endl;
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
  SelectorCA sel;

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
  vector<GCoord> bb = g1.boundingBox();
  cout << bb[0] << " x " << bb[1] << endl;

  cout << "-------------------\n";
  cout << "Radius = " << g1.radius() << endl;
  cout << "Rgyr = " << g1.radiusOfGyration() << endl;

  cout << "-------------------\n";
  cout << g1.numberOfResidues() << " Residues in group.\n";
  AtomicGroup g2 = g1.getResidue(b);
  cout << "Size = " << g2.size() << endl;
  cout << g2 << endl;

  cout << "-------------------\n";
  AtomicGroup g4 = g1.copy();
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

  cout << "-------------------\n";
  cout << "Box test:\n";
  AtomicGroup g7 = g1;
  g1.periodicBox(GCoord(13, 26, 39));
  cout << g1 << endl;

  cout << "-------------------\n";
  cout << "Box inheritance test:\n";
  cout << g3 << endl;
  cout << "Test updating of box to (25,26,29):\n";
  g1.periodicBox(25,26,29);
  cout << g1.periodicBox() << endl;
  cout << g3.periodicBox() << endl;
  cout << "Test updating derived box (pre-box) to (7,8,9):\n";
  g7.periodicBox(7,8,9);
  cout << g1.periodicBox() << endl;
  cout << "Testing copy (should be (7,8,9):\n";
  AtomicGroup g6 = g1.copy();
  g1.periodicBox(1,2,3);
  cout << g6.periodicBox() << endl;

}
