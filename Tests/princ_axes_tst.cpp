#include <iostream>
#include <string>
#include <vector>
#include <tr1/memory>


using namespace std;
using namespace tr1;

#include <loos.hpp>
#include <Atom.hpp>
#include <AtomicGroup.hpp>

int main() {
  pAtom a0(new Atom(1, "CA", GCoord(0, 0, 0)));
  pAtom a1(new Atom(2, "CA", GCoord(8, 2, 0)));
  pAtom a2(new Atom(3, "CA", GCoord(-12, -1, 0)));
  pAtom a3(new Atom(4, "CA", GCoord(11, 1, 0)));
  pAtom a4(new Atom(5, "CA", GCoord(3, -2, 0)));
  pAtom a5(new Atom(6, "CA", GCoord(22, 2, 0)));
  pAtom a6(new Atom(7, "CA", GCoord(-16, -1, 0)));
  pAtom a7(new Atom(8, "CA", GCoord(8, 1, 0)));

  AtomicGroup g;

  g.append(a0);
  g.append(a1);
  g.append(a2);
  g.append(a3);
  g.append(a4);
  g.append(a5);
  g.append(a6);
  g.append(a7);

  cout << g << endl;

  vector<GCoord> foo = g.principalAxes();
  cout << foo[0] << endl;
  cout << foo[1] << endl;
  cout << foo[2] << endl;
  cout << foo[3] << endl;

}
