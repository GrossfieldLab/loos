#include <iostream>
#include <vector>
#include <stdexcept>

#include "Atom.hpp"

#include "Kernel.hpp"
#include "Parser.hpp"

#include "Selectors.hpp"
#include "pdb.hpp"

using namespace std;
using namespace loos;

void test(PDB& pdb, const string s, bool expected = false) {
  Parser p;
  try {
    cout << "\n--------------------------------------\n";
    cout << "Parsing '" << s << "'\n";
    p.parse(s);
  }
  catch (...) {
    if (expected)
      cout << "Expected exception caught.\n";
    else
      cout << "===============================> UNEXPECTED EXCEPTION\n";
    return;
  }
  if (expected) {
    cout << "===============================> EXPECTED EXCEPTION NOT FOUND\n";
    return;
  }
  cout << p.kernel() << endl;
  KernelSelector ksel(p.kernel());
  AtomicGroup grp = pdb.select(ksel);
  cout << "Selected " << grp.size() << " @ " << grp.centroid();
  vector<GCoord> box = grp.boundingBox();
  cout << " & " << box[0] << " x " << box[1] << endl;
}


int main(int argc, char *argv[]) {

  PDB pdb(argv[1]);


  test(pdb, "name == 'CA'");
  test(pdb, "resid =~ '1\\d+'", true);
  test(pdb, "!(name == 'CA')");
  test(pdb, "!(name == 'CA'", true);
  test(pdb, "segid -> 'L(\\d+)' < 3");
  test(pdb, "(segid -> '(L|P)(\\d+)') <= 3");
  test(pdb, "(segid -> '(L|P)(\\d+)') <= 10 && name == 'CA'");
  test(pdb, "name =~ 'C' && (resid >= 10 && resid <= 63) && segid != 'SOLV'");
  test(pdb, "!(name =~ 'C' && (resid >= 10 && resid <= 63) && segid != 'SOLV')");
}


