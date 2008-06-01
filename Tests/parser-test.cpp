#include <iostream>
#include <vector>
#include <stdexcept>

#include "Atom.hpp"

#include "Kernel.hpp"
#include "Parser.hpp"

using namespace std;
using namespace loos;

void test(const string s, bool expected = false) {
  Parser p;
  try {
    cout << "Parsing '" << s << "'\n";
    p.parse(s);
  }
  catch (...) {
    if (expected)
      cout << "Expected exception caught.\n";
    else
      cout << "=====> UNEXPECTED EXCEPTION <=====\n";
    return;
  }
  if (expected)
    cout << "=====> EXPECTED EXCEPTION NOT FOUND <=====\n";
  cout << p.kernel() << endl;
}


int main(int argc, char *argv[]) {

  test("name == 'CA'");
  test("resid =~ '1\\d+'", true);
}


