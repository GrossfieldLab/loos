#include <iostream>
#include <vector>
#include <stdexcept>

#include "Atom.hpp"

#include "Kernel.hpp"
#include "Parser.hpp"

using namespace std;
using namespace loos;

int main(int argc, char *argv[]) {
  Parser p(argv[1]);
  cout << p.kernel() << endl;
}
