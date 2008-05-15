#include <iostream>
#include <vector>
#include <stdexcept>

#include "Atom.hpp"

#include "Kernel.hpp"
#include "Tokenizer.hpp"
#include "Parser.hpp"

using namespace std;
using namespace loos;

int main(int argc, char *argv[]) {
  loos::Kernel kern;

  cout << "lex: '" << argv[1] << "'\n";
  parser::parse(argv[1], kern);

  cout << kern;

  pAtom pa(new Atom(1, "CA", GCoord()));

  pa->resid(1);
  pa->resname("HIS");
  pa->segid("PROT");

  kern.execute(pa);
  if (kern.stack().size() != 1)
    throw(runtime_error("Execution error"));

  loos::Value results = kern.stack().pop();
  cout << results << endl;

  
  

}
