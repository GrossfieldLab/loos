#include <loos.hpp>
#include <pdb.hpp>


int main(int argc, char *argv[]) {

  PDB p(argv[1]);

  cerr << "*** Read in " << p.size() << " atoms...\n";
  cout << p;

  AtomicGroup g = p;
  cout << g << endl;

  PDB q(g);

  cout << q;

  PDB r = *(q.clone());
  q[0]->resid(-999);
  cout << "One:\n" << q << endl;
  cout << "Two:\n" << r << endl;

}

