#include <loos.hpp>
#include <pdb.hpp>


#define ERROR_CHECK


int main(int argc, char *argv[]) {
  long n = atol(argv[1]);
  PDB pdb(argv[2]);

  int a = pdb.minId();
  int b = pdb.maxId();
  int d = b-a;

  int *indices = new int[n];
  int i;

  cerr << "Generating indices...\n";
  srand48(time(0));
  for (i=0; i<n; i++)
    indices[i] = (int)(drand48() * d) + a;

  cerr << "Searching...";
  pAtom pa;

  for (i=0; i<n; i++) {
    pa = pdb.findById(indices[i]);
#if defined(ERROR_CHECK)
    if (pa->id() != indices[i]) {
      cerr << "***ERROR*** at index " << i << "(" << indices[i] << ")\n";
      exit(-20);
    } else if (pa == 0)
      exit(-10);
#endif
  }

  cerr << "done\n";
}
