//
// Quick test of fitting
// (c) 2010 Tod D. Romo, Grossfield Lab, URMC
//



#include <loos.hpp>
#include <Simplex.hpp>
#include "fitter.hpp"
#include "vsa-lib.hpp"


using namespace loos;
using namespace std;


int main(int argc, char *argv[]) {
  if (argc == 1) {
    cout << "Usage- quickie tag model subsystem environment eigvals eigvecs [tag model sub env eigvals eigvecs ...]\n";
    exit(0);
  }

  string hdr = invocationHeader(argc, argv);


  int k=1;

  FitAggregator uberfit;

  // Track allocations for cleanup...
  vector<ENMFitter*> fits;
  vector<SuperBlock*> blocks;
  SpringFunction *spring = springFactory(argv[k++]);
  uint nargs = spring->paramSize();
  vector<double> seeds;
  cout << "Expecting " << nargs << " seeds.\n";
  for (uint i=0; i<nargs; ++i)
    seeds.push_back(strtod(argv[k++], 0));

  while (k < argc) {
    string tag(argv[k++]);
    AtomicGroup model = createSystem(argv[k++]);
    AtomicGroup subsystem = selectAtoms(model, argv[k++]);
    AtomicGroup environment = selectAtoms(model, argv[k++]);
    AtomicGroup combined = subsystem + environment;
    
    DoubleMatrix s;
    readAsciiMatrix(argv[k++], s);
    
    DoubleMatrix U;
    readAsciiMatrix(argv[k++], U);
    
    // Now setup blocker & springs...
    SuperBlock *blocker = new SuperBlock(spring, combined);
    
    VSA* vsa = new VSA(blocker, subsystem.size());
  
    ENMFitter* fitter = new ENMFitter(vsa, s, U);
    fitter->name(tag);
    fitter->verbose(true);
    fitter->normalize(true);

    uberfit.push_back(fitter);

    fits.push_back(fitter);
    blocks.push_back(blocker);
  }


  
  Simplex<double> simp(nargs);
  simp.tolerance(1e-4);
  
  vector<double> lengths;
  for (vector<double>::iterator vi = seeds.begin(); vi != seeds.end(); ++vi)
    lengths.push_back(*vi/2.0);

  simp.seedLengths(lengths);

  // Do a quick check first...
  cout << "----INITIAL----\n";
  double check = uberfit(seeds);
  cout << "----INITIAL----\n";
  uberfit.resetCount();


  vector<double> fit = simp.optimize(seeds, uberfit);

  cout << "----FINAL----\n";
  cout << simp.finalValue() << "\t= ";
  copy(fit.begin(), fit.end(), ostream_iterator<double>(cout, "\t"));
  cout << endl;
  uberfit.resetCount();
  check = uberfit(fit);
  cout << "----FINAL----\n";
  

  // Cleanup (make valgrind happy)
  for (uint i=0; i<fits.size(); ++i) {
    delete fits[i];
    delete blocks[i];
  }
  delete spring;

}
