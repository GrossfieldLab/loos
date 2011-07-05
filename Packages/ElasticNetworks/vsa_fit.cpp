/*
  vsa_fit

  (c) 2010 Tod D. Romo, Grossfield Lab, URMC


  Fits a basic VSA to a set of PCA results
*/



#include <loos.hpp>
#include <Simplex.hpp>
#include "fitter.hpp"
#include "vsa-lib.hpp"

#include <boost/program_options.hpp>


using namespace loos;
using namespace std;
using namespace ENM;

namespace po = boost::program_options;


// Globals...

typedef vector<string>    vString;

// Track these for cleanup...
vector<Fitter*> fits;
vector<SuperBlock*> blocks;
vector<VSA*> vsas;
SpringFunction* spring;

int verbosity;
bool mass_flag;

vector<double> initial_seeds, initial_lengths;


FitAggregator* parseOptions(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);

  FitAggregator* uberfit = new FitAggregator;

  try {
    string config_file;

    vString tags;
    vString models;
    vString pcas;
    vString subs;
    vString envs;

    string spring_name;
    double seed_scale;

    po::options_description visible("Allowed options", 120);
    visible.add_options()
      ("help", "Produce this help message")
      ("verbosity,v", po::value<int>(&verbosity)->default_value(0), "Verbosity level")
      ("mass,m", po::value<bool>(&mass_flag)->default_value(false), "Enable use of mass in VSA")
      ("config,C", po::value<string>(&config_file), "Options config file");


    po::options_description optimization("Optimization Settings");
    optimization.add_options()
      ("spring", po::value<string>(&spring_name), "Spring function to use")
      ("length", po::value<double>(&seed_scale), "Scale for seed lengths")
      ("seeds", po::value< vector<double> >(&initial_seeds), "Seed values");


    po::options_description system("System Description");
    system.add_options()
      ("tag", po::value<vString>(&tags)->composing(), "Name to associate with system")
      ("model", po::value<vString>(&models)->composing(), "Model coordinates")
      ("sub", po::value<vString>(&subs)->composing(), "Subsystem selection")
      ("env", po::value<vString>(&envs)->composing(), "Environment selection")
      ("pca", po::value<vString>(&pcas)->composing(), "PCA file prefix");
    
    po::options_description command_line;
    command_line.add(visible).add(system).add(optimization);

    po::positional_options_description p;
    p.add("spring", 1);
    p.add("length", 1);
    p.add("seeds", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(command_line).positional(p).run(), vm);
    po::notify(vm);

    // Now handle config file...
    if (!config_file.empty()) {
      ifstream ifs(config_file.c_str());
      if (!ifs) {
        cerr << "Cannot open config file " << config_file << endl;
        exit(-1);
      }
      store(parse_config_file(ifs, command_line), vm);
      notify(vm);
    }

    if (vm.count("help") || !(vm.count("spring") && vm.count("length") && vm.count("seeds"))) {
      cerr << "Usage- vsa_fit [options] <systems> spring length seeds [seeds]\n";
      cerr << visible;
      cerr << system;
      cerr << optimization;
      exit(-1);
    }

    // Set up global spring function...
    
    spring = springFactory(spring_name);
    uint nargs = spring->paramSize();
    if (initial_seeds.size() != nargs) {
      cerr << "Error- spring wanted " << nargs << " seed values.\n";
      exit(-2);
    }

    for (vector<double>::iterator i = initial_seeds.begin(); i != initial_seeds.end(); ++i)
      initial_lengths.push_back(*i * seed_scale);



    for (uint i=0; i<tags.size(); ++i) {
      AtomicGroup model = createSystem(models[i]);
      if (mass_flag)
        massFromOccupancy(model);
      AtomicGroup subsystem = selectAtoms(model, subs[i]);
      AtomicGroup environment = selectAtoms(model, envs[i]);
      AtomicGroup combined = subsystem + environment;
    
      DoubleMatrix s;
      readAsciiMatrix(pcas[i] + "_s.asc", s);
    
      DoubleMatrix U;
      readAsciiMatrix(pcas[i] + "_U.asc", U);

      if (s.rows() < U.cols()) {
        DoubleMatrix ss(U.cols(), 1);
        for (uint i=0; i<s.rows(); ++i)
          ss[i] = s[i];

        s = ss;
      }
    
      // Now setup blocker & springs...
      SuperBlock *blocker = new SuperBlock(spring, combined);
    
      VSA* vsa = new VSA(blocker, subsystem.size());

      if (mass_flag) {
        DoubleMatrix M = getMasses(combined);
        vsa->setMasses(M);
      }
  
      Fitter* fitter = new Fitter(vsa, s, U);
      fitter->name(tags[i]);
      fitter->verbose(true);
      fitter->normalize(true);

      uberfit->push_back(fitter);

      fits.push_back(fitter);
      blocks.push_back(blocker);
      vsas.push_back(vsa);
    }

    // options dump temporarily removed
    cout << "# " << hdr << endl;
    
  }
  catch(exception& e) {
    cerr << "Error - " << e.what() << endl;
    exit(-1);
  }

  return(uberfit);
}




void showSprings(ostream& os) {
  vector<string> v = springNames();
  os << "Valid springs: ";
  for (vector<string>::iterator i = v.begin(); i != v.end(); ++i)
    os << *i << ( i == v.end() - 1 ? "" : ", ");
  os << endl;
}



int main(int argc, char *argv[]) {
  
  FitAggregator* uberfit = parseOptions(argc, argv);
  
  Simplex<double> simp(spring->paramSize());
  simp.tolerance(1e-4);
  
  simp.seedLengths(initial_lengths);

  // Do a quick check first...
  cout << "----INITIAL----\n";
  double check = (*uberfit)(initial_seeds);
  cout << "----INITIAL----\n";
  uberfit->resetCount();


  vector<double> fit = simp.optimize(initial_seeds, *uberfit);

  cout << "----FINAL----\n";
  cout << simp.finalValue() << "\t= ";
  copy(fit.begin(), fit.end(), ostream_iterator<double>(cout, "\t"));
  cout << endl;
  uberfit->resetCount();
  check = (*uberfit)(fit);
  cout << "----FINAL----\n";
  

  // Cleanup (make valgrind happy)
  for (uint i=0; i<fits.size(); ++i) {
    delete fits[i];
    delete blocks[i];
    delete vsas[i];
  }
  delete spring;

}
