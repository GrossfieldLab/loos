/*
  vsa_fit

  (c) 2010 Tod D. Romo, Grossfield Lab, URMC


  Fits a basic VSA to a set of PCA results
*/



#include <loos.hpp>
#include <Simplex.hpp>
#include "fitter.hpp"
#include "anm-lib.hpp"

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
vector<ANM*> anms;
SpringFunction* spring;

int verbosity;
bool mass_flag;

vector<double> initial_seeds, initial_bound_lengths, initial_unbound_lengths;
vector<double> bound_seeds, unbound_seeds;

FitAggregator* parseOptions(int argc, char *argv[]) {

  string hdr = invocationHeader(argc, argv);

  FitAggregator* uberfit = new FitAggregator;

  try {
    string config_file;

    vString tags;
    vString models;
    vString pcas;
    vString subs;
    //    vString envs;

    vector<string> spring_name;
    vector<double> seed_scale;

    po::options_description visible("Allowed options", 120);
    visible.add_options()
      ("help", "Produce this help message")
      ("verbosity,v", po::value<int>(&verbosity)->default_value(0), "Verbosity level")
      ("mass,m", po::value<bool>(&mass_flag)->default_value(false), "Enable use of mass in VSA")
      ("config,C", po::value<string>(&config_file), "Options config file");


    po::options_description optimization("Optimization Settings");
    optimization.add_options()
      ("spring", po::value< vector<string> >(&spring_name), "Spring function to use")
      ("length", po::value< vector<double> >(&seed_scale), "Scale for seed lengths")
      ("seeds", po::value< vector<double> >(&initial_seeds), "Seed values");


    po::options_description system("System Description");
    system.add_options()
      ("tag", po::value<vString>(&tags)->composing(), "Name to associate with system")
      ("model", po::value<vString>(&models)->composing(), "Model coordinates")
      ("sub", po::value<vString>(&subs)->composing(), "Subsystem selection")
      //      ("env", po::value<vString>(&envs)->composing(), "Environment selection")
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
      cerr << "Usage- vsa_fit_bonds [options] <systems> spring length seeds [seeds]\n";
      cerr << visible;
      cerr << system;
      cerr << optimization;
      exit(-1);
    }

    // Set up global spring function...
    
    bound_spring = springFactory(spring_name[0]);
    unbound_spring = springFactory(spring_name[1]);
    uint nargs = bound_spring->paramSize();
    nargs += unbound_spring->paramSize();
    if (initial_seeds.size() != nargs) {
      cerr << "Error- Your springs wanted " << nargs << " total seed values.\n";
      exit(-2);
    }

    
    for (uint i = 0; i <= bound_spring->paramSize(); ++i){
      initial_bound_lengths.push_back(initial_seeds[i] * seed_scale[0]);
      bound_seeds.pushback(initial_seeds[i]);
    }
    for (uint i = bound_spring->paramSize(); i <= nargs; ++i)[
      initial_unbound_lengths.push_back(initial_seeds[i] * seed_scale[1]);
      unbound_seeds.pushback(initial_seeds[i]);
    }


    for (uint i=0; i<tags.size(); ++i) {
      AtomicGroup model = createSystem(models[i]);
      if (mass_flag)
        massFromOccupancy(model);
      AtomicGroup subset = selectAtoms(model, subs[i]);
      //AtomicGroup environment = selectAtoms(model, envs[i]);
      //AtomicGroup combined = subsystem + environment;
    
      DoubleMatrix s;
      readAsciiMatrix(pcas[i] + "_s.asc", s);
    
      DoubleMatrix U;
      readAsciiMatrix(pcas[i] + "_U.asc", U);
      ///////////////////////////////////////////
      
      //   Filling the connectivity map
      //   Decides which spring function to use..
      loos::Math::Matrix<int> connectivity_map(subset.size(), subset.size());
      if (subset.hasBonds()){
	for (int j = 0; j < subset.size(); ++j){
	  if (subset[j]->hasBonds()){
	    for (int k = 0; k < subset.size(); ++k) {
	      if (subset[j]->isBoundTo(subset[k]->id()))
		connectivity_map(j,k) = 1;
	      else
		connectivity_map(j,k) = 0;
	    }
	  }
	}
      }
      //////////////////////////////////////////////////////

      // Now setup blocker & springs...
      SuperBlock *blocker = new SuperBlock(unbound_spring, subset);
      BoundSuperBlock *decBlocker = new BoundSuperBlock(blocker, bound_spring, connectivity_map);
    
      ANM* anm = new ANM(decBlocker);//, subsystem.size());

      // if (mass_flag) {
      //   DoubleMatrix M = getMasses(combined);
      //   anm->setMasses(M);
      // }
  
      Fitter* fitter = new Fitter(anm, s, U);
      fitter->name(tags[i]);
      fitter->verbose(true);
      fitter->normalize(true);

      uberfit->push_back(fitter);

      fits.push_back(fitter);
      blocks.push_back(decBlocker);
      anms.push_back(anm);
    }

    // Now dump the options for logging...
    cout << "# " << hdr << endl;
    vector<string> opts = optionsValues(vm);
    copy(opts.begin(), opts.end(), ostream_iterator<string>(cout, "\n"));
    
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
  
  Simplex<double> bound_simp(bound_spring->paramSize());
  simp.tolerance(1e-4);
  Simplex<double> unbound_simp(bound_spring->paramSize());
  simp.tolerance(1e-4);
  
  bound_simp.seedLengths(initial_bound_lengths);
  unbound_simp.seedLengths(initial_unbound_lengths);

  // Do a quick check first...
  cout << "----INITIAL----\n";
  double check = (*uberfit)(bound_seeds);
  double check = (*uberfit)(unbound_seeds);
  cout << "----INITIAL----\n";
  uberfit->resetCount();


  vector<double> bound_fit = bound_simp.optimize(bound_seeds, *uberfit);
vector<double> unbound_fit = unbound_simp.optimize(unbound_seeds, *uberfit);

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
    delete anms[i];
  }
  delete unbound_spring;
  delete bound_spring;

}
