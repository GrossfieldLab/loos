#if !defined(CONVERGENCE_OPTIONS_HPP)
#define CONVERGENCE_OPTIONS_HPP



#include <boost/format.hpp>


#include <loos.hpp>


namespace loos {
  namespace OptionsFramework {

    class BasicConvergence : public OptionsPackage {
    public:
      BasicConvergence() : seed(0) { }

      void addGeneric(po::options_description& o) {
        o.add_options()
          ("seed", po::value<uint>(&seed)->default_value(seed), "Random number seed (0 = auto)");
      };

      bool postConditions(po::variables_map& vm) {
        if (seed == 0)
          seed = randomSeedRNG();
        else
          rng_singleton().seed(seed);

        return(true);
      }

      std::string print() const {
        std::ostringstream oss;
        oss << boost::format("seed=%d") % seed;
        return(oss.str());
      }

      uint seed;
    };

  };

};



#endif
