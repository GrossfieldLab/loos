

#if !defined(OPTIONS_FRAMEWORK_HPP)
#define OPTIONS_FRAMEWORK_HPP

#include <loos.hpp>
#include <boost/program_options.hpp>


namespace loos {
  namespace DensityTools {
    namespace OptionsFramework {

      namespace po = boost::program_options;


      class OptionsPackage {
      public:
        virtual ~OptionsPackage() { }

        virtual void addGeneric(po::options_description& opts) { }
        virtual void addHidden(po::options_description& opts) { }
        virtual void addPositional(po::positional_options_description& opts) { }
        virtual string print() const { return(""); }

        // true = problem with options, false = ok
        virtual bool check(po::variables_map& map) { return(false); }
  
        virtual bool postConditions(po::variables_map& map) { return(true); }

        virtual string help() const { return(""); }
      };


      // -------------------------------------------------

      class BasicOptions : public OptionsPackage {
      public:
        BasicOptions() : verbosity(0) { }

        void addGeneric(po::options_description& opts);
        string print() const;

        int verbosity;
      };

      // -------------------------------------------------

      class OutputPrefixOptions : public OptionsPackage {
      public:
        void addGeneric(po::options_description& opts);
        string print() const;

        string prefix;
      };
      
      // -------------------------------------------------

      class BasicSelectionOptions : public OptionsPackage {
      public:
        BasicSelectionOptions() : selection("all") { }

        void addGeneric(po::options_description& opts);
        string print() const;

        string selection;
      };

      // -------------------------------------------------

      class BasicTrajectoryOptions : public OptionsPackage {
      public:
        BasicTrajectoryOptions() : skip(0), frame_index_spec("") { }

        void addGeneric(po::options_description& opts);
        void addHidden(po::options_description& opts);

        void addPositional(po::positional_options_description& pos);

        bool check(po::variables_map& map);

        bool postConditions(po::variables_map& map);

        string help() const;
        string print() const;

        string frame_index_spec;
        string model_name, traj_name;
        unsigned int skip;
      };



      // ----------------------------------------------------------------------

      typedef vector<OptionsPackage *> vOpts;


      class AggregateOptions {
      public:
        AggregateOptions& addOptions(OptionsPackage* pack);

        bool parseOptions(int argc, char *argv[]);
        string print() const;

      private:
        vector<OptionsPackage *> options;
      };


    };
  };
};


#endif
