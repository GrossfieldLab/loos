

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
        virtual std::string print() const { return(""); }

        // true = problem with options, false = ok
        virtual bool check(po::variables_map& map) { return(false); }
  
        virtual bool postConditions(po::variables_map& map) { return(true); }

        virtual std::string help() const { return(""); }
      };


      // -------------------------------------------------

      class BasicOptions : public OptionsPackage {
      public:
        BasicOptions() : verbosity(0) { }

        void addGeneric(po::options_description& opts);
        std::string print() const;

        int verbosity;
      };

      // -------------------------------------------------

      class OutputPrefixOptions : public OptionsPackage {
      public:
        void addGeneric(po::options_description& opts);
        std::string print() const;

        std::string prefix;
      };
      
      // -------------------------------------------------

      class BasicSelectionOptions : public OptionsPackage {
      public:
        BasicSelectionOptions() : selection("all") { }

        void addGeneric(po::options_description& opts);
        std::string print() const;

        std::string selection;
      };

      // -------------------------------------------------

      class ModelWithCoordsOptions : public OptionsPackage {
      public:
        ModelWithCoordsOptions() : coords_name("") { }

        void addGeneric(po::options_description& opts);
        void addHidden(po::options_description& opts);
        void addPositional(po::positional_options_description& pos);

        bool check(po::variables_map& map);

        std::string help() const;
        std::string print() const;
        
        std::string model_name, coords_name;
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

        std::string help() const;
        std::string print() const;

        unsigned int skip;
        std::string frame_index_spec;
        std::string model_name, traj_name;
      };



      // ----------------------------------------------------------------------

      typedef std::vector<OptionsPackage *> vOpts;


      class AggregateOptions {
      public:
        AggregateOptions() : program_name("unknown_tool"),
                             generic("Allowed Options"),
                             hidden("Hidden Options")
        { }

        AggregateOptions& add(OptionsPackage* pack);

        bool parse(int argc, char *argv[]);
        std::string print() const;
        void showHelp();

      private:
        std::string program_name;

        po::options_description generic;
        po::options_description hidden;
        po::options_description command_line;
        po::positional_options_description pos;
        po::variables_map vm;

        std::vector<OptionsPackage *> options;


        void setupOptions();
      };

      // ----------------------------------------------------------------------

      std::vector<uint> assignFrameIndices(pTraj& traj, const std::string& desc, const uint skip);
      AtomicGroup loadStructureWithCoords(const std::string model_name, const std::string optional_coords_name);


    };
  };
};


#endif
