

#if !defined(OPTIONS_FRAMEWORK_HPP)
#define OPTIONS_FRAMEWORK_HPP

#include <loos.hpp>
#include <boost/program_options.hpp>
#include <tr1/unordered_set>


namespace loos {

  //! Namespace for encapsulating options processing
  namespace OptionsFramework {

    namespace po = boost::program_options;


    //! Base class for options
    /**
     * Options may have a short (1-letter) equivalent.  The
     * convention is that core options (i.e. those declared in
     * OptionsFramework.cpp) should all be lower-case.  Package
     * options should be in upper-case.  Tool-specific options
     * should all be long-form, unless they are commonly used, in
     * which case it's recommended that they be upper-case.
     **/
    class OptionsPackage {
    public:
      virtual ~OptionsPackage() { }

      //! Appends generic options (those that the user can see)
      virtual void addGeneric(po::options_description& opts) { }

      //! Appends hidden options (these generally match positional)
      virtual void addHidden(po::options_description& opts) { }

      //! Appends positional options
      virtual void addPositional(po::positional_options_description& opts) { }

      //! Returns a string listing the encapsulated options suitable for logging
      virtual std::string print() const { return(""); }

      //! Validates passed options, returning true if there is a problem
      virtual bool check(po::variables_map& map) { return(false); }

      //! Post-processing of options
      virtual bool postConditions(po::variables_map& map) { return(true); }

      //! Returns a slice of the example command-line in the help output
      virtual std::string help() const { return(""); }
    };


    // -------------------------------------------------

    //! Options that SHOULD be common to all tools
    class BasicOptions : public OptionsPackage {
    public:
      BasicOptions() : verbosity(0) { }
      BasicOptions(const int i) : verbosity(i) { }

      void addGeneric(po::options_description& opts);
      std::string print() const;

      int verbosity;
    };

    // -------------------------------------------------

    //! Options related to specifying an output prefix
    class OutputPrefixOptions : public OptionsPackage {
    public:
      OutputPrefixOptions() : prefix("output") { }
      OutputPrefixOptions(const std::string& s) : prefix(s) { }

      void addGeneric(po::options_description& opts);
      std::string print() const;

      std::string prefix;
    };
      
    // -------------------------------------------------

    //! Provide a single selection
    class BasicSelectionOptions : public OptionsPackage {
    public:
      BasicSelectionOptions() : selection("all") { }
      BasicSelectionOptions(const std::string& sel) : selection(sel) { }

      void addGeneric(po::options_description& opts);
      std::string print() const;

      std::string selection;
    };

    // -------------------------------------------------

    //! Request a model with coordinates
    /**
     * Since not all formats have coordinates (i.e. PSF),
     * the coordinates can be taken from an alternate file using the
     * -c or --coordinates option.  Also adds a positional argument
     * for the model description.
     **/
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

    //! Basic trajectory options
    /**
     * Adds a model and trajectory argument to the command line, and
     * provides --skip (-k) and --range (-r) options for specifying
     * which frames of the trajectory to operate over.
     **/
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


    // -------------------------------------------------

    class RequiredOptions : public OptionsPackage {
      typedef std::tr1::unordered_map<std::string,std::string>    Hash;
    public:
      void addHidden(po::options_description& o);
      void addPositional(po::positional_options_description& pos);
      bool check(po::variables_map& map);
      bool postConditions(po::variables_map& map);
      std::string help() const;
      std::string print() const;

      void addOption(const std::string& name, const std::string& description);
      std::string value(const std::string& s);

    private:
      Hash keys;
      Hash values;
    };


    // ----------------------------------------------------------------------

    typedef std::vector<OptionsPackage *> vOpts;

    //! Combines a set of OptionsPackages
    class AggregateOptions {
    public:
      //! Name is taken from argv[0] when AggregateOptions::parse() is called
      AggregateOptions() : program_name(""),
                           generic("Allowed Options"),
                           hidden("Hidden Options")
      { }

      AggregateOptions(const std::string& name) : program_name(name),
                                                  generic("Allowed Options"),
                                                  hidden("Hidden Options")
      { }

      //! Add a pointer to an OptionsPackage that will be used for options
      /**
       * Takes a pointer to an OptionsPackage, and appends this to
       * the list of options that will be used to build up the
       * command-line.  Returns a reference to itself so that add()
       * calls can be chained.
       **/
      AggregateOptions& add(OptionsPackage* pack);

      //! Parses a command line, returning true if parsing was ok
      bool parse(int argc, char *argv[]);

      //! Prints out the values of the options in all contained OptionsPackages
      std::string print() const;

      //! Displays the help for this tool
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

    //! Generate a vector of frame indices to operate over
    /**
     * This utility function takes a string and skip and will
     * generate a vector containing the frames of the trajectory to
     * use.  If the string is nonempty, it takes priority over the
     * skip.  The string can be a comma-separated list of
     * octave/matlab-style ranges.
     **/
    std::vector<uint> assignFrameIndices(pTraj& traj, const std::string& desc, const uint skip);

    //! Load a model with an optional coordinates file
    AtomicGroup loadStructureWithCoords(const std::string model_name, const std::string optional_coords_name);


  };
};


#endif
