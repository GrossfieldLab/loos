

#if !defined(OPTIONS_FRAMEWORK_HPP)
#define OPTIONS_FRAMEWORK_HPP

#include <loos.hpp>
#include <boost/program_options.hpp>


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
      virtual bool postConditions(po::variables_map& map) {
        return(true);
      }

      //! Returns a slice of the example command-line in the help output
      virtual std::string help() const { return(""); }
    };


    // -------------------------------------------------

    //! Options that SHOULD be common to all tools
    class BasicOptions : public OptionsPackage {
    public:
      BasicOptions() : verbosity(0) { }
      BasicOptions(const int i) : verbosity(i) { }
      BasicOptions(const std::string& s) : verbosity(0), full_help(s) { }
      BasicOptions(const int i, const std::string& s) : verbosity(i), full_help(s) { }

      void addGeneric(po::options_description& opts);
      bool check(po::variables_map& map);

      std::string print() const;

      void setFullHelp(const std::string& s);

      int verbosity;
      std::string full_help;
    };

    // -------------------------------------------------

    //! Options related to specifying an output prefix
    class OutputPrefix : public OptionsPackage {
    public:
      OutputPrefix() : prefix("output") { }
      OutputPrefix(const std::string& s) : prefix(s) { }

      void addGeneric(po::options_description& opts);
      std::string print() const;

      std::string prefix;
    };
      
    // -------------------------------------------------

    //! Provide a single selection
    class BasicSelection : public OptionsPackage {
    public:
      BasicSelection() : selection("all") { }
      BasicSelection(const std::string& sel) : selection(sel) { }

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
    class ModelWithCoords : public OptionsPackage {
    public:
      ModelWithCoords() : coords_name("") { }

      void addGeneric(po::options_description& opts);
      void addHidden(po::options_description& opts);
      void addPositional(po::positional_options_description& pos);

      bool check(po::variables_map& map);

      bool postConditions(po::variables_map& map);

      std::string help() const;
      std::string print() const;

      std::string model_name, coords_name;

      AtomicGroup model;

    };



    // -------------------------------------------------

    //! Basic trajectory options
    /**
     * Adds a model and trajectory argument to the command line, and
     * provides --skip (-k) option for skipping the first n-frames.
     * 
     * The contained trajectory object will already be skipped to the
     * correct frame by postConditions().
     **/
    class BasicTrajectory : public OptionsPackage {
    public:
      BasicTrajectory() : skip(0) { }

      void addGeneric(po::options_description& opts);
      void addHidden(po::options_description& opts);

      void addPositional(po::positional_options_description& pos);

      bool check(po::variables_map& map);

      bool postConditions(po::variables_map& map);

      std::string help() const;
      std::string print() const;

      std::vector<uint> frameList() const;

      unsigned int skip;
      std::string model_name, traj_name;

      AtomicGroup model;
      pTraj trajectory;
    };



    // -------------------------------------------------

    //! Trajectory with range or skip
    /**
     * Adds a model and trajectory argument to the command line, and
     * provides --skip (-k) and --range (-r) options for specifying
     * which frames of the trajectory to operate over.
     **/
    class TrajectoryWithFrameIndices : public OptionsPackage {
    public:
      TrajectoryWithFrameIndices() : skip(0), frame_index_spec("") { }

      void addGeneric(po::options_description& opts);
      void addHidden(po::options_description& opts);

      void addPositional(po::positional_options_description& pos);

      bool check(po::variables_map& map);

      bool postConditions(po::variables_map& map);

      std::string help() const;
      std::string print() const;

      std::vector<uint> frameList() const;

      unsigned int skip;
      std::string frame_index_spec;
      std::string model_name, traj_name;

      AtomicGroup model;
      pTraj trajectory;
    };


    // -------------------------------------------------

    //! Provides simple way to add command-line arguments (required options)
    /**
     * This class handles required command-line options (also known as
     * command line arguments).  Each argument is defined by a string
     * tag and a description and is parsed from the command line as a
     * string.  Arguments are added via the addOption() method and the
     * values are retrieved using value().
     *
     * Since these are required options, the class will automatically
     * generate a parsing error if any argument is unset.
    **/
    class RequiredArguments : public OptionsPackage {
      typedef std::pair<std::string, std::string>    StringPair;
    public:
      RequiredArguments() : vargs_set(false) { }

      void addHidden(po::options_description& o);
      void addPositional(po::positional_options_description& pos);
      bool check(po::variables_map& map);
      bool postConditions(po::variables_map& map);
      std::string help() const;
      std::string print() const;

      //! Add a required argument given a name (tag) and a description (currently unused)
      void addArgument(const std::string& name, const std::string& description);

      //! Add a required argument that can be an arbitrary number of items
      /**
       * This argument will always appear at the end of the command
       * line, after all other required arguments.  It also means that
       * the RequiredOptions object should be the last one added to
       * the AggregateOptions object, otherwise any subsequence
       * positional arguments will be missed...
       */
      void addVariableArguments(const std::string& name, const std::string& description);

      //! Retrieve the value for an argument
      std::string value(const std::string& s) const;

      //! Retreive the variable-number argument
      std::vector<std::string> variableValues(const std::string& s) const;

    private:
      bool vargs_set;
      std::vector<StringPair> arguments;
      StringPair variable_arguments;
      po::variables_map held_map;
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
      std::vector<std::string> print() const;

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



    std::string stringVectorAsStringWithCommas(const std::vector<std::string>& v);

  };
};


#endif
