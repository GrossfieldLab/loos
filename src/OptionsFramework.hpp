/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#if !defined(LOOS_OPTIONS_FRAMEWORK_HPP)
#define LOOS_OPTIONS_FRAMEWORK_HPP

#include <loos_defs.hpp>
#include <boost/program_options.hpp>
#include <AtomicGroup.hpp>
#include <Trajectory.hpp>
#include <MultiTraj.hpp>
#include <sfactories.hpp>
#include <boost/algorithm/string.hpp>
#include <exceptions.hpp>
#include <Weights.hpp>
#include <UniformWeight.hpp>


namespace loos {

  //! Namespace for encapsulating options processing
  /**
   * The OptionsFramework provides a consistent set of "common" options
   * across all LOOS tools.  By mixing and matching the subclasses of
   * OptionsPackage, a tool can decide which set of "common" options
   * it will use.  Packages may also define their own package-wide set
   * of common options.  In addition, subclassing OptionsPackage
   * within a tool is an easy mechanism for providing command-line
   * options that are specific to an individual tool without having to
   * write the support-code that using boost::program_options
   * requires.
   *
   * The basic pattern for using OptionsFramework is to create a new
   * object that represents the options required by the tool.  Once
   * all objects are created, they are added to an AggregateOptions
   * object that is responsible for parsing and validating the
   * command-line as well as generating help messages and serializing
   * the options' state for inclusion in the tool's output.  For
   * example, a tool that requires a model (with coordinates) and an
   * output prefix would setup its options as follows:
   * \code
   * BasicOptions* bopts = new BasicOptions;        // Create basic option (i.e. verbosity & help)
   * OutputPrefix* popts = new OutputPrefix;        // Allow --prefix
   * ModelWithCoords* mopts = new ModelWithCoords;  // Add a model (and ensure coordinates are present)
   *
   * AggregateOptions options;
   * options.add(bopts).add(popts).add(mopts);      // Combine all options
   * if (!options.parse(argc, argv))                // If an error occurs, the help message
   *                                                // will already be printed
   *    exit(-1);
   * \endcode
   *
   * The values of the options are generally stored in the member
   * variables of the appropriate object.  Simply use these, or copy
   * them out into the tool.  For example, to find the verbosity level
   * in the above example, check \c bopts->verbosity, and to access
   * the model read in, use \c mopts->model.
   *
   * One feature of boost::program_options is that the long-name
   * options may have a single letter short-cut.  Since
   * OptionsFramework integrates options from many different sources,
   * there is a possibility of short-cut collisions (i.e. two classes
   * using the same single letter code).  We therefore recommend the
   * following convention when providing these short-cuts:
   *   - LOOS wide common options should use lower-case single letters
   *   - Package-wide common options should use upper-case single letters
   *   - Tool specific options should always be long-names, unless the
   *     options are frequently used, in which case an upper-case
   *     letter should be used.
   *
   * The full set of command-line options is created by the
   * AggregateOptions class.  Using AggregateOptions::add(), different
   * OptionsPackage instances can be combined to build up the full set
   * of command-line options.  The order that the OptionsPackage
   * objects are added is important, as it will determine the order of
   * "positional" options as well as how the options are listed in the
   * help.  We recommend the following convention in options order:
   *    - BasicOptions
   *    - OutputPrefix
   *    - BasicSelection
   *    - Model/Trajectory (i.e. ModelWithCoords, BasicTrajectory,
   *      TrajectoryWithFrameIndices)
   *    - Tool specific options
   *    - RequiredArguments (see below)
   *
   * Frequently, a tool requires a number of command-line arguments
   * that are not optional.  The RequiredArguments class can be used
   * to incorporate these, rather than handling it within a
   * tool-specific subclass of OptionsPackage.  A RequiredArguments
   * object can build up a command line by passing a string tag and
   * description to RequiredArguments::addArgument() in the same order
   * that the arguments would appear on the command line.
   *
   * Another common command line format is to have an argument that
   * can appear one or more times.  This can also be handled by
   * RequiredArguments by using the
   * RequiredArguments::addVariableArguments() method.  If this
   * feature is used, the specified argument will appear
   * <em>after</em> all other required arguments.  Since the final
   * argument can appear one or more times, it will consume all
   * remaining command-line arguments.  This means that if you use
   * RequiredArguments, it should be the <em>last</em> OptionsPackage
   * included in the AggregateOptions object.
   *
   * Which OptionsPackage subclasses to use in a tool depends on what
   * the tool needs.  All tools, however, should include BasicOptions
   * to provide a help message and optional full-help message.  While
   * verbosity is defined in BasicOptions, tools do not need to
   * support it.  If a tool needs an output prefix, then include
   * OutputPrefix.  Similarly, if a tool needs a single selection to
   * determine which atoms to operate on, include a BasicSelection.
   * If the tool needs multiple selections, then you will need to
   * handle this explicitly with a tool-specific subclass of
   * OptionsPackage or via RequiredArguments.
   *
   * For tools that work with a model only, we recommend that you use
   * ModelWithCoords.  This will guarantee that the generatic
   * AtomicGroup has coordinates (e.g. the user specified a PSF file
   * rather than a PDF).  Tools that iterate over all frames in a
   * trajectory should use BasicTrajectory.  This class does provide a
   * --skip option that lets the user skip the first n-frames of a
   * trajectory (presumably, the equilibration stage).  The trajectory
   * it creates will be "primed" so that the first read will be the
   * n+1'th frame.  Finally, tools that can operate on a range of a
   * trajectory, or an arbitrary set of frames, should use
   * TrajectoryWithFrameIndices.  The
   * TrajectoryWithFrameIndices::frameList() method returns a vector
   * of unsigned integers that specify which frames of the trajectory
   * a tool should use.  For example:
   * \code
   * pTraj traj = trajectory_options->trajectory;
   * vector<uint> indices = trajectory_options->frameList();
   * for (vector<uint>::iterator i = indices.begin(); i != indices.end(); ++i)
   *    processTrajectoryFrame(traj->readFrame(*i));
   * \endcode
   *
   * <b>Writing a ToolOption class</b>
   *
   * The common idiom to add tool-specific options \e not covered by
   * the existing objects is to subclass \c OptionsPackage.  The
   * virtual memberfunctions (e.g. \c addGeneric(), \c addHidden(),
   * etc) mimic the steps you would normally take to use \c
   * boost::program_options.  It is only necessary to override the
   * functions required specifically for the tool.
   *
   * As a trivial example, here is a class that will add a bool option
   * to the command line,
   *
   * \code
   * namespace opts = loos::OptionsFramework;
   * namespace po = loos::OptionsFramework::po;
   *
   * class ToolOptions : public opts::OptionsPackage {
   * public:
   *
   *   void addGeneric(po::options_description& o) {
   *      o.add_options()
   *         ("flag", po::value<bool>(&flag)->default_value(false), "Set a flag");
   *   }
   *
   *   // print() included for completeness.  It's not necessary, but good form to
   *   // include...
   *   string print() const {
   *      std::ostringstream oss;
   *
   *      oss << boost::format("flag=%d") % flag;
   *      return(oss.str());
   *   }
   *
   *   bool flag;
   * };
   * \endcode
   *
   * This class can then be used as follows,
   * \code
   *
   * BasicOptions* bopts = new BasicOptions;
   * ModelWithCoords* mopts = new ModelWithCoords;
   * ToolOptions* topts = new ToolOptions;
   *
   * AggregateOptions options;
   * options.add(bopts).add(mopts).add(topts);
   *
   * if (!options.parse(argc, argv))
   *    exit(-1);
   * \endcode
   *
   * In more complex cases, \c check() can be used to validate the
   * command-line and \c postConditions() can be used to implement any
   * necessary post-processing of command-line options but still
   * within the context of the parser.  This means that an error will
   * result in a help message being printed.  For example, if a file
   * is required for input, the \c ifstream can be opened here.  If it
   * fails, then \c postConditions() will report the error and a help
   * message will be printed out.
   *
   * <hr>
   *
   * Notes:
   *   - Model and Trajectory options classes will create the
   *     appropriate model and trajectory objects which can be copied out
   *     of the respective options object for use in a tool.
   *
   *   - Pointers to OptionsPackage subclasses are used within the
   *     OptionsFramework.  Unlike most of LOOS, these are
   *     <em>not</em> shared pointers so the tool is responsible for
   *     managing them.  However, we anticipate that most will only
   *     use a small amount of memory and will need to exist for the
   *     life of a tool anyway, so deletion of the objects should not
   *     be a problem.
   *
   *   - OptionsPackages cannot be combined multiple times, e.g. using
   *     two BasicTrajectory objects because the tool needs to read
   *     from two different trajectories.  We anticipate these cases
   *     as being infrequent, and as such, there is no direct support
   *     for it in OptionsPackages.  Use a tool-specific
   *     OptionsPackage that handles the required args (matching as
   *     closely as possible the existing OptionsPackage options).
   *
   */
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

      //! Returns a string listing the encapsulated options, suitable for logging
      virtual std::string print() const { return(""); }

      //! Validates passed options, returning true if there is a problem or false if not
      /**
       * check() is typically used to validate positional options.
       * For example, if you tool needs a command line that looks
       * like:
       * \verbatim
       * tool [options] min-value max-value number-of-bins
       * \endverbatim
       * then check() would verify that min-value, max-value, and
       * number-of-bins were passed on the command-line.  Another
       * example is where there are mutually-exclusive options.  This
       * would be checked by check().
       */
      virtual bool check(po::variables_map& map) { return(false); }

      //! Post-processing of options returning true if there were no problems, otherwise false.
      /**
       * postConditions() is called after options parsing and
       * validation is complete.  This is a mechanism for a subclass
       * to do additional processing with the options it has been
       * provided.  For example, a model option subclass might read in
       * the specified model and copy coordinates from an optionally
       * specified file.
       *
       * Note that the return value from postConditions() is the
       * opposite of check().  Here, a true is returned if there are
       * no problems.
       */
      virtual bool postConditions(po::variables_map& map) {
        return(true);
      }

      //! Returns a slice of the example command-line in the help output
      /**
       * This is used specifically for positional options.  If your
       * tool has a command line that looks like:
       * \verbatim
       * tool [options] min-value max-value selection
       * \endverbatim
       * The required arguments min-value, max-value, and selection
       * are actually "hidden" options and will not be printed out
       * when boost::program_options generates its help message.  The
       * AggregateOptions class will build the full command line for
       * the help message by using the help() methods from each
       * OptionsPackage, printing out something like the above on the
       * command line as part of the help output.
       */

      virtual std::string help() const { return(""); }
    };


    // -------------------------------------------------

    //! Options common to all tools (including --fullhelp)
    class BasicOptions : public OptionsPackage {
    public:
      BasicOptions() : verbosity(0) { }
      BasicOptions(const int i) : verbosity(i) { }
      BasicOptions(const std::string& s) : verbosity(0), full_help(s) { }
      BasicOptions(const int i, const std::string& s) : verbosity(i), full_help(s) { }

      void setFullHelp(const std::string& s);

      int verbosity;
      std::string full_help;

    private:
      void addGeneric(po::options_description& opts);
      bool check(po::variables_map& map);

      std::string print() const;
    };

    // -------------------------------------------------

    //! Gets a string as prefix for output files (--prefix)
    class OutputPrefix : public OptionsPackage {
    public:
      OutputPrefix() : prefix("output"), label("Output prefix") { }
      OutputPrefix(const std::string& s) : prefix(s), label("Output prefix") { }
      OutputPrefix(const std::string& s, const std::string& t) : prefix(s), label(t) { }


      std::string prefix;
      std::string label;

    private:
      void addGeneric(po::options_description& opts);
      std::string print() const;
    };

    // -------------------------------------------------

    //! Provides a single LOOS selection (--selection)
    class BasicSelection : public OptionsPackage {
    public:
      BasicSelection() : selection("all"), label("Which atoms to use") { }
      BasicSelection(const std::string& sel) : selection(sel), label("Which atoms to use") { }
      BasicSelection(const std::string& sel, const std::string& lbl) :
        selection(sel),
        label(lbl)
      { }

      std::string selection;
      std::string label;

    private:
      void addGeneric(po::options_description& opts);
      std::string print() const;
    };


    // -------------------------------------------------

    //! Provides a mechanism for controlling how to split an AtomicGroup
    class BasicSplitBy : public OptionsPackage {
    public:
      enum SplitType { NONE, MOLECULE, SEGID, RESIDUE };


      BasicSplitBy() : split_method("mol"), label("Split selection by (none, mol, segid, res)"), split_type(MOLECULE) { }
      BasicSplitBy(const std::string& method) : split_method(method), label("Split selection by (none, mol, segid, res)"), split_type(methodToType(method)) { }
      BasicSplitBy(const std::string& method, const std::string& lbl) :
        split_method(method),
        label(lbl),
	split_type(methodToType(method))
      { }

      std::string split_method;
      std::string label;
      SplitType split_type;

      std::vector<AtomicGroup> split(const AtomicGroup& grp);

    private:
      SplitType methodToType(const std::string& name);
      void addGeneric(po::options_description& opts);
      bool postConditions(po::variables_map& map);
      std::string print() const;
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

      std::string model_name, coords_name, model_type;

      AtomicGroup model;

    private:
      void addGeneric(po::options_description& opts);
      void addHidden(po::options_description& opts);
      void addPositional(po::positional_options_description& pos);

      bool check(po::variables_map& map);

      bool postConditions(po::variables_map& map);

      std::string help() const;
      std::string print() const;
    };

    // -------------------------------------------------

    //! Request Two models with coordinates
    /**
     * Since not all formats have coordinates (i.e. PSF),
     * the coordinates can be taken from an alternate file using the
     * -c or --coordinates option.  Also adds a positional argument
     * for the model description.
     **/
    class TwoModelsWithCoords : public OptionsPackage {
    public:
      TwoModelsWithCoords() : coords1_name(""), desc1("model1"),
                              coords2_name(""), desc2("model2") { }

      TwoModelsWithCoords(const std::string& d1,
                          const std::string& d2) :
        coords1_name(""), desc1(d1),
        coords2_name(""), desc2(d2)
      { }

      std::string model1_name, model1_type, coords1_name, desc1;
      std::string model2_name, model2_type, coords2_name, desc2;

      AtomicGroup model1;
      AtomicGroup model2;

    private:
      void addGeneric(po::options_description& opts);
      void addHidden(po::options_description& opts);
      void addPositional(po::positional_options_description& pos);

      bool check(po::variables_map& map);

      bool postConditions(po::variables_map& map);

      std::string help() const;
      std::string print() const;
    };




    // -------------------------------------------------

    //! Basic trajectory with a --skip option
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


      unsigned int skip;
      std::string model_name, model_type, traj_name, traj_type;

      //! Model that describes the trajectory
      AtomicGroup model;

      //! The trajectory, primed by the --skip value (if specified)
      pTraj trajectory;


    private:
      void addGeneric(po::options_description& opts);
      void addHidden(po::options_description& opts);

      void addPositional(po::positional_options_description& pos);

      bool check(po::variables_map& map);

      bool postConditions(po::variables_map& map);

      std::string help() const;
      std::string print() const;
    };



    // -------------------------------------------------

    //! Trajectory with either a --range or --skip
    /**
     * Adds a model and trajectory argument to the command line, and
     * provides --skip (-k) and --range (-r) options for specifying
     * which frames of the trajectory to operate over.
     *
     * Use TrajectoryWithFrameIndices::frameList() to get a vector of
     * unsigned ints representing which frames the user requested.
     **/
    class TrajectoryWithFrameIndices : public OptionsPackage {
    public:
      TrajectoryWithFrameIndices() : skip(0), stride(1), frame_index_spec("") { }

      //! Returns the list of frames the user requested
      std::vector<uint> frameList() const;

      unsigned int skip, stride;
      std::string frame_index_spec;
      std::string model_name, model_type, traj_name, traj_type;


      //! Model that describes the trajectory
      AtomicGroup model;

      //! The trajectory
      pTraj trajectory;


    private:
      void addGeneric(po::options_description& opts);
      void addHidden(po::options_description& opts);

      void addPositional(po::positional_options_description& pos);

      bool check(po::variables_map& map);

      bool postConditions(po::variables_map& map);

      std::string help() const;
      std::string print() const;
    };



    // -------------------------------------------------




    //! Multiple trajectories as one trajectory via MultiTrajectory
    /**
     *
     * Handles multiple trajectories as one large virtual trajectory via
     * a MultiTrajectory object.  Skip and stride are applied to each
     * contained sub-trajectory, not the overall trajectory.
     *
     * Internally, it instantiates a MultiTrajectory object in mtraj, and
     * then wraps it in a pTraj that is set to not deallocate the source
     * object.
     *
     **/
    class MultiTrajOptions : public OptionsPackage {
    public:
      MultiTrajOptions() : skip(0), stride(1) { }


      uint skip;
      uint stride;
      std::vector< std::string > traj_names;
      std::string model_name, model_type, frame_index_spec;

      std::vector<uint> frameList() const;

      //! Model that describes the trajectory
      AtomicGroup model;

      //! The original collection of trajectories
      MultiTrajectory mtraj;

      //! The trajectory
      pTraj trajectory;

      //! A table of information about the contained trajectories
      std::string trajectoryTable() const;


    private:
      void addGeneric(po::options_description& opts);
      void addHidden(po::options_description& opts);

      void addPositional(po::positional_options_description& pos);

      bool check(po::variables_map& map);

      bool postConditions(po::variables_map& map);

      std::string help() const;
      std::string print() const;
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
     *
     * The value returned for options are strings and must be parsed
     * into the appropriate type.  A simple way to do this is to use
     * loos::parseStringAs<>.
     *
     * Example:
     * \code
     * RequiredArguments* ropts = new RequiredArguments;
     * ropts->addArgument("name", "Name of data file");
     * ropts->addArgument("scale", "Scaling to apply to data");
     * ...
     * string name = ropts->value("name");
     * double scale = parseStringAs<double>(ropts->value("scale"));
     * \endcode
     *
     * The case where one or more arguments are required is supported
     * via the RequiredArguments::addVariableArguments() and
     * RequiredArguments::variableValues().  If this feature is used,
     * then the RequiredArguments object <em>must</em> be the last
     * OptionsPackage chained together in the AggregateOptions object.
    **/
    class RequiredArguments : public OptionsPackage {
      typedef std::pair<std::string, std::string>    StringPair;
    public:
      RequiredArguments() : vargs_set(false) { }

      RequiredArguments(const std::string& argname, const std::string& argdesc) : vargs_set(false)
      {
        addArgument(argname, argdesc);
      }


      //! Add a required argument given a name (tag) and a description (currently unused)
      void addArgument(const std::string& name, const std::string& description);

      //! Add a required argument that can be an arbitrary number of items
      /**
       * This argument will always appear at the end of the command
       * line, after all other required arguments.  It also means that
       * the RequiredOptions object should be the last one added to
       * the AggregateOptions object, otherwise any subsequence
       * positional arguments will be missed.
       *
       * Example:
       * \code
       * RequiredArguments* ropts = new RequiredArguments;
       * ropts->addVariableArguments("selection", "selection");
       * ...
       * vector<string> selections = ropts->variableValues("selection");
       * \endcode
       */
      void addVariableArguments(const std::string& name, const std::string& description);

      //! Retrieve the value for an argument
      std::string value(const std::string& s) const;

      //! Retrieve the variable-number argument
      std::vector<std::string> variableValues(const std::string& s) const;

    private:
      void addHidden(po::options_description& o);
      void addPositional(po::positional_options_description& pos);
      bool check(po::variables_map& map);
      bool postConditions(po::variables_map& map);
      std::string help() const;
      std::string print() const;

    private:
      bool vargs_set;
      std::vector<StringPair> arguments;
      StringPair variable_arguments;
      po::variables_map held_map;
    };


    // ----------------------------------------------------------------------
    class OutputTrajectoryOptions : public OptionsPackage {
    public:
      OutputTrajectoryOptions() : name("output.dcd"), label("Output Trajectory"), append(false) {}
      OutputTrajectoryOptions(const std::string& s) : name(s), label("Output Trajectory"), append(false) {}
      OutputTrajectoryOptions(const std::string& s, const bool appending) : name(s), label("Output Trajectory"), append(appending) {}


      std::string name;
      std::string label;
      bool append;
      std::string type;
      std::string basename;
      pTrajectoryWriter outraj;


    private:
      void addGeneric(po::options_description& opts);
      void addHidden(po::options_description& opts);

      void addPositional(po::positional_options_description& map);

      bool check(po::variables_map& map);

      bool postConditions(po::variables_map& map);

      std::string help() const;
      std::string print() const;

    };

    // ----------------------------------------------------------------------

    class OutputTrajectoryTypeOptions : public OptionsPackage {
    public:
      OutputTrajectoryTypeOptions() :
	label("Output Trajectory Type"),
	append(false),
	type("dcd") {}

      OutputTrajectoryTypeOptions(const std::string& s) :
	label("Output Trajectory Type"),
	append(false),
	type(s) {}

      OutputTrajectoryTypeOptions(const std::string& s, const bool appending) :
	label("Output Trajectory Type"),
	append(appending), type(s) {}

      pTrajectoryWriter createTrajectory(const std::string& prefix);


      std::string label;
      bool append;
      std::string type;

    private:
      void addGeneric(po::options_description& opts);
      std::string print() const;
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

      //! Explicitly set the program name (for help message and printing)
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

      //! Returns a string representing the option values in all
      //! contained OptionPackage objects
      std::vector<std::string> print() const;

      //! Displays the help for this tool
      void showHelp();

    private:
      std::string program_name;
      std::string config_name;

      po::options_description generic;
      po::options_description hidden;
      po::options_description command_line;
      po::positional_options_description pos;
      po::variables_map vm;

      std::vector<OptionsPackage *> options;


      void setupOptions();
    };

    class WeightsOptions : public OptionsPackage {
    public:
      WeightsOptions(): has_weights(false) { }

      std::string weights_name;
      std::string list_name;
      Weights* weights;
      bool has_weights;

    private:
      void addGeneric(po::options_description& opts);
      bool postConditions(po::variables_map& map);
    };


  };
};


#endif
