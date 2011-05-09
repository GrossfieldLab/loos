#include <OptionsFramework.hpp>

namespace loos {
  namespace OptionsFramework {

    void BasicOptions::addGeneric(po::options_description& opts) {
      opts.add_options()
        ("help", "Produce this message")
        ("verbosity,v", po::value<int>(&verbosity)->default_value(verbosity), "Verbosity");
    }

    std::string BasicOptions::print() const {
      std::ostringstream oss;
      oss << "verbosity=" << verbosity ;
      return(oss.str());
    }

    // -------------------------------------------------------

    void OutputPrefixOptions::addGeneric(po::options_description& opts) {
      opts.add_options()
        ("prefix,p", po::value<std::string>(&prefix)->default_value(prefix), "Output prefix");
    }

    std::string OutputPrefixOptions::print() const {
      std::ostringstream oss;
      oss << "prefix='" << prefix << "'";
      return(oss.str());
    }
      
    // -------------------------------------------------------

    void BasicSelectionOptions::addGeneric(po::options_description& opts) {
      opts.add_options()
        ("selection,s", po::value<std::string>(&selection)->default_value(selection), "Which atoms to use");
    }

    std::string BasicSelectionOptions::print() const {
      std::ostringstream oss;
      oss << "selection='" << selection << "'";
      return(oss.str());
    }

    // -------------------------------------------------------

    void ModelWithCoordsOptions::addGeneric(po::options_description& opts) {
      opts.add_options()
        ("coordinates,c", po::value<std::string>(&coords_name)->default_value(coords_name), "File to use for coordinates");
    }

    void ModelWithCoordsOptions::addHidden(po::options_description& opts) {
      opts.add_options()
        ("model", po::value<std::string>(&model_name), "Model Filename");
    }


    void ModelWithCoordsOptions::addPositional(po::positional_options_description& pos) {
      pos.add("model", 1);
    }


    bool ModelWithCoordsOptions::check(po::variables_map& map) {
      return(!map.count("model"));
    }

    std::string ModelWithCoordsOptions::help() const { return("model"); }


    std::string ModelWithCoordsOptions::print() const {
      std::ostringstream oss;

      oss << boost::format("model='%s'") % model_name;
      if (!coords_name.empty())
        oss << boost::format(", coords='%s'") % coords_name;

      return(oss.str());
    }



    // -------------------------------------------------------


    void BasicTrajectoryOptions::addGeneric(po::options_description& opts) {
      opts.add_options()
        ("skip,k", po::value<unsigned int>(&skip)->default_value(skip), "Number of frames to skip")
        ("range,r", po::value<std::string>(&frame_index_spec), "Which frames to use (matlab style range)");
    };

    void BasicTrajectoryOptions::addHidden(po::options_description& opts) {
      opts.add_options()
        ("model", po::value<std::string>(&model_name), "Model filename")
        ("traj", po::value<std::string>(&traj_name), "Trajectory filename");
    }

    void BasicTrajectoryOptions::addPositional(po::positional_options_description& pos) {
      pos.add("model", 1);
      pos.add("traj", 1);
    }

    bool BasicTrajectoryOptions::check(po::variables_map& map) {
      return(! (map.count("model") && map.count("traj")));
    }

    bool BasicTrajectoryOptions::postConditions(po::variables_map& map) {
      if (skip > 0 && !frame_index_spec.empty()) {
        std::cerr << "Error- you cannot specify both a skip and a frame range...I might get confused!\n";
        return(false);
      }

      return(true);
    }

    std::string BasicTrajectoryOptions::help() const { return("model trajectory"); }
    std::string BasicTrajectoryOptions::print() const {
      std::ostringstream oss;
      oss << boost::format("model='%s', traj='%s', ") % model_name % traj_name;
      if (skip > 0)
        oss << "skip=" << skip;
      else
        oss << "range=" << frame_index_spec;

      return(oss.str());
    }


    // -------------------------------------------------------


    AggregateOptions& AggregateOptions::add(OptionsPackage* pack) {
      options.push_back(pack); return(*this);
    }


    void AggregateOptions::setupOptions() {
      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        (*i)->addGeneric(generic);

      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        (*i)->addHidden(hidden);

      command_line.add(generic).add(hidden);

      pos = po::positional_options_description();
      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        (*i)->addPositional(pos);
    }

    void AggregateOptions::showHelp() {
      std::cout << "Usage- " << program_name << " [options] ";
      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        std::cout << (*i)->help() << " ";
      std::cout << std::endl;
      std::cout << generic;
    }


    bool AggregateOptions::parse(int argc, char *argv[]) {
      if (program_name.empty())
        program_name = std::string(argv[0]);

      setupOptions();
      bool show_help = false;

      try {
        po::store(po::command_line_parser(argc, argv).
                  options(command_line).positional(pos).run(), vm);
        po::notify(vm);
      }
      catch (std::exception& e) {
        std::cerr << "Error- " << e.what() << std::endl;
        show_help = true;
      }

      if (!show_help)
        show_help = vm.count("help");

      if (!show_help)
        for (vOpts::iterator i = options.begin(); i != options.end() && !show_help; ++i)
          show_help = (*i)->check(vm);

      if (show_help) {
        showHelp();
        return(false);
      }

      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        if (!(*i)->postConditions(vm)) {
          showHelp();
          return(false);
        }

      return(true);
      
    }

    std::string AggregateOptions::print() const {
      std::string result(program_name);

      result += ": ";
    
      for (vOpts::const_iterator i = options.begin(); i != options.end(); ++i) {
        result += (*i)->print();
        if (i != options.end() - 1)
          result += ", ";
      }

      return(result);
    }



    std::vector<uint> assignFrameIndices(pTraj& traj, const std::string& desc, const uint skip = 0) {
      std::vector<uint> frames;

      if (desc.empty())
        for (uint i=skip; i<traj->nframes(); ++i)
          frames.push_back(i);
      else
        frames = parseRangeList<uint>(desc);

      return(frames);
          
    }

    AtomicGroup loadStructureWithCoords(const std::string model_name, const std::string coord_name = std::string("")) {
      AtomicGroup model = createSystem(model_name);
      if (!coord_name.empty()) {
        AtomicGroup coords = createSystem(coord_name);
        model.copyCoordinates(coords);
      }

      if (! model.hasCoords())
        throw(LOOSError("Error- no coordinates found in specified model(s)"));

      return(model);
    }

  };
};

