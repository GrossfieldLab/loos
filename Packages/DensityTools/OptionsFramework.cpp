#include <OptionsFramework.hpp>

namespace loos {
  namespace DensityTools {
    namespace OptionsFramework {

      void BasicOptions::addGeneric(po::options_description& opts) {
        opts.add_options()
          ("help", "Produce this message")
          ("verbosity,v", po::value<int>(&verbosity)->default_value(verbosity), "Verbosity");
      }

      std::string BasicOptions::print() const {
        std::ostringstream oss;
        oss << "# verbosity=" << verbosity ;
        return(oss.str());
      }

      // -------------------------------------------------------

      void OutputPrefixOptions::addGeneric(po::options_description& opts) {
        opts.add_options()
          ("prefix,p", po::value<std::string>(&prefix)->default_value(prefix), "Output prefix");
      }

      std::string OutputPrefixOptions::print() const {
        std::ostringstream oss;
        oss << "# prefix='" << prefix << "'\n";
        return(oss.str());
      }
      
      // -------------------------------------------------------

      void BasicSelectionOptions::addGeneric(po::options_description& opts) {
        opts.add_options()
          ("selection,s", po::value<std::string>(&selection)->default_value(selection), "Which atoms to use");
      }

      std::string BasicSelectionOptions::print() const {
        std::ostringstream oss;
        oss << "# selection='" << selection << "'\n";
        return(oss.str());
      }


      // -------------------------------------------------------


      void BasicTrajectoryOptions::addGeneric(po::options_description& opts) {
        opts.add_options()
          ("skip,S", po::value<unsigned int>(&skip)->default_value(skip), "Number of frames to skip")
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
        oss << "# model='" << model_name << "', traj='" << traj_name << "', ";
        if (skip > 0)
          oss << "skip=" << skip;
        else
          oss << "range=" << frame_index_spec;
        oss << std::endl;

        return(oss.str());
      }


      // -------------------------------------------------------


      AggregateOptions& AggregateOptions::addOptions(OptionsPackage* pack) {
        options.push_back(pack); return(*this);
      }


      bool AggregateOptions::parseOptions(int argc, char *argv[]) {

        po::options_description generic("Allowed options");
        for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
          (*i)->addGeneric(generic);

        po::options_description hidden("Hidden options");
        for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
          (*i)->addHidden(hidden);

        po::options_description command_line;
        command_line.add(generic).add(hidden);

        po::positional_options_description pos;
        for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
          (*i)->addPositional(pos);

        bool show_help = false;

        po::variables_map vm;
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
          std::cout << "Usage- " << argv[0] << " [options] ";
          for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
            std::cout << (*i)->help() << " ";
          std::cout << std::endl;
          std::cout << generic;
          return(false);
        }

        for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
          if (!(*i)->postConditions(vm))
            return(false);

        return(true);
      
      }

      std::string AggregateOptions::print() const {
        std::string result;
    
        for (vOpts::const_iterator i = options.begin(); i != options.end(); ++i)
          result += (*i)->print();

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

    };
  };
};
