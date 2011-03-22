#include <OptionsFramework.hpp>


namespace loos {
  namespace DensityTools {
    namespace OptionsFramework {

      void BasicOptions::addGeneric(po::options_description& opts) {
        opts.add_options()
          ("help", "Produce this message")
          ("verbosity,v", po::value<int>(&verbosity)->default_value(verbosity), "Verbosity");
      }

      string BasicOptions::print() const {
        ostringstream oss;
        oss << "# verbosity=" << verbosity ;
        return(oss.str());
      }

      // -------------------------------------------------------

      void OutputPrefixOptions::addGeneric(po::options_description& opts) {
        opts.add_options()
          ("prefix,p", po::value<string>(&prefix)->default_value(prefix), "Output prefix");
      }

      string OutputPrefixOptions::print() const {
        ostringstream oss;
        oss << "# prefix='" << prefix << "'\n";
        return(oss.str());
      }
      
      // -------------------------------------------------------

      void BasicSelectionOptions::addGeneric(po::options_description& opts) {
        opts.add_options()
          ("selection,s", po::value<string>(&selection)->default_value(selection), "Which atoms to use");
      }

      string BasicSelectionOptions::print() const {
        ostringstream oss;
        oss << "# selection='" << selection << "'\n";
        return(oss.str());
      }


      // -------------------------------------------------------


      void BasicTrajectoryOptions::addGeneric(po::options_description& opts) {
        opts.add_options()
          ("skip,S", po::value<unsigned int>(&skip)->default_value(skip), "Number of frames to skip")
          ("range,r", po::value<string>(&frame_index_spec), "Which frames to use (matlab style range)");
      };

      void BasicTrajectoryOptions::addHidden(po::options_description& opts) {
        opts.add_options()
          ("model", po::value<string>(&model_name), "Model filename")
          ("traj", po::value<string>(&traj_name), "Trajectory filename");
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
          cerr << "Error- you cannot specify both a skip and a frame range...I might get confused!\n";
          return(false);
        }

        return(true);
      }

      string BasicTrajectoryOptions::help() const { return("model trajectory"); }
      string BasicTrajectoryOptions::print() const {
        ostringstream oss;
        oss << "# model='" << model_name << "', traj='" << traj_name << "', ";
        if (skip > 0)
          oss << "skip=" << skip;
        else
          oss << "range=" << frame_index_spec;
        oss << endl;

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

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(command_line).positional(pos).run(), vm);
        po::notify(vm);

        bool show_help = vm.count("help");
        if (!show_help)
          for (vOpts::iterator i = options.begin(); i != options.end() && !show_help; ++i)
            show_help = (*i)->check(vm);

        if (show_help) {
          cout << "Usage- " << argv[0] << " [options] ";
          for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
            cout << (*i)->help() << " ";
          cout << endl;
          cout << generic;
          return(false);
        }

        for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
          if (!(*i)->postConditions(vm))
            return(false);

        return(true);
      
      }

      string AggregateOptions::print() const {
        string result;
    
        for (vOpts::const_iterator i = options.begin(); i != options.end(); ++i)
          result += (*i)->print();

        return(result);
      }





    };
  };
};
