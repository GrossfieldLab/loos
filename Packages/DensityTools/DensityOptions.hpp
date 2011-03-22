#if !defined(DENSITY_OPTIONS_HPP)
#define DENSITY_OPTIONS_HPP


#include <loos.hpp>
#include <DensityGrid.hpp>
#include <internal-water-filter.hpp>
#include <OptionsFramework.hpp>

namespace loos {
  namespace DensityTools {

    namespace OptionsFramework {

      class BasicWaterOptions : public OptionsPackage {
      public:
        BasicWaterOptions() :
          pad(1.0),
          radius(10.0),
          water_string("name == 'OH2'"),
          prot_string("name == 'CA'"),
          filter_mode("axis") { }
        
        void addGeneric(po::options_description& opts) {
          opts.add_options()
            ("pad,P", po::value<double>(&pad)->default_value(pad), "Pad (for bounding box)")
            ("radius,r", po::value<double>(&radius)->default_value(radius), "Radius (for principal axis filter)")
            ("zrange", po::value<string>(), "Clamp the volume to integrate over in Z (min:max)")
            ("water,w", po::value<string>(&water_string)->default_value(water_string), "Water selection")
            ("prot,p", po::value<string>(&prot_string)->default_value(prot_string), "Protein selection")
            ("grid,g", po::value<string>(), "Name of grid to use in grid-mode")
            ("mode,m", po::value<string>(&mode)->default_value(filter_mode), "Mode (axis|box|grid)");
        }

        bool postConditions(po::variables_map& map) {
          if (mode == "axis") {
            filter_func = new WaterFilterAxis(radius);
          } else if (mode == "box") {
            filter_func = new WaterFilterBox(pad);
          } else if (mode == "grid") {
            if (! vm.count("grid")) {
              cerr << "ERROR - you must specify a grid to use when using grid-mode\n";
              return(false);
            }

            string grid_name = vm["grid"].as<string>();
            ifstream ifs(grid_name.c_str());
            ifs >> the_grid;
            cerr << "Read in grid with size " << the_grid.gridDims() << endl;
      
            filter_func = new WaterFilterBlob(the_grid);

          } else {
            cerr << "ERROR - unknown mode " << mode << endl;
            return(false);
          }
          
          // Handle "decoration"
          if (vm.count("zrange")) {
            double zmin, zmax;
            string s = vm["zrange"].as<string>();
            int i = sscanf(s.c_str(), "%lf:%lf", &zmin, &zmax);
            if (i != 2) {
              cerr << boost::format("ERROR - unable to parse range '%s'\n") % s;
              return(false);
            }

            filter_func = new ZClippedWaterFilter(filter_func, zmin, zmax);
          }
          
          return(true);
        }

        double zmin, zmax;
        double pad;
        double radius;
        string water_string, prot_string, grid_name, filter_mode;
        DensityGrid<int> the_grid;
        WaterFilterBase* filter_func;
      };
      

    };

  
  };
};


#endif
