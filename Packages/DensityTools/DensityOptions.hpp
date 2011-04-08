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
            ("water,W", po::value<std::string>(&water_string)->default_value(water_string), "Water selection")
            ("prot,P", po::value<std::string>(&prot_string)->default_value(prot_string), "Protein selection")
            ("pad", po::value<double>(&pad)->default_value(pad), "Pad (for bounding box)")
            ("bulked", po::value<std::string>(), "Add bulk water (z-slices between cutoff and bounding box) [pad,zmin:zmax]")
            ("radius,R", po::value<double>(&radius)->default_value(radius), "Radius (for principal axis filter)")
            ("zrange", po::value<std::string>(), "Clamp the volume to integrate over in Z (min:max)")
            ("grid,G", po::value<std::string>(), "Name of grid to use in grid-mode (for internal waters)")
            ("mode,M", po::value<std::string>(&filter_mode)->default_value(filter_mode), "Mode (axis|box|grid)");
        }

        bool postConditions(po::variables_map& map) {
          if (filter_mode == "axis") {
            filter_func = new WaterFilterAxis(radius);
          } else if (filter_mode == "box") {
            filter_func = new WaterFilterBox(pad);
          } else if (filter_mode == "grid") {
            if (! map.count("grid")) {
              std::cerr << "ERROR - you must specify a grid to use when using grid-mode\n";
              return(false);
            }

            std::string grid_name = map["grid"].as<std::string>();
            std::ifstream ifs(grid_name.c_str());
            ifs >> the_grid;
            std::cerr << "Read in grid with size " << the_grid.gridDims() << std::endl;
      
            filter_func = new WaterFilterBlob(the_grid);

          } else {
            std::cerr << "ERROR - unknown mode " << filter_mode << std::endl;
            return(false);
          }
          
          // Handle "decoration"
          if (map.count("zrange")) {
            double zmin, zmax;
            std::string s = map["zrange"].as<std::string>();
            int i = sscanf(s.c_str(), "%lf:%lf", &zmin, &zmax);
            if (i != 2) {
              std::cerr << boost::format("ERROR - unable to parse range '%s'\n") % s;
              return(false);
            }

            filter_func = new ZClippedWaterFilter(filter_func, zmin, zmax);
          }

          if (map.count("bulked")) {
            double zmin, zmax, pad;
            std::string s = map["water"].as<std::string>();
            int i = sscanf(s.c_str(), "%lf,%lf:%lf", &pad, &zmin, &zmax);
            if (i != 3) {
              std::cerr << boost::format("ERROR - unable to parse bulk range '%s'\n") % s;
              return(false);
            }
            
            filter_func = new BulkedWaterFilter(filter_func, pad, zmin, zmax);
          }
          
          
          return(true);
        }

        double zmin, zmax;
        double pad;
        double radius;
        std::string water_string, prot_string, grid_name, filter_mode;
        DensityGrid<int> the_grid;
        WaterFilterBase* filter_func;
      };
      

    };

  
  };
};


#endif
