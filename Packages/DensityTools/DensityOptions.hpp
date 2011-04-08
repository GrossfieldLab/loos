#if !defined(DENSITY_OPTIONS_HPP)
#define DENSITY_OPTIONS_HPP

#include <boost/format.hpp>


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
          grid_name(""),
          filter_mode("axis"),
          bulked_spec(""),
          zrange_spec("")
        { }
        
        void addGeneric(po::options_description& opts) {
          opts.add_options()
            ("water,W", po::value<std::string>(&water_string)->default_value(water_string), "Water selection")
            ("prot,P", po::value<std::string>(&prot_string)->default_value(prot_string), "Protein selection")
            ("pad", po::value<double>(&pad)->default_value(pad), "Pad (for bounding box)")
            ("bulked", po::value<std::string>(&bulked_spec)->default_value(bulked_spec), "Add bulk water (z-slices between cutoff and bounding box) [pad,zmin:zmax]")
            ("radius,R", po::value<double>(&radius)->default_value(radius), "Radius (for principal axis filter)")
            ("zrange", po::value<std::string>(&zrange_spec)->default_value(zrange_spec), "Clamp the volume to integrate over in Z (min:max)")
            ("grid,G", po::value<std::string>(&grid_name)->default_value(grid_name), "Name of grid to use in grid-mode (for internal waters)")
            ("mode,M", po::value<std::string>(&filter_mode)->default_value(filter_mode), "Mode (axis|box|grid)");
        }

        bool postConditions(po::variables_map& map) {
          if (filter_mode == "axis") {
            filter_func = new WaterFilterAxis(radius);
          } else if (filter_mode == "box") {
            filter_func = new WaterFilterBox(pad);
          } else if (filter_mode == "grid") {
            if (grid_name.empty()) {
              std::cerr << "ERROR - you must specify a grid to use when using grid-mode\n";
              return(false);
            }

            std::ifstream ifs(grid_name.c_str());
            ifs >> the_grid;
            std::cerr << "Read in grid with size " << the_grid.gridDims() << std::endl;
      
            filter_func = new WaterFilterBlob(the_grid);

          } else {
            std::cerr << "ERROR - unknown mode " << filter_mode << std::endl;
            return(false);
          }
          
          // Handle "decoration"
          if (!zrange_spec.empty()) {
            double zmin, zmax;
            int i = sscanf(zrange_spec.c_str(), "%lf:%lf", &zmin, &zmax);
            if (i != 2) {
              std::cerr << boost::format("ERROR - unable to parse range '%s'\n") % zrange_spec;
              return(false);
            }

            filter_func = new ZClippedWaterFilter(filter_func, zmin, zmax);
          }

          if (!bulked_spec.empty()) {
            double zmin, zmax, pad;
            int i = sscanf(bulked_spec.c_str(), "%lf,%lf:%lf", &pad, &zmin, &zmax);
            if (i != 3) {
              std::cerr << boost::format("ERROR - unable to parse bulk range '%s'\n") % bulked_spec;
              return(false);
            }
            
            filter_func = new BulkedWaterFilter(filter_func, pad, zmin, zmax);
          }
          
          
          return(true);
        }

        std::string print() const {
          std::ostringstream oss;

          oss << boost::format("water='%s', prot='%s', pad=%f, bulked='%s', radius=%f, zrange='%s', grid='%s', mode='%s'")
            % water_string
            % prot_string
            % pad
            % (bulked_spec.empty() ? "undef" : bulked_spec)
            % radius
            % (zrange_spec.empty() ? "undef" : zrange_spec)
            % (grid_name.empty() ? "undef" : grid_name)
            % filter_mode;

          return(oss.str());
        }
        

        double zmin, zmax;
        double pad;
        double radius;
        std::string water_string, prot_string, grid_name, filter_mode;
        std::string bulked_spec, zrange_spec;
        DensityGrid<int> the_grid;
        WaterFilterBase* filter_func;
      };
      

    };

  
  };
};


#endif
