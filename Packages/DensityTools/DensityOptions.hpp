// -------------------------------------------------
// Density Package specific tool options
// -------------------------------------------------

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



#if !defined(LOOS_DENSITY_OPTIONS_HPP)
#define LOOS_DENSITY_OPTIONS_HPP

#include <boost/format.hpp>


#include <loos.hpp>
#include <DensityGrid.hpp>
#include <internal-water-filter.hpp>

namespace loos {

  namespace OptionsFramework {

      //! Options specific to tools that work with water/internal-water
      /**
       * Allows definition of water, protein, padding, among others.
       * Also includes a "factory" for setting how internal waters are defined
       * (selected).
       *
       * Mode can be: axis (distance from principal axis), box (a
       * bounding box), and grid (grid-mask).
       *
       * Decorations include Z-range (clamping to a set-zrange),
       * and "bulked", which adds back in the entire system above and
       * below a z-plane.
       **/
      class BasicWater : public OptionsPackage {
        /*
          Grids used as masks are integer grids.  They should only contain values from
          [0,max_grid_mask_value].
         */
        static const int max_grid_mask_value = 0xff;
      public:
        BasicWater() :
          pad(1.0),
          radius(10.0),
	  contacts(5),
          water_string("name == 'OH2'"),
          prot_string("name == 'CA'"),
          grid_name(""),
          filter_mode("axis"),
          bulked_spec(""),
          zrange_spec(""),
          filter_func(0)
        { }
        
        void addGeneric(po::options_description& opts) {
          opts.add_options()
            ("water,W", po::value<std::string>(&water_string)->default_value(water_string), "Water selection")
            ("prot,P", po::value<std::string>(&prot_string)->default_value(prot_string), "Protein selection")
            ("pad", po::value<double>(&pad)->default_value(pad), "Pad (for bounding box)")
            ("bulked", po::value<std::string>(&bulked_spec)->default_value(bulked_spec), "Add bulk water (z-slices between cutoff and bounding box) [pad,zmin:zmax]")
            ("radius,R", po::value<double>(&radius)->default_value(radius), "Radius (for principal axis filter and radius filter)")
            ("zrange", po::value<std::string>(&zrange_spec)->default_value(zrange_spec), "Clamp the volume to integrate over in Z (min:max)")
            ("grid,G", po::value<std::string>(&grid_name)->default_value(grid_name), "Name of grid to use in grid-mode (for internal waters)")
	    ("contacts,C", po::value<uint>(&contacts)->default_value(contacts), "Minimum number of contacts required for contacts filter")
            ("mode,M", po::value<std::string>(&filter_mode)->default_value(filter_mode), "Mode (axis|box|radius|contacts|core|grid)");
        }

        bool postConditions(po::variables_map& map) {
          if (filter_mode == "axis") {
            filter_func = new DensityTools::WaterFilterAxis(radius);
          } else if (filter_mode == "box") {
            filter_func = new DensityTools::WaterFilterBox(pad);
          } else if (filter_mode == "radius") {
            filter_func = new DensityTools::WaterFilterRadius(radius);
	  } else if (filter_mode == "contacts") {
	    filter_func = new DensityTools::WaterFilterContacts(radius, contacts);
	  } else if (filter_mode == "core") {
	    filter_func = new DensityTools::WaterFilterCore(radius);
          } else if (filter_mode == "grid") {
            if (grid_name.empty()) {
              std::cerr << "ERROR - you must specify a grid to use when using grid-mode\n";
              return(false);
            }

            std::ifstream ifs(grid_name.c_str());
	    if (!ifs) {
		std::cerr << "Error- cannot open grid file '" << grid_name << "' for reading." << std::endl;
		return(false);
	    }
	    
            ifs >> the_grid;
            std::cerr << "Read in grid with size " << the_grid.gridDims() << std::endl;
            if (!validateGridMask(the_grid)) {
              std::cerr << "\n***WARNING***WARNING***WARNING***\n\n"
                        << "The grid '" << grid_name
                        << "' does not appear to be a grid mask.  Your output will be suspect!\n"
                        << "Make sure the grid is an integer grid (e.g. from blobid or pick_block).\n"
                        << "If there are many unique blobs in the grid, then you may see this warning\n"
                        << "in which case you should reconsider how you ran blobid.\n"
                        << "Use --fullhelp for more information.\n"
                        << "\n***WARNING***WARNING***WARNING***\n\n";
            }

            
      
            filter_func = new DensityTools::WaterFilterBlob(the_grid);

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

            filter_func = new DensityTools::ZClippedWaterFilter(filter_func, zmin, zmax);
          }

          if (!bulked_spec.empty()) {
            double zmin, zmax, pad;
            int i = sscanf(bulked_spec.c_str(), "%lf,%lf:%lf", &pad, &zmin, &zmax);
            if (i != 3) {
              std::cerr << boost::format("ERROR - unable to parse bulk range '%s'\n") % bulked_spec;
              return(false);
            }
            
            filter_func = new DensityTools::BulkedWaterFilter(filter_func, pad, zmin, zmax);
          }
          
          
          return(true);
        }

        std::string print() const {
          std::ostringstream oss;

          oss << boost::format("water='%s', prot='%s', pad=%f, bulked='%s', radius=%f, contacts=%d, zrange='%s', grid='%s', mode='%s'")
            % water_string
            % prot_string
            % pad
            % (bulked_spec.empty() ? "undef" : bulked_spec)
            % radius
	    % contacts
            % (zrange_spec.empty() ? "undef" : zrange_spec)
            % (grid_name.empty() ? "undef" : grid_name)
            % filter_mode;

          return(oss.str());
        }
        

        ~BasicWater() {
          if (filter_func)
            delete filter_func;
          filter_func = 0;
        }

        //! Parameters sent to various decorators
        double zmin, zmax;
        //! Extra padding for water
        double pad;
        //! Optional parameter used in by the WaterFilter
        double radius;
	uint contacts;
        //! User-specified strings
        std::string water_string, prot_string, grid_name, filter_mode;
        //! User-specified strings (used internally by the WaterFilter decorators)
        std::string bulked_spec, zrange_spec;
        //! Filter for determining internal waters
        DensityTools::WaterFilterBase* filter_func;
        //! Grid mask for internal waters
        DensityTools::DensityGrid<int> the_grid;

      private:

        bool validateGridMask(const DensityTools::DensityGrid<int>& grid) {
          for (long i=0; i<grid.size(); ++i)
            if (grid(i) < 0 || grid(i) > max_grid_mask_value)
              return(false);

          return(true);
        }

      };
      

  };

};


#endif
