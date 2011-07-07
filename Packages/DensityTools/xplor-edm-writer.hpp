// -------------------------------------------------
// ASCII Xplor-formatted Electron Density Map writer
// -------------------------------------------------

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008 Tod D. Romo, Alan Grossfield
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




#if !defined(XPLOREDMWRITER_HPP)
#define XPLOREDMWRITER_HPP

#include <loos.hpp>

#include <DensityGrid.hpp>
#include <SimpleMeta.hpp>

namespace loos {

  namespace DensityTools {

    //! Functor for writing out ASCII formatted X-Plor electron density maps
    template<class T>
    struct XEDMWriter {
      XEDMWriter(std::ostream& os) : i(0), mos(os), fmt(5)
      {
        fmt.scientific().width(12).right();
      }

      void operator()(const T d) {
        mos << fmt(d);
        if (i++ >= 5) {
          mos << std::endl;
          i = 0;
        }
      }

      //! Special m.f. for starting a new frame
      void frame(const int k) {
        if (i != 0) {
          i = 0;
          mos << std::endl;
        }
        mos << std::setw(8) << k << std::endl;
      }


      int i;
      std::ostream& mos;
      loos::Fmt fmt;
    };



    //! Write out an DensityTools::DensityGrid as an ASCII formatted X-PLOR electron density map
    template<class T> void writeXplorEDM(std::ostream& os, DensityTools::DensityGrid<T>& grid) {
      loos::GCoord gridmin = grid.minCoord();
      loos::GCoord gridmax = grid.maxCoord();
      loos::GCoord delta = grid.gridDelta();
      loos::GCoord gridsize;
      DensityTools::DensityGridpoint dims = grid.gridDims();

      DensityTools::DensityGridpoint mins, maxs, nas;

      // Handle the header...
      for (int i=0; i<3; i++) {
        mins[i] = static_cast<int>(floor(gridmin[i] * delta[i]));
        maxs[i] = static_cast<int>(floor(gridmax[i] * delta[i]));
        gridsize[i] = dims[i] / delta[i];
      }
      nas = dims;

      // Special handling for grid meta-data...
      SimpleMeta meta = grid.metadata();
      os << std::endl << std::setw(8) << meta.size() << std::endl;
      for (SimpleMeta::iterator i = meta.begin(); i != meta.end(); ++i)
        os << *i << std::endl;

      for (int i=0; i<3; i++)
        os << std::setw(8) << nas[i] << std::setw(8) << mins[i] << std::setw(8) << maxs[i];
      os << std::endl;

      loos::Fmt fc(5);
      fc.width(12).scientific();

      // Assume our "crystal" is orthonormal...
      os << fc(gridsize[0]) << fc(gridsize[1]) << fc(gridsize[2]) << fc(90.0) << fc(90.0) << fc(90.0) << std::endl;
      os << "ZYX\n";

      // Instantiate the writing functor...
      XEDMWriter<T> writer(os);

      // The format writes out the map a plane at a time, so we extract
      // a plane via operator[] and operate on that...

      for (int k=0; k<dims[2]; k++) {
        DensityTools::DensityGridPlane<T> plane = grid[k];

        // Prime the output
        writer.frame(k);

        for (int j=0; j<dims[1]; j++)
          for (int i=0; i<dims[0]; i++)
            writer(plane[j][i]);
      }

      os << std::endl << std::endl;
    }

  };

};


#endif
