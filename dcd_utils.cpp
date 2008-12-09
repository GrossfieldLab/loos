/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
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





#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pwd.h>

#include <string>
#include <sstream>

#include <utils.hpp>
#include <dcd_utils.hpp>


namespace loos {

  //! Map DCD coordinates onto a grid given a window and a range.
  /*!
    Returns a raw double array representing a 3D grid of densities
    averaging over the specified window.  Each coordinate in a DCD frame
    is considered a point of unit mass, for purposes of computing
    densities...

    \param dcd The DCD object to use for reading frames
    \param avg_box Stores the average box size (extracted from crystal params)
    \param avg_unitvol Stores the average box-volume
    \param gridsizes Array giving the i, j, k dimensions of the grid.
    \param indices Indices into the DCD frame to operate over
    (i.e. which atoms)
    \param frameno Which frame to start on
    \param window How many frames to operate over
    \param scale Scales up the density

  */
  double *gridify(DCD& dcd, double *avg_box, double *avg_unitvol, int gridsizes[], const std::vector<int> indices, int frameno = 0, int window = 1, double scale = 1.0) {
    int grid_dim = gridsizes[0] * gridsizes[1] * gridsizes[2];
    int* grid = new int[grid_dim];
    double* density = new double[grid_dim];
    int i;
    for (i=0; i<grid_dim; i++)
      density[i] = 0.;

    avg_box[0] = avg_box[1] = avg_box[2] = 0.0;
    *avg_unitvol = 0.0;

    // Iterate over all frames from the initially specified one through
    // the window... 
    int offset;
    for (offset = 0; offset < window; offset++) {

      for (i=0; i<grid_dim; i++)
        grid[i] = 0;

      bool b = dcd.readFrame(frameno + offset);

      if (!b) {
        return(0);
      }

      // Since the box-size may changed, and we're expecting periodic
      // boundary conditiones, extract the crystal params and then
      // modulo them back into the unit cell... 

      std::vector<dcd_double> xtal = dcd.crystalParams();
      GCoord box(xtal[0], xtal[1], xtal[2]);

      std::vector<GCoord> crds = dcd.mappedCoords(indices);

      for (i=0; i<3; i++)
        avg_box[i] += box[i];

      double delta[3];
      for (i=0; i<3; i++)
        delta[i] = (double)gridsizes[i] / box[i];

      // Iterate over all coordinates binning them into the grid...
      int n = crds.size();
      for (i=0; i<n; i++) {
        GCoord cc = crds[i];
        cc.reimage(box);
        double x;
        int u[3];
        int j;
        for (j=0; j<3; j++) {
          x = cc[j] + box[j] / 2.0;
          u[j] = (int)floor(x * delta[j]);
          assert (u[j] >= 0 && u[j] < gridsizes[j]);   // Safety check
        }

        // Get the grid indices
        int ii = (u[2] * gridsizes[1] + u[1]) * gridsizes[0] + u[0];
        assert(ii < grid_dim);    // Again with the safety check
        grid[ii] += 1;
      }

    
      // Now convert into density values and accumulate...
      double unit_volume = 1.0;
      for (i=0; i<3; i++)
        unit_volume *= (box[i] / gridsizes[i]);
      *avg_unitvol += unit_volume;

      for (i=0; i<grid_dim; i++)
        density[i] += (scale * grid[i]) / unit_volume;

    }

    // Convert to averages...

    for (i=0; i<grid_dim; i++)
      density[i] /= window;
    for (i=0; i<3; i++)
      avg_box[i] /= window;
    *avg_unitvol /= window;

    return(density);
  }



}
