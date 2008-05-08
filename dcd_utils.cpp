/*
  dcd_utils.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Util functions for DCD manipulation...
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


double *gridify(DCD& dcd, double *avg_box, double *avg_unitvol, int gridsizes[], const vector<int> indices, int frameno = 0, int window = 1, double scale = 1.0) {
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

    vector<dcd_double> xtal = dcd.crystalParams();
    GCoord box(xtal[0], xtal[1], xtal[2]);

    vector<GCoord> crds = dcd.mappedCoords(indices);

    for (i=0; i<3; i++)
      avg_box[i] += box[i];

    double delta[3];
    for (i=0; i<3; i++)
      delta[i] = (double)gridsizes[i] / box[i];

    // Iterate over all coordinates binning them into the grid...
    int n = crds.size();
    for (i=0; i<n; i++) {
      GCoord cc = crds[i];
      cc.canonical(box);
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



