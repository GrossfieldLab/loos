/*
  Utility functions/classes for LOOS grids
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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


#include <GridUtils.hpp>


using namespace std;



namespace loos {
  namespace DensityTools {

    
    std::vector<double> gaussian1d(const int w, const double sigma) {
    
      double a = 1.0 / sqrt(2.0 * M_PI * sigma);
      double b = -1.0 / (2.0 * sigma);
    
      std::vector<double> kernel;
      for (int i=0; i<=w; ++i) {
        double x = 2.0*i/w-1.0;
        kernel.push_back(a*exp(b*x*x));
      }
    
      return(kernel);
    }


  };
};



