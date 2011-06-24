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


#if !defined(GRID_UTILS_HPP)
#define GRID_UTILS_HPP

#include <DensityGrid.hpp>

namespace loos {

  namespace DensityTools {
  
    template<typename T>
    class Threshold {
    public:
      Threshold(const T& t) : thresh(t) { }
      bool operator()(const T& t) const { return(t >= thresh); }
    private:
      T thresh;
    };

    template<typename T>
    class ThresholdRange {
    public:
      ThresholdRange(const T& l, const T& h) : lo(l), hi(h) { }
      bool operator()(const T& t) const { return( t >= lo && t <= hi ); }
    private:
      T lo, hi;
    };


    template<typename T>
    class NonzeroDensity {
    public:
      bool operator()(const T& t) const { return(t > 0.0); }
    };



    template<typename T, class Functor>
    std::vector<DensityGridpoint> floodFill(const DensityGridpoint seed, const DensityGrid<T>& data_grid,
                                      const int id, DensityGrid<int>& blob_grid, const Functor& op)
    {
      std::vector<DensityGridpoint> stack;
      stack.push_back(seed);
      std::vector<DensityGridpoint> list;
      list.push_back(seed);
      blob_grid(seed) = id;

      while(!stack.empty()) {
        DensityGridpoint point = stack.back();
        stack.pop_back();
        for (int k=-1; k<=1; ++k)
          for (int j=-1; j<=1; ++j)
            for (int i=-1; i<=1; ++i) {
              if (i == j && j == k)
                continue;
              DensityGridpoint probe = point + DensityGridpoint(i, j, k);
              if (!data_grid.inRange(probe))
                continue;
              if (blob_grid(probe) == 0 && op(data_grid(probe))) {
                blob_grid(probe) = id;
                stack.push_back(probe);
                list.push_back(probe);
              }
            }
      
      }

      return(list);
    }

    template<typename T, class Functor>
    int floodFill(const DensityGridpoint seed, const DensityGrid<T>& data_grid, const Functor& op) {
      DensityGrid<int> blob_grid(data_grid.minCoord(), data_grid.maxCoord(), data_grid.gridDims());
    
      std::vector<DensityGridpoint> result = floodFill(seed, data_grid, 1, blob_grid, op);
      return(result.size());
    }


    template<typename T, class Functor>
    std::vector<loos::GCoord> findPeaks(const DensityGrid<T>& grid, DensityGrid<int>& blobs, const Functor& op) {
      std::vector<loos::GCoord> peaks;

      DensityGridpoint dims = grid.gridDims();

      int id = 0;
      for (int k=0; k<dims.z(); ++k)
        for (int j=0; j<dims.y(); ++j)
          for (int i=0; i<dims.x(); ++i) {
            DensityGridpoint p (i, j, k);
            if (!blobs(p) && op(grid(p))) {
              std::vector<DensityGridpoint> points = floodFill(p, grid, ++id, blobs, op);
              if (!points.empty()) {
                loos::GCoord center(0,0,0);
                double mass = 0.0;
                for (std::vector<DensityGridpoint>::iterator i = points.begin(); i != points.end(); ++i) {
                  double m = grid(*i);
                  center += m * grid.gridToWorld(*i);
                  mass += m;
                }
                center /= mass;
                peaks.push_back(center);
              }
            }
          }
    
      return(peaks);
    }



    template<typename T, class Functor>
    std::vector<loos::GCoord> findPeaks(const DensityGrid<T>& grid, const Functor& op) {

      DensityGridpoint dims = grid.gridDims();
      DensityGrid<int> blobs(grid.minCoord(), grid.maxCoord(), dims);
      return(findPeaks(grid, blobs, op));
    }


    template<class T, class Functor>
    loos::AtomicGroup gridToAtomicGroup(const DensityGrid<T>& grid, const Functor& op) {
      loos::AtomicGroup group;
      DensityGridpoint dims = grid.gridDims();

      int id = 0;
      for (int k=0; k<dims.z(); ++k)
        for (int j=0; j<dims.y(); ++j)
          for (int i=0; i<dims.x(); ++i) {
            DensityGridpoint p(i, j, k);
            if (op(grid(p))) {
              loos::pAtom atom(new loos::Atom(++id, "UNK", grid.gridToWorld(p)));
              atom->resid(id);
              atom->resname("GRD");
              atom->mass(grid(p));
              group.append(atom);
            }
          }
      return(group);
    }



    template<class T>
    void gridConvolve(DensityGrid<T>& grid, DensityGrid<T>& kernel) {
      DensityGrid<T> tmp(grid);
      DensityGridpoint gdim = grid.gridDims();
      DensityGridpoint kdim = kernel.gridDims();
    
      int kkc = kdim.z() / 2;
      int kjc = kdim.y() / 2;
      int kic = kdim.x() / 2;

      for (int k=0; k<gdim.z(); ++k)
        for (int j=0; j<gdim.y(); ++j)
          for (int i=0; i<gdim.x(); ++i) {
            T sum = 0;

            for (int kk=0; kk<kdim.z(); ++kk)
              for (int jj=0; jj<kdim.y(); ++jj)
                for (int ii=0; ii<kdim.x(); ++ii) {
                  int gk = k + (kk - kkc);
                  if (gk < 0 || gk > gdim.z())
                    continue;

                  int gj = j + (jj - kjc);
                  if (gj < 0 || gj > gdim.y())
                    continue;

                  int gi = i + (ii - kic);
                  if (gi < 0 || gi > gdim.x())
                    continue;

                  sum += grid(k,j,i) * kernel(kk, jj, ii);
                }

            tmp(k, j, i) = sum;
          }

      grid = tmp;
    }


    template<class T>
    void gridConvolve(DensityGrid<T>& grid, std::vector<T>& kernel) {
      DensityGridpoint gdim = grid.gridDims();
      int kn = kernel.size();
      int kc = kn / 2;

      DensityGrid<T> tmp(grid.minCoord(), grid.maxCoord(), grid.gridDims());
      tmp.metadata(grid.metadata());

      // First convolve along the k axis...
      for (int j=0; j<gdim.y(); ++j)
        for (int i=0; i<gdim.x(); ++i)
          for (int k = 0; k<gdim.z(); ++k) {
            T sum = 0;
            for (int ii = 0; ii<kn; ++ii) {
              int idx = k + ii - kc;
              if (idx < 0 || idx >= gdim.z())
                continue;
              sum += grid(idx, j, i) * kernel[ii];
            }
            tmp(k, j, i) = sum;
          }

      // Convolve along the j axis
      DensityGrid<T> tmp2(grid.minCoord(), grid.maxCoord(), grid.gridDims());
      tmp2.metadata(grid.metadata());
      for (int k=0; k<gdim.z(); ++k)
        for (int i=0; i<gdim.x(); ++i)
          for (int j=0; j<gdim.y(); ++j) {
            T sum = 0;
            for (int ii=0; ii<kn; ++ii) {
              int idx = j + ii - kc;
              if (idx < 0 || idx >= gdim.y())
                continue;
              sum += tmp(k, idx, i) * kernel[ii];
            }
            tmp2(k, j, i) = sum;
          }

      // Finally, convolve along the i axis
      for (int k=0; k<gdim.z(); ++k)
        for (int j=0; j<gdim.y(); ++j)
          for (int i=0; i<gdim.x(); ++i) {
            T sum = 0;
            for (int ii=0; ii<kn; ++ii) {
              int idx = i + ii - kc;
              if (idx < 0 || idx >= gdim.x())
                continue;
              sum += tmp2(k, j, idx) * kernel[ii];
            }
            tmp(k, j, i) = sum;
          }

      grid = tmp;
    }

    std::vector<double> gaussian1d(const int, const double);

  };

};



#endif
