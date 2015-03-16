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



#if !defined(LOOS_SORTING_HPP)
#define LOOS_SORTING_HPP


#include <vector>
#include <algorithm>
#include <stdexcept>

#include <loos/loos_defs.hpp>

namespace loos {

  //! Policy class for sorting in ascending sequence
  template<typename T>
  class AscendingSort {
  public:
    AscendingSort(const T& A) : A_(A) { }
    bool operator()(const uint i, const uint j) const {
      return(A_[i] < A_[j]);
    }
    
  private:
    const T& A_;
  };

  //! Policy class for sorting in descending sequence
  template<typename T>
  class DescendingSort {
  public:
    DescendingSort(const T& A) : A_(A) { }
    bool operator()(const uint i, const uint j) const {
      return(A_[i] > A_[j]);
    }
    
  private:
    const T& A_;
  };


  //! Sort a container using the given sort policy, returning the
  //! indices that permutes the container into the sorted order.
  /**
   * The container to be sorted must support T::size() and
   * T::operator[].  What is returned is a vector of unsigned ints
   * that represent the index into the container when it is sorted.
   */
  template<typename T, class SortPredicate>
  std::vector<uint> sortedIndex(const T& A) {
    std::vector<uint> indices(A.size());
    
    for (uint i = 0; i<A.size(); ++i)
      indices[i] = i;
    
    std::sort(indices.begin(), indices.end(), SortPredicate(A));
    return(indices);
  }

  //! Sort a container in ascending sequence
  template<typename T>
  std::vector<uint> sortedIndex(const T& A) {
    return(sortedIndex<T, AscendingSort<T> >(A));
  }

}



#endif
