/*
  MatrixWrite.hpp

  Class for handling writing of Matrix objects
*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008-2009, Tod D. Romo, Alan Grossfield
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


#if !defined(MATRIXWRITE_HPP)
#define MATRIXWRITE_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>   // ?
#include <stdexcept>
#include <cassert>
#include <iterator>

#include <utility>
#include <boost/format.hpp>

#include <Matrix.hpp>


namespace loos {

  // Convenience typdef for specifying end-points of sub-matrices
  typedef std::pair<int, int> MDuple;

  // Forward declaration for matrix writing implementation
  template<class T, class P, template<typename> class S, class F>
  struct MatrixWriteImpl;

  namespace internal {
    template<typename T>
    struct BasicMatrixFormatter {
      std::string operator()(const T& t) {
        std::stringstream ss;
        ss << t;
        return(ss.str());
      }
    };
  }

  // The following are the templated global functions.  Do not
  // overload/specialize them.  Instead, specialize the
  // MatrixWriteImpl class...

  //! Write a submatrix to a stream
  /**
   * This family of functions write a matrix in ASCII format suitable
   * for loading into Octave/Matlab or gnuplot.  The \a meta
   * information is written as part of the comment at the start of the
   * file.  The MDuple \a start and \a end are just pairs that give an
   * \a (j,i) starting and ending point within the matrix to write.
   * Note that these arguments are not always honored (such as with a
   * triangular matrix).  The \a trans flag causes the output matrix
   * to be the transpose of the stored matrix.  The \a fmt arg is a
   * functor that is expected for format each element of the matrix as
   * a string.  You can use this to adjust the precision of the output
   * or delimit it, etc.  The default is to use whatever the default
   * operator<<() would be for type T.
   */
  template<class T, class P, template<typename> class S, class F>
  std::ostream& writeAsciiMatrix(std::ostream& os, const Math::Matrix<T,P,S>& M,
                                 const std::string& meta, const MDuple& start,
                                 const MDuple& end, const bool trans = false, F fmt = F()) {
    return(MatrixWriteImpl<T,P,S,F>::write(os, M, meta, start, end, trans, fmt));
  }


  //! Write a submatrix to a stream
  template<class T, class P, template<typename> class S>
  std::ostream& writeAsciiMatrix(std::ostream& os, const Math::Matrix<T,P,S>& M,
                                 const std::string& meta, const MDuple& start,
                                 const MDuple& end, const bool trans = false) {
    return(MatrixWriteImpl<T,P,S,internal::BasicMatrixFormatter<T> >::write(os, M, meta, start, end, trans));
  }

  //! Write an entire matrix to a stream
  template<class T, class P, template<typename> class S, class F>
  std::ostream& writeAsciiMatrix(std::ostream& os, const Math::Matrix<T,P,S>& M,
                                 const std::string& meta, const bool trans = false, F fmt = F()) {
    MDuple start(0,0);
    MDuple end(M.rows(), M.cols());
    return(MatrixWriteImpl<T,P,S,F>::write(os, M, meta, start, end, trans, fmt));
  }

  //! Write an entire matrix to a stream
  template<class T, class P, template<typename> class S>
  std::ostream& writeAsciiMatrix(std::ostream& os, const Math::Matrix<T,P,S>& M,
                                 const std::string& meta, const bool trans = false) {
    MDuple start(0,0);
    MDuple end(M.rows(), M.cols());
    return(MatrixWriteImpl<T,P,S,internal::BasicMatrixFormatter<T> >::write(os, M, meta, start, end, trans));
  }

  //! Write a submatrix to a file
  template<class T, class P, template<typename> class S, class F>
  void writeAsciiMatrix(const std::string& fname, const Math::Matrix<T,P,S>& M,
                                 const std::string& meta, const MDuple& start,
                        const MDuple& end, const bool trans = false, F fmt = F()) {
    std::ofstream ofs(fname.c_str());
    if (ofs == 0)
      throw(std::runtime_error("Cannot open " + fname + " for writing."));
    MatrixWriteImpl<T,P,S,F>::write(ofs, M, meta, start, end, trans, fmt);
  }

  //! Write a submatrix to a file
  template<class T, class P, template<typename> class S>
  void writeAsciiMatrix(const std::string& fname, const Math::Matrix<T,P,S>& M,
                                 const std::string& meta, const MDuple& start,
                        const MDuple& end, const bool trans = false) {
    std::ofstream ofs(fname.c_str());
    if (ofs == 0)
      throw(std::runtime_error("Cannot open " + fname + " for writing."));
    MatrixWriteImpl<T,P,S,internal::BasicMatrixFormatter<T> >::write(ofs, M, meta, start, end, trans);
  }


  //! Write an entire matrix to a file
  template<class T, class P, template<typename> class S, class F>
  void writeAsciiMatrix(const std::string& fname, const Math::Matrix<T,P,S>& M,
                        const std::string& meta, const bool trans = false, F fmt = F()) {
    MDuple start(0,0);
    MDuple end(M.rows(), M.cols());

    std::ofstream ofs(fname.c_str());
    if (ofs == 0)
      throw(std::runtime_error("Cannot open " + fname + " for writing."));
    MatrixWriteImpl<T,P,S,F>::write(ofs, M, meta, start, end, trans, fmt);
  }


  //! Write an entire matrix to a file
  template<class T, class P, template<typename> class S>
  void writeAsciiMatrix(const std::string& fname, const Math::Matrix<T,P,S>& M,
                        const std::string& meta, const bool trans = false) {
    MDuple start(0,0);
    MDuple end(M.rows(), M.cols());

    std::ofstream ofs(fname.c_str());
    if (ofs == 0)
      throw(std::runtime_error("Cannot open " + fname + " for writing."));
    MatrixWriteImpl<T,P,S,internal::BasicMatrixFormatter<T> >::write(ofs, M, meta, start, end, trans);
  }


  // Writing implementation and specializations...

  template<class T, class P, template<typename> class S, class F>
  struct MatrixWriteImpl {
    static std::ostream& write(std::ostream& os,
                               const Math::Matrix<T,P,S>& M,
                               const std::string& meta,
                               const MDuple& start, const MDuple& end,
                               const bool trans, F fmt = F()) {
      os << "# " << meta << std::endl;
      os << boost::format("# %d %d (%d)\n") % M.rows() % M.cols() % trans;
      for (int j=start.first; j<end.first; j++) {
        for (int i=start.second; i<end.second; i++)
          if (trans)
            os << fmt(M(i, j)) << " ";
          else
            os << fmt(M(j, i)) << " ";
        os << std::endl;
      }
      return(os);
    }
  };

  //! Write out a sparse matrix.
  /** Ignores \a start, \a end, and \a trans */
  template<class T, class P, class F>
  struct MatrixWriteImpl<T, P, Math::SparseArray, F> {
    static std::ostream& write(std::ostream& os,
                               const Math::Matrix<T,P,Math::SparseArray>& M,
                               const std::string& meta,
                               const MDuple& start, const MDuple& end,
                               const bool trans, F fmt = F()) {
      os << "# " << meta << std::endl;
      os << boost::format("# %d %d %d SPARSE\n") % M.actualSize() % M.rows() % M.cols();
      typename Math::Matrix<T,P,Math::SparseArray>::const_iterator ci;
      
      for (ci = M.begin(); ci != M.end(); ++ci)
        os << (*ci).first << "\t" << fmt((*ci).second) << std::endl;
      
      return(os);
    }
  };


  //! Write out a triangular matrix
  /** Ignores \a start, \a end, and \a trans */
  template<class T, template<typename> class S, class F>
  struct MatrixWriteImpl<T, Math::Triangular, S, F> {
    static std::ostream& write(std::ostream& os,
                               const Math::Matrix<T,Math::Triangular,S>& M,
                               const std::string& meta,
                               const MDuple& start, const MDuple& end,
                               const bool trans, F fmt = F()) {
      os << "# " << meta << std::endl;
      os << boost::format("# %d TRIANGULAR\n") % M.rows();
      long s = M.size();
      for (long i=0; i<s; i++)
        os << fmt(M[i]) << std::endl;
      
      return(os);
    }
  };

}


#endif
