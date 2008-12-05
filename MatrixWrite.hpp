/*
  MatrixWrite.hpp

  Class for handling writing of Matrix objects
*/


/*
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester
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
  template<class T, class P, template<typename> class S>
  struct MatrixWriteImpl;

  // The following are the templated global functions.  Do not
  // overload/specialize them.  Instead, specialize the
  // MatrixWriteImpl class...

  //! Write a submatrix to a stream
  template<class T, class P, template<typename> class S>
  std::ostream& writeAsciiMatrix(std::ostream& os, const Matrix<T,P,S>& M,
                                 const std::string& meta, const MDuple& start,
                                 const MDuple& end, const bool trans = false) {
    return(MatrixWriteImpl<T,P,S>::write(os, M, meta, start, end, trans));
  }

  //! Write an entire matrix to a stream
  template<class T, class P, template<typename> class S>
  std::ostream& writeAsciiMatrix(std::ostream& os, const Matrix<T,P,S>& M,
                                 const std::string& meta, const bool trans = false) {
    MDuple start(0,0);
    MDuple end(M.rows(), M.cols());
    return(MatrixWriteImpl<T,P,S>::write(os, M, meta, start, end, trans));
  }

  //! Write a submatrix to a file
  template<class T, class P, template<typename> class S>
  void writeAsciiMatrix(const std::string& fname, const Matrix<T,P,S>& M,
                                 const std::string& meta, const MDuple& start,
                                 const MDuple& end, const bool trans = false) {
    std::ofstream ofs(fname.c_str());
    if (ofs == 0)
      throw(std::runtime_error("Cannot open " + fname + " for writing."));
    MatrixWriteImpl<T,P,S>::write(ofs, M, meta, start, end, trans);
  }


  //! Write an entire matrix to a file
  template<class T, class P, template<typename> class S>
  void writeAsciiMatrix(const std::string& fname, const Matrix<T,P,S>& M,
                                 const std::string& meta, const bool trans = false) {
    MDuple start(0,0);
    MDuple end(M.rows(), M.cols());

    std::ofstream ofs(fname.c_str());
    if (ofs == 0)
      throw(std::runtime_error("Cannot open " + fname + " for writing."));
    MatrixWriteImpl<T,P,S>::write(ofs, M, meta, start, end, trans);
  }


  // Writing implementation and specializations...

  template<class T, class P, template<typename> class S>
  struct MatrixWriteImpl {
    static std::ostream& write(std::ostream& os,
                               const Matrix<T,P,S>& M,
                               const std::string& meta,
                               const MDuple& start, const MDuple& end,
                               const bool trans) {
      os << "# " << meta << std::endl;
      os << boost::format("# %d %d (%d)\n") % M.rows() % M.cols() % trans;
      for (int j=start.first; j<end.first; j++) {
        for (int i=start.second; i<end.second; i++)
          if (trans)
            os << M(i, j) << " ";
          else
            os << M(j, i) << " ";
        os << std::endl;
      }
      return(os);
    }
  };

  //! Write out a sparse matrix.
  /** Ignores \a start, \a end, and \a trans */
  template<class T, class P>
  struct MatrixWriteImpl<T, P, SparseArray > {
    static std::ostream& write(std::ostream& os,
                               const Matrix<T,P,SparseArray>& M,
                               const std::string& meta,
                               const MDuple& start, const MDuple& end,
                               const bool trans) {
      os << "# " << meta << std::endl;
      os << boost::format("# %d %d %d SPARSE\n") % M.actualSize() % M.rows() % M.cols();
      typename Matrix<T,P,SparseArray>::const_iterator ci;
      
      for (ci = M.begin(); ci != M.end(); ++ci)
        os << (*ci).first << "\t" << (*ci).second << std::endl;
      
      return(os);
    }
  };

  //! Write out a triangular matrix
  /** Ignores \a start, \a end, and \a trans */
  template<class T, template<typename> class S>
  struct MatrixWriteImpl<T, Triangular, S> {
    static std::ostream& write(std::ostream& os,
                               const Matrix<T,Triangular,S>& M,
                               const std::string& meta,
                               const MDuple& start, const MDuple& end,
                               const bool trans) {
      os << "# " << meta << std::endl;
      os << boost::format("# %d TRIANGULAR\n") % M.rows();
      long s = M.size();
      for (long i=0; i<s; i++)
        os << M[i] << std::endl;
      
      return(os);
    }
  };

}


#endif
