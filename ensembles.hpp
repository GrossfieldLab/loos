/*
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Ensemble calculations...
*/

#if !defined(ENSEMBLES_HPP)
#define ENSEMBLES_HPP


#include <loos.hpp>
#include <vector>

#include <boost/tuple/tuple.hpp>

namespace loos {
  //! Compute the average structure of a set of AtomicGroup objects
  /**Is xform-aware */
  AtomicGroup averageStructure(const vector<AtomicGroup>& ensemble);

  //! Compute an iterative superposition (a la Alan)
  /**Is xform-aware */
  boost::tuple<vector<XForm>, greal, int> iterativeAlignment(vector<AtomicGroup>& ensemble, greal threshold, int maxiter=1000);
};



#endif
