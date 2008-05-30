/*
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Ensemble calculations...
*/

#if !defined(ENSEMBLES_HPP)
#define ENSEMBLES_HPP

namespace loos {
  //! Compute the average structure of a set of AtomicGroup objects
  /**Is xform-aware */
  AtomicGroup averageStructure(const vector<AtomicGroup>&);

  //! Compute an iterative superposition (a la Alan)
  /**Is xform-aware */
  greal iterativeAlignment(vector<AtomicGroup>& ensmble, greal threshold, int maxiter=1000);
};



#endif
