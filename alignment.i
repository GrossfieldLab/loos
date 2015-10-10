/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2015, Tod D. Romo, Alan Grossfield
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

%rename(cpp_iterativeAlignment)    loos::iterativeAlignment;

%header %{
#include <alignment.hpp>





#if defined(SWIGPYTHON)

  namespace loos {
    struct AlignmentResult {
      std::vector<XForm> transforms;
      double rmsd;
      int iterations;
    };


    AlignmentResult iterativeAlignmentPy(std::vector< std::vector<double> >& ensemble, greal threshold = 1e-6, int maxiter=1000) {
      AlignmentResult res;

      boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble, threshold, maxiter);
      res.transforms = boost::get<0>(ares);
      res.rmsd = boost::get<1>(ares);
      res.iterations = boost::get<2>(ares);

      return(res);
    }



    AlignmentResult iterativeAlignmentPy(std::vector<AtomicGroup>& ensemble, greal threshold = 1e-6, int maxiter=1000) {
      AlignmentResult res;

      boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(ensemble, threshold, maxiter);
      res.transforms = boost::get<0>(ares);
      res.rmsd = boost::get<1>(ares);
      res.iterations = boost::get<2>(ares);

      return(res);
    }

    AlignmentResult iterativeAlignmentPy(const AtomicGroup& g, pTraj& traj, std::vector<uint>& frame_indices, greal threshold = 1e-6, int maxiter=1000) {
      AlignmentResult res;

      boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(g, traj, frame_indices, threshold, maxiter);
      res.transforms = boost::get<0>(ares);
      res.rmsd = boost::get<1>(ares);
      res.iterations = boost::get<2>(ares);

      return(res);
    }


    AlignmentResult iterativeAlignmentPy(const AtomicGroup& g, pTraj& traj, greal threshold = 1e-6, int maxiter=1000) {
      AlignmentResult res;

      boost::tuple<std::vector<XForm>, greal, int> ares = iterativeAlignment(g, traj, threshold, maxiter);
      res.transforms = boost::get<0>(ares);
      res.rmsd = boost::get<1>(ares);
      res.iterations = boost::get<2>(ares);

      return(res);
    }


  }

#endif

%}

%include "alignment.hpp"


namespace loos {

#if defined(SWIGPYTHON)

  struct AlignmentResult {
    std::vector<XForm> transforms;
    double rmsd;
    int iterations;
  };


  loos::AlignmentResult iterativeAlignmentPy(std::vector< std::vector<double> >& ensemble, greal threshold = 1e-6, int maxiter = 1000);
  
  loos::AlignmentResult iterativeAlignmentPy(std::vector<AtomicGroup>& ensemble, greal threshold = 1e-6, int maxiter = 1000);


  loos::AlignmentResult iterativeAlignmentPy(const AtomicGroup& g, pTraj& traj, std::vector<uint>& frames, greal threshold = 1e-6, int maxiter = 1000);

  loos::AlignmentResult iterativeAlignmentPy(const AtomicGroup& g, pTraj& traj, greal threshold = 1e-6, int maxiter = 1000);

};


%pythoncode %{
def xformVectorToList(v):
    l = []
    for x in v:
        l.append(XForm(x))
    return(l)


def iterativeAlignEnsemble(ensemble, threshold=1e-8, maxiter=1000):
    # Convert to vector<AtomicGroup>...smart pointers should make overhead ok...
    enlist = loos.DoubleVector()
    for e in ensemble:
       enlist.push_back(e.coordsAsVector())
    result = iterativeAlignmentPy(enlist, threshold, maxiter)
    return(xformVectorToList(result.transforms), result.rmsd, result.iterations)



# Optional 'framelist' argument specifies indices of frames to use
def iterativeAlignment(model, traj, threshold=1e-8, maxiter=1000, **kwargs):
    # If traj is not a loos.Trajectory, assume it supports the trajectory()
    # method to access the underlying loos one
    if not isinstance(traj, loos.Trajectory):
        traj = traj.trajectory()

    # Handle framelist request
    if 'framelist' in kwargs:
        framelist = kwargs['framelist']
        # If it's not a vector<uint>, assume it's iterable
        if not isinstance(framelist, loos.UIntVector):
           flist = loos.UIntVector(kwargs['framelist'])
        else:
           flist = framelist
        result = iterativeAlignmentPy(model, traj, flist, threshold, maxiter)

    else:
        result = iterativeAlignmentPy(model, traj, threshold, maxiter)

    return(xformVectorToList(result.transforms), result.rmsd, result.iterations)

%}

#endif

%template(XFormVector) std::vector<loos::XForm>;
