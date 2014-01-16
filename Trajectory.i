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

%shared_ptr(loos::Trajectory)

%header %{
#include <Trajectory.hpp>
#include <StreamWrapper.hpp>
%}



namespace loos {

  // Note, this is actually derived from boost::noncopyable
  class Trajectory;
  typedef boost::shared_ptr<Trajectory> pTraj;
  class Trajectory {
  public:
    Trajectory();
    Trajectory(const std::string& s);
    Trajectory(const char* s);
    virtual ~Trajectory();

    virtual uint natoms(void) const =0;
    virtual float timestep(void) const =0;
    virtual uint nframes(void) const =0;
    void rewind(void);
    virtual bool hasPeriodicBox(void) const =0;
    virtual GCoord periodicBox(void) const =0;
    virtual std::vector<GCoord> coords(void) =0;
    virtual void updateGroupCoords(AtomicGroup& g) =0;
    void seekNextFrame(void);
    void seekFrame(const uint i);
    virtual bool parseFrame(void) =0;
    bool readFrame(void);
    bool readFrame(const int i);
  };


  %extend Trajectory {

    ulong __len__() const {
      return($self->nframes());
    }


  };


  %pythoncode %{
    
class TrajectoryIterator:
    """
    This class wraps a LOOS Trajectory and an AtomicGroup so
    that the trajectory can be used as a Python iterator.  A
    range of frames to iterate over can be set, as well as a
    stride (step).  Any AtomicGroup that shares atoms with
    the one bound into the Iterator will be updated at each
    iteration.

    Basic usage:
      model = loos.createSystem(model_name)
      calphas = selectAtoms(model, 'name == "CA"')
      traj = loos.createTrajectory(traj_name, model)

      # Take every tenth frame...
      itraj = TrajectoryIterator(traj, model, stride = 10)
      for frame in itraj:
         computeSomething(frame)
         computeSomethingElse(calphas)

      # Take every other frame, skipping the first 100
      itraj = TrajectoryIterator(traj, model, stride = 2, skip = 100)
      for frame in itraj:
         computeSomething(frame)
         computeSomethingElse(calphas)

    """
    def __init__(self, traj, frame, skip = 0, stride = 1, iterator = None):
        self.traj = traj
        self.frame = frame
        self.index = None  

        if (iterator is None):
            self.iterator = iter(range(skip, traj.nframes(), stride))
        else: 
            self.iterator = iter(iterator)


    def __iter__(self):
        return(self)

    def currentIndex(self):
        return(self.index)

    def next(self):
        self.index = next(self.iterator)
        self.traj.readFrame(self.index)
        self.traj.updateGroupCoords(self.frame)
        return(self.frame)






class AlignedTrajectory:
    """
    This class provides an iterator over a trajectory that has
    been iteratively aligned (see loos::iterativeAlignment()
    in the C++ LOOS documentation).


    Basic usage:

      calphas = selectAtoms(model, 'name == "CA"')

      # Align and iterate over same set of atoms
      atraj = AlignedTrajectory(traj, calphas)
      for frame in atraj:
         ...

      # Align using C-alphas but iterate over all atoms
      atraj = AlignedTrajectory(traj, model, alignwith = calphas)
      for frame in atraj:
         ...

      # Align using C-alphas but iterate over all atoms, skipping
      # every other frame and the first 100 frames
      atraj = AlignedTrajectory(traj, model, alignwith = calphas, skip = 100, stride = 2)
      for frame in atraj:
         ...


    """


class AlignedTrajectory:

    def __init__(self, traj, frame, skip = 0, stride = 1, iterator = None, alignwith = None):
        self.frame = frame
        self.traj = traj

        if (alignwith is None):
            alignwith = self.frame
        
        if (iterator is None):
            it = iter(range(skip, traj.nframes(), stride))
        else:
            it = iter(iterator)

        self.framelist = UIntVector()
        for i in it:
            self.framelist.push_back(i)
        
        res = iterativeAlignmentPy(alignwith, traj, self.framelist)
        # Must copy out of tuple, otherwise will lose them...
        self.xforms = []
        for x in res.transforms:
            self.xforms.append(x)

        self.index = 0
 


    def __iter__(self):
        return(self)

    def next(self):
        if (self.index >= len(self.framelist)):
            raise StopIteration
        self.traj.readFrame(self.framelist[self.index])
        self.traj.updateGroupCoords(self.frame)
        self.frame.applyTransform(self.xforms[self.index])

        self.index += 1
        return(self.frame)

    def currentIndex(self):
        return(self.framelist[self.index-1])



	      %}
  

}

