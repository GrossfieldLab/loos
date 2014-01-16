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

      itraj = TrajectoryIterator(traj, model)
      itraj.setStride(10)    # Use only every 10th frame
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

    def setIterator(self, it):
        self.iterator = iter(it)

    def setRange(self, start, stop, stride = 1):
        self.iterator = iter(xrange(start, stop, stride))

    def setSkip(self, skip): 
        self.iterator = iter(xrange(skip, self.traj.nframes()))  

    def setStride(self, stride): 
        self.iterator = iter(xrange(0, self.traj.nframes(), stride))

    def setFrame(self, frame):
        self.frame = frame

    def currentIndex(self):
        return(self.index)


    def next(self):
        self.index = next(self.iterator)
        self.traj.readFrame(self.index)
        self.traj.updateGroupCoords(self.frame)
        return(self.frame)






class AlignedTrajectory:

    def __init__(self, trajiter, alignwith = None):
        self.frame = trajiter.frame
        self.traj = trajiter.traj
        self.framelist = []
        ensemble = AtomicGroupVector()

        if (alignwith is None):
            alignwith = self.frame

        for f in trajiter:
            subset = alignwith.copy()
            ensemble.push_back(subset)
            self.framelist.append(trajiter.currentIndex())

        res = iterativeAlignmentPy(ensemble)
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

