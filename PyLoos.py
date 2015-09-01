import loos    


class PyTraj:
    """
    This class wraps a LOOS Trajectory and an AtomicGroup so
    that the trajectory can be used as a Python iterator.  A
    range of frames to iterate over can be set, as well as a
    stride (step).  Any AtomicGroup that shares atoms with
    the one bound into the Iterator will be updated at each
    iteration.


    """

    def __init__(self, fname, model, **dict):

        skip = 0
        stride = 1
        iterator = None
        
        if 'skip' in dict:
            skip = dict['skip']
        if 'stride' in dict:
            stride = dict['stride']
        if 'iterator' in dict:
            iterator = dict['iterator']
        
        self.frame = model
        self.fname = fname
        self.traj = loos.createTrajectory(fname, model)

        self.framelist = []
        if iterator is None:
            it = range(skip, self.traj.nframes(), stride)
        else:
            it = iter(iterator)
        
        for i in it:
            self.framelist.append(i)

        self.index = 0


    def __iter__(self):
        return(self)

    def __len__(self):
        return(len(self.framelist))

    def reset(self):
        self.index = 0

    def next(self):
        if (self.index >= len(self.framelist)):
            raise StopIteration
        frame = self.__getitem__(self.index)

        self.index += 1
        return(self.frame)


    def readFrame(self, i):
        if (i < 0 or i >= len(self.framelist)):
            raise IndexError
        self.traj.readFrame(i)
        return(self.frame)

    def currentFrame(self):
        return(self.frame)
    
    def currentIndexInTrajectory(self):
        return(self.framelist[self.index-1])

    def currentIndexInFramelist(self):
        return(self.index-1)

    def averageStructure(self):
        flist = loos.UIntVector(self.framelist)
        return(loos.averageStructure(self.frame, self.traj, flist))



    def getSlice(self, s):
        indices = list(range(*s.indices(self.__len__())))
        ensemble = []
        for i in indices:
            self.traj.readFrame(self.framelist[i])
            self.traj.updateGroupCoords(self.frame)
            dup = self.frame.copy()
            ensemble.append(dup)
        return(ensemble)


    def __getitem__(self, i):

        if isinstance(i, slice):
            return(self.getSlice(i))

        if (i < 0):
            i += len(self.framelist)
        if (i >= len(self.framelist) or i < 0):
            raise IndexError
        self.traj.readFrame(self.framelist[i])
        self.traj.updateGroupCoords(self.frame)
        return(self.frame)






class PyAlignedTraj:
    """
    This class provides an iterator over a trajectory that has
    been iteratively aligned (see loos::iterativeAlignment()
    in the C++ LOOS documentation).


    Basic usage:

      calphas = selectAtoms(model, 'name == "CA"')

      # Align and iterate over same set of atoms
      atraj = PyAlignedTraj(traj, calphas)
      for frame in atraj:
         ...

      # Align using C-alphas but iterate over all atoms
      atraj = PyAlignedTraj(traj, model, alignwith = calphas)
      for frame in atraj:
         ...

      # Align using C-alphas but iterate over all atoms, skipping
      # every other frame and the first 100 frames
      atraj = PyAlignedTraj(traj, model, alignwith = calphas, skip = 100, stride = 2)
      for frame in atraj:
         ...


    """
class PyAlignedTraj(PyTraj):

    def __init__(self, traj, frame, skip = 0, stride = 1, iterator = None, alignwith = None):
        PyTraj.__init__(self, traj, frame, skip, stride, iterator)

        if (alignwith is None):
            alignwith = self.frame

        res = loos.iterativeAlignmentPy(alignwith, traj, self.framelist)
        self.xforms = loos.XFormVector()
        for x in res.transforms:
            self.xforms.push_back(x)


    def getSlice(self, s):
        indices = list(range(*s.indices(self.__len__())))
        ensemble = []
        for i in indices:
            self.traj.readFrame(self.framelist[i])
            self.traj.updateGroupCoords(self.frame)
            dup = self.frame.copy()
            dup.applyTransform(self.xforms[i])
            ensemble.append(dup)
        return(ensemble)


    def __getitem__(self, i):
        if isinstance(i, slice):
            return(self.getSlice(i))

        f = PyTraj.__getitem__(self, i)
        f.applyTransform(self.xforms[i])
        return(f)

    def transform(self, i):
        if (i < 0):
            i += len(self.framelist)
        if (i >= len(self.framelist) or i < 0):
            raise IndexError
        return(self.xforms[i])


    def averageStructure(self):
        return(loos.averageStructure(self.frame, self.xforms, self.traj, self.framelist))

    def currentTransform(self):
        return(self.xforms[self.index-1])



class VirtualTrajectory:

    def __init__(self, *trajs, **dict):
        self.skip = 0
        self.stride = 1
        self.nframes = 0
        self.iterator = None
        self.trajectories = list(trajs)

        self.index = 0
        self.framelist = []
        self.trajlist = [] 
        self.stale = 1

        if 'skip' in dict:
            self.skip = dict['skip']
        if 'stride' in dict:
            self.stride = dict['stride']
        if 'iterator' in dict:
            self.iterator = dict['iterator']
        
    def addTrajectory(self, traj):
        self.trajectories.append(traj)
        self.stale = 1

    def countTrajectoryFrames(self):
        n = 0
        for t in self.trajectories:
            n += len(t)
        return(n)

    def verifyModels(self):
        if not self.trajectories:
            return
        n = 0
        for t in self.trajectories:
            if n == 0:
                n = len(t.frame)
            elif n != len(t.frame):
                raise RuntimeError('Inconsistant models or subsets inside a virtual trajectory')

    def initFrameList(self):
        n = self.countTrajectoryFrames()
        if (self.iterator is None):
            it = iter(range(self.skip, n, self.stride))
        else:
            it = iter(iterator)

        frames = []
        trajs = []
        for t in self.trajectories:
            for i in range(len(t)):
                frames.append(i)
                trajs.append(t)

        self.framelist = frames
        self.trajlist = trajs
            
    def __len__(self):
        if self.stale:
            self.initFrameList()
        return(len(self.framelist))

                
