import loos    


class PyTraj:
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

      # Take every tenth frame...
      itraj = PyTraj(traj, model, stride = 10)
      for frame in itraj:
         computeSomething(frame)
         computeSomethingElse(calphas)

      # Take every other frame, skipping the first 100
      itraj = PyTraj(traj, model, stride = 2, skip = 100)
      for frame in itraj:
         computeSomething(frame)
         computeSomethingElse(calphas)

    """
    def __init__(self, traj, frame, skip = 0, stride = 1, iterator = None):
        self.frame = frame
        self.traj = traj

        if (iterator is None):
            it = iter(range(skip, traj.nframes(), stride))
        else:
            it = iter(iterator)

        self.framelist = loos.UIntVector()
        for i in it:
            self.framelist.push_back(i)

        self.index = 0
 


    def __iter__(self):
        return(self)

    def __len__(self):
        return(len(self.framelist))

    def next(self):
        if (self.index >= len(self.framelist)):
            raise StopIteration
        frame = self.__getitem__(self.index)

        self.index += 1
        return(self.frame)

    def currentIndex(self):
        return(self.framelist[self.index-1])

    def averageStructure(self):
        return(averageStructure(self.frame, self.traj, self.framelist))


    def __getitem__(self, i):
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


    def __getitem__(self, i):
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
        return(averageStructure(self.frame, self.xforms, self.traj, self.framelist))

    def currentTransform(self):
        return(self.xforms[self.index-1])



