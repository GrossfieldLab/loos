# (c) 2015 Tod D. Romo, Grossfield Lab, URMC
import loos

class Trajectory(object):
    """
    This class turns a loos Trajectory into something more
    python-like.  Behind the scenes, it wraps a loos::AtomicGroup and
    a loos::Trajectory.

    """

    def __init__(self, fname, model, **kwargs):
        """
        Instantiate a Trajectory object.  Takes a filename and an
        AtomicGroup model.

        keyword arguments:
        skip -- # of frames to skip from the start
        stride -- # of frames to step through by
        iterator -- Python iterator used to pick frames
                    (overrides skip and stride)
        """
        skip = 0
        stride = 1
        iterator = None
        
        if 'skip' in kwargs:
            skip = kwargs['skip']
        if 'stride' in kwargs:
            stride = kwargs['stride']
        if 'iterator' in kwargs:
            iterator = kwargs['iterator']
        
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
        self.index = 0
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
        return(frame)


    def trajectory(self):
        return(self.traj)

    def readFrame(self, i):
        if (i < 0 or i >= len(self.framelist)):
            raise IndexError
        self.traj.readFrame(self.framelist[i])
        self.traj.updateGroupCoords(self.frame)     # ???

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








class VirtualTrajectory(object):

    def __init__(self, *trajs, **kwargs):
        self.skip = 0
        self.stride = 1
        self.nframes = 0
        self.iterator = None
        self.trajectories = list(trajs)

        self.index = 0
        self.framelist = []
        self.trajlist = [] 
        self.stale = 1

        if 'skip' in kwargs:
            self.skip = kwargs['skip']
        if 'stride' in kwargs:
            self.stride = kwargs['stride']
        if 'iterator' in kwargs:
            self.iterator = kwargs['iterator']
        
    def addTrajectory(self, *traj):
        self.trajectories.extend(traj)
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
        frames = []
        trajs = []
        for t in self.trajectories:
            for i in range(len(t)):
                frames.append(i)
                trajs.append(t)

        self.framelist = []
        self.trajlist = []
        n = len(frames)
        if (self.iterator is None):
            it = iter(range(self.skip, n, self.stride))
        else:
            it = iter(iterator)
        for i in it:
            self.framelist.append(frames[i])
            self.trajlist.append(trajs[i])

        self.index = 0
        self.stale = 0
            
    def __len__(self):
        if self.stale:
            self.initFrameList()
        return(len(self.framelist))

                
    def __getitem__(self, i):
        if self.stale:
            self.initFrameList()

        if isinstance(i, slice):
            return(self.getSlice(i))

        if (i < 0):
            i += len(self)
        if (i >= len(self)):
            raise IndexError

        return(self.trajlist[i][self.framelist[i]])


    def __iter__(self):
        if self.stale:
            self.initFrameList()
        self.index = 0
        return(self)

    def reset(self):
        self.index = 0

    def next(self):
        if self.stale:
            self.initFrameList()
        if (self.index >= len(self.framelist)):
            raise StopIteration
        frame = self.__getitem__(self.index)
        self.index += 1
        return(frame)

    def getSlice(self, s):
        indices = list(range(*s.indices(self.__len__())))
        ensemble = []
        for i in indices:
            frame = self.trajlist[i][self.framelist[i]].copy()
            ensemble.append(frame)
        return(ensemble)




class AlignedVirtualTrajectory(VirtualTrajectory):

    def __init__(self, *trajs, **kwargs):
        super(AlignedVirtualTrajectory, self).__init__(*trajs, **kwargs)
        self.aligned = False
        self.xformlist = []
        self.rmsd = 0
        self.iters = 0
        if ('alignwith' in kwargs):
            self.alignwith = kwargs['alignwith']
        else:
            self.alignwith = 'name == "CA"'

    def addTrajectory(self, *traj):
        self.trajectories.extend(traj)
        self.stale = True
        self.aligned = False
            
    def align(self):
        current_traj = None
        current_subset = None
        ensemble = []

        if self.stale:
            self.initFrameList()

        for i in range(len(self.framelist)):
            t = self.trajlist[i]
            if t != current_traj:
                current_traj = t
                current_subset = loos.selectAtoms(t.currentFrame(), self.alignwith)
            t.readFrame(self.framelist[i])
            ensemble.append(current_subset.copy())

        (self.xformlist, self.rmsd, self.iters) = loos.iterativeAlignEnsemble(ensemble)
        self.aligned = True


        
    def __getitem__(self, i):
        if not self.aligned:
            self.align()

        if (i < 0):
            i += len(self.framelist)
        if (i >= len(self.framelist)):
            raise IndexError

        frame = self.trajlist[i][self.framelist[i]]
        frame.applyTransform(self.xformlist[i])
        return(frame)
