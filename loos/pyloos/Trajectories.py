# (c) 2015 Tod D. Romo, Grossfield Lab, URMC
import loos
import copy

class Trajectory(object):
    """
    This class turns a loos Trajectory into something more
    python-like.  Behind the scenes, it wraps a loos::AtomicGroup and
    a loos::Trajectory.

    Remember that all atoms are shared.  If you want to decouple the
    trajectory from other groups, pass it a copy of the model.


    Examples:
        model = loos.createSystem('foo.pdb')
        traj = loos.pyloos.Trajectory('foo.dcd', model)
        calphas = loos.selectAtoms(model, 'name == "CA"')
        for frame in traj:
            print calphas.centroid()



        # The same thing but skipping the first 50 frames
        # and taking every other frame
        traj = loos.pyloos.Trajectory('foo.dcd', model, skip=50, stride=2)

        # Only use frames 19-39 (i.e. the 20th through 40th frames)
        traj = loos.pyloos.Trajectory('foo.dcd', model, iterator=range(19,40))

        # An alternative way of only iterating over a subset...
        model = loos.createSystem('foo.pdb')
        traj = loos.pyloos.Trajectory('foo.dcd', model, subset='name == "CA"')


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
          subset -- Use this selection as the subset to return for frames...
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
        if 'subset' in kwargs:
            self._subset = loos.selectAtoms(model, kwargs['subset'])
        else:
            self._subset = model
            
        self._model = model
        self._fname = fname
        self._traj = loos.createTrajectory(fname, model)

        self._framelist = []
        if iterator is None:
            it = range(skip, self._traj.nframes(), stride)
        else:
            it = iter(iterator)
        
        for i in it:
            self._framelist.append(i)

        self._index = 0

    def fileName(self):
        """
        Return the filename that this Trajectory represents
        """
        return(self._fname)
        
    def setSubset(self, selection):
        """
        Set the subset used when iterating over a trajectory.
        The selection is a LOOS selection string.
        """
        self._subset = loos.selectAtoms(self._model, selection)

        
    def __iter__(self):
        self._index = 0
        return(self)

    def __len__(self):
        """
        Number of frames in the trajectory
        """
        return(len(self._framelist))


    def reset(self):
        """Reset the iterator"""
        self._index = 0

    def next(self):
        if (self._index >= len(self._framelist)):
            raise StopIteration
        frame = self.__getitem__(self._index)

        self._index += 1
        return(frame)


    def trajectory(self):
        """Access the wrapped loos.Trajectory"""
        return(self._traj)

    def model(self):
        """Return the current model"""
        return(self._model)
    
    def readFrame(self, i):
        """Read a frame and update the model"""
        if (i < 0 or i >= len(self._framelist)):
            raise IndexError
        self._traj.readFrame(self._framelist[i])
        self._traj.updateGroupCoords(self._model)
        return(self._subset)

    def currentFrame(self):
        """Return the current frame (subset)"""
        return(self._subset)

    def currentRealIndex(self):
        """The 'real' frame in the trajectory for this index"""
        return(self._framelist[self._index-1])

    def currentIndex(self):
        """The state of the iterator"""
        return(self._index-1)


    def frameNumber(self, i):
        """
        Returns the real frame numbers corresponding to the passed indices.  Can accept
        either an integer or a list of integers.

        For example,
            t = loos.pyloos.Trajectory('foo.dcd', model, skip=50)

            t.frameNumber(0)           == 50
            t.frameNumber(range(0,2))  == [50,51]
        """
        if type(i) is int:
            if (i < 0):
                i += len(self._framelist)
            return(self._framelist[i])
        
        indices = [x if x >=0 else len(self._framelist)+x for x in i]
        framenos = [self._framelist[x] for x in indices]
        return(framenos)


    def getSlice(self, s):
        indices = list(range(*s.indices(self.__len__())))
        ensemble = []
        for i in indices:
            self._traj.readFrame(self._framelist[i])
            self._traj.updateGroupCoords(self._model)
            dup = self._subset.copy()
            ensemble.append(dup)
        return(ensemble)


    def __getitem__(self, i):
        """Handle array indexing and slicing.  Negative indices are
        relative to the end of the trajectory"""
        if isinstance(i, slice):
            return(self.getSlice(i))

        if (i < 0):
            i += len(self._framelist)
        if (i >= len(self._framelist) or i < 0):
            raise IndexError
        self._traj.readFrame(self._framelist[i])
        self._traj.updateGroupCoords(self._model)
        return(self._subset)



class VirtualTrajectory(object):
    """
    This class can combine multiple loos.pyloos.Trajectory objects
    into one big "virtual" trajectory.  Any skips or strides set in
    the contained trajectories will be honored.  In addition, a skip
    and a stride for the whole meta-trajectory are available.

    There is no requirement that the subsets used for all trajectories
    must be the same.  Ideally, the frame (subset) that is returned
    should be compatible (e.g. same atoms in the same order), but the
    models used for each trajectory (and the corresponding subset
    selection) can vary.  Why would you want to do this?  Imagine
    combining three different GPCRs where the subsets are the common
    trans-membrane C-alphas.  This makes processing all of the
    ensembles together easier.

    WARNING...
    Since each contained trajectory can have a different set of shared
    atoms it updates, care must be taken when pre-selecting atoms.


    Examples:
      model = loos.createSystem('foo.pdb')
      traj1 = loos.pyloos.Trajectory('foo-1.dcd', model)
      traj2 = loos.pyloos.Trajectory('foo-2.dcd', model)
      vtraj = loos.pyloos.VirtualTrajectory(traj1, traj2)

      for frame in vtraj:
          print frame.centroid()


      # Adding another trajectory in with its own stride and skip
      traj3 = loos.pyloos.Trajectory('foo-3.dcd', skip=50, stride=2)
      vtraj.append(traj3)

      # Same as above but stride through the combined trajectory
      vtraj10 = loos.pyloos.VirtualTrajectory(traj1, traj2, stride=10)

    

    """
    def __init__(self, *trajs, **kwargs):
        """
        Instantiate a VirtualTrajectory object.

        Keyword arguments:
            skip -- # of frames to skip at start of composite traj
          stride -- # of frames to step through in the composite traj
        iterator -- Python iterator used to pick frames from the composite traj

        """

        self._skip = 0
        self._stride = 1
        self._nframes = 0
        self._iterator = None
        self._trajectories = list(trajs)

        self._index = 0
        self._framelist = []
        self._trajlist = [] 
        self._stale = 1

        if 'skip' in kwargs:
            self._skip = kwargs['skip']
        if 'stride' in kwargs:
            self._stride = kwargs['stride']
        if 'iterator' in kwargs:
            self._iterator = kwargs['iterator']
        if 'subset' in kwargs:
            self.setSubset(kwargs['subset'])

        
    def append(self, *traj):
        """
        Add a trajectory to the end of the virtual trajectory.  Resets
        the iterator state
        """
        self._trajectories.extend(traj)
        self._stale = 1


    def setSubset(self, selection):
        """
        Set the subset selection for all managed trajectories
        """
        for t in self._trajectories:
            t.setSubset(selection)

    def currentFrame(self):
        """
        Return the current frame/model.  If the iterator is past the
        end of the trajectory list, return the last valid frame.
        """
        if self._stale:
            self.initFrameList()

        if self._index >= len(self._framelist):
            i = len(self._framelist) - 1
        else:
            i = self._index
            
        return(self._trajectories[self._trajlist[i]].currentFrame())

    def currentIndex(self):
        """
        Return index into composite trajectory for current frame
        """
        return(self._index-1)
    

    def currentTrajectoryIndex(self):
        """
        Returns the index into the list of trajectories that the
        current frame is from
        """
        return(self._trajlist[self._index-1])

    def currentTrajectory(self):
        """
        Returns the loos.pyloos.Trajectory object that the current
        frame is from
        """
        return(self._trajectories[self._trajlist[self._index-1]])

    def frameLocation(self, i):
        """
        Returns a tuple containing information about where a frame in the virtual trajectory 
        comes from.

        (frame-index, traj-index, trajectory, real-frame-within-trajectory)

        Consider the following:
          t1 = loos.pyloos.Trajectory('foo.dcd', model)   # 50 frames
          t2 = loos.pyloos.Trajectory('bar.dcd', model)   # 25 frames
          vt = loos.pyloos.VirtualTrajectory(t1, t2)

        The frame-index is the index into the corresponding trajectory object.  For example,
        frameLocation(50) would have a frame-index of 0 because vt[50] would return the first
        frame from t2.

        The traj-index is the index into the list of managed trajectories for the frame.
        In the above example, the traj-index will be 1.

        The trajectory is the actual loos.pyloos.Trajectory object that contains the frame.

        The real-frame-within-trajectory is the same as calling trajectory.frameNumber(frame-index).
        Instead of the t1 above, imagine it was setup this way,
          t1 = loos.pyloos.Trajectory('foo.dcd', model, skip=25)
        Now, vt.frameLocation(0) will return (0, 0, t1, 25),
        and vt.frameLocation(25) will return (25, 1, t2, 0)
        
        """
        if (self._stale):
            self.initFrameList()
            
        if (i < 0):
            i += len(self._framelist)

        t = self._trajectories[self._trajlist[i]]
        return( self._framelist[i], self._trajlist[i], t, t.frameNumber(self._framelist[i]))
    
    def initFrameList(self):
        frames = []
        trajs = []
        for j in range(len(self._trajectories)):
            t = self._trajectories[j]
            for i in range(len(t)):
                frames.append(i)
                trajs.append(j)

        self._framelist = []
        self._trajlist = []
        n = len(frames)
        if (self._iterator is None):
            it = iter(range(self._skip, n, self._stride))
        else:
            it = iter(iterator)
        for i in it:
            self._framelist.append(frames[i])
            self._trajlist.append(trajs[i])

        self._index = 0
        self._stale = 0
            
    def __len__(self):
        """
        Total number of frames
        """
        if self._stale:
            self.initFrameList()
        return(len(self._framelist))

                
    def __getitem__(self, i):
        """
        Return the ith frame in the composite trajectory.  Supports
        Python slicing.  Negative indices are relative to the end of
        the composite trajectory.
        """
        if self._stale:
            self.initFrameList()

        if isinstance(i, slice):
            return(self.getSlice(i))

        if (i < 0):
            i += len(self)
        if (i >= len(self)):
            raise IndexError

        return(self._trajectories[self._trajlist[i]][self._framelist[i]])


    def __iter__(self):
        if self._stale:
            self.initFrameList()
        self._index = 0
        return(self)

    def reset(self):
        self._index = 0

    def next(self):
        if self._stale:
            self.initFrameList()
        if (self._index >= len(self._framelist)):
            raise StopIteration
        frame = self.__getitem__(self._index)
        self._index += 1
        return(frame)

    def getSlice(self, s):
        indices = list(range(*s.indices(self.__len__())))
        ensemble = []
        for i in indices:
            frame = self._trajectories[self._trajlist[i]][self._framelist[i]].copy()
            ensemble.append(frame)
        return(ensemble)




class AlignedVirtualTrajectory(VirtualTrajectory):
    """
    A virtual trajectory that supports iterative alignment.  Only the
    transformation needed to align each frame is stored.  When a frame
    is accessed, it is automatically transformed into the aligned
    orientation.  The selection used for aligning can be set with the
    'alignwith' argument to the constructor (or the alignWith()
    method).


    There are two ways that a trajectory can be aligned.  The first
    uses in iterative alignment method (the same used in LOOS).  This
    is the default method.  In order to do the alignment, the
    alignwith subset must be read into memory and temporarily stored.
    This can potentially use a lot of memory and create delays in
    execution.  Once the alignment is complete, however, those cached
    frames are released and subsequent frame accesses will be quick.

    The second method is to align each frame to a reference
    structure.  This method is selected when a reference structure is
    passed to the constructor (with the 'reference' keyword), or when
    setReference() is called.  Note that you can pass None to
    setReference() which will return the AlignedVirtualTrajectory to
    the iterative method.  Also note that the reference structure is
    copied into the AVT object as a deep copy (i.e. it does not share
    any atoms).

    See VirtualTrajectory for some basic examples in addition to
    below:

    # Align using only C-alphas (the default)
    vtraj = loos.pyloos.AlignedVirtualTrajectory(traj1, traj2)

    # Align using only backbone atoms
    vtraj = loos.pyloos.AlignedVirtualTrajectory(traj1, traj2, alignwith='name =~ "^(C|N|O|CA)$"')

    # Add another trajectory
    vtraj.append(traj3)

    # Align using only C-alphas and a reference structure
    refmodel = loos.createSystem('foo-ref.pdb')
    refsubset = loos.selectAtoms(refmodel, 'name == "CA"')
    vtraj = loos.pyloos.AlignedVirtualTrajectory(traj1, traj2, reference = refsubset)

    """

    def __init__(self, *trajs, **kwargs):
        super(AlignedVirtualTrajectory, self).__init__(*trajs, **kwargs)
        self._aligned = False
        self._xformlist = []
        self._rmsd = 0
        self._iters = 0
        if 'alignwith' in kwargs:
            self._alignwith = kwargs['alignwith']
        else:
            self._alignwith = 'name == "CA"'

        if 'reference' in kwargs:
            self._reference = copy.deepcopy(kwargs['reference'])
        else:
            self._reference = None


    def append(self, *traj):
        """
        Add another trajectory at the end.  Requires re-aligning
        """
        self._trajectories.extend(traj)
        self._stale = True
        self._aligned = False

    def alignWith(self, selection):
        """
        Change the selection used to align with.  Requires re-aligning
        """
        self._alignwith = selection
        self._aligned = False


    def __iter__(self):
        self.align()
        self._index = 0
        return(self)


    def setReference(self, reference):
        self._reference = copy.deepcopy(reference)
        self._aligned = False
    
    def align(self):
        """
        Align the frames (called implicitly on iterator or array access)
        """
        current_traj = None
        current_subset = None
        ensemble = []

        if self._stale:
            self.initFrameList()

        if self._reference:
            self._xformlist = []
            for i in range(len(self._framelist)):
                t = self._trajectories[self._trajlist[i]]
                if t != current_traj:
                    current_traj = t
                    current_subset = loos.selectAtoms(t.model(), self._alignwith)
                t.readFrame(self._framelist[i])
                m = current_subset.superposition(self._reference)
                x = loos.XForm()
                x.load(m)
                self._xformlist.append(x)

            self._rmsd = 0.0
            self._iters = 0

        else:
            
            for i in range(len(self._framelist)):
                t = self._trajectories[self._trajlist[i]]
                if t != current_traj:
                    current_traj = t
                    current_subset = loos.selectAtoms(t.model(), self._alignwith)
                t.readFrame(self._framelist[i])
                ensemble.append(current_subset.copy())

            (self._xformlist, self._rmsd, self._iters) = loos.iterativeAlignEnsemble(ensemble)

        self._aligned = True

        
    def getSlice(self, s):
        indices = list(range(*s.indices(self.__len__())))
        ensemble = []
        for i in indices:
            frame = self._trajectories[self._trajlist[i]][self._framelist[i]].copy()
            frame.applyTransform(self._xformlist[i])
            ensemble.append(frame)
        return(ensemble)

        
    def __getitem__(self, i):
        """
        Returns the ith frame aligned.  Supports Python slices.  Negative indices are relative
        to the end of the composite trajectory.
        """
        if not self._aligned:
            self.align()

        if isinstance(i, slice):
            return(self.getSlice(i))
        
        if (i < 0):
            i += len(self._framelist)
        if (i >= len(self._framelist)):
            raise IndexError

        frame = self._trajectories[self._trajlist[i]][self._framelist[i]]
        frame.applyTransform(self._xformlist[i])
        return(frame)
