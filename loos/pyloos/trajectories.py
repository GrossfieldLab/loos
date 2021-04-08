"""
Python-based trajectory classes that wrap loos.Trajectory objects

"""
import loos
import copy

## Python-based wrapper for LOOS Trajectories
# This class turns a loos Trajectory into something more
# python-like.  Behind the scenes, it wraps a loos::AtomicGroup and
# a loos::Trajectory.  The behavior of the trajectory can be controlled
# through passed keywords,
#
# Keyword    | Description
# -----------|------------------------------------------------------------------------------
# skip=n     | Skip the first n-frames of the wrapped trajectory
# stride=n   | Step through the wrapped trajectory n-frames at a time
# iterator=i | Use the python iterator object i to select frames from the wrapped trajectory
# subset=s   | Use 's' to select a subset of the model to use for each frame
#
# Remember that all atoms are shared.  If you want to decouple the
# trajectory from other groups, pass it a copy of the model.
#
# examples:
# \code
#  model = loos.createSystem('foo.pdb')
#  traj = loos.pyloos.Trajectory('foo.dcd', model)
#  calphas = loos.selectAtoms(model, 'name == "CA"')
#  for frame in traj:
#      print calphas.centroid()
# \endcode
#
# The same thing but skipping the first 50 frames and taking every other frame
# \code
# traj = loos.pyloos.Trajectory('foo.dcd', model, skip=50, stride=2)
# \endcode
#
# Only use frames 19-39 (i.e. the 20th through 40th frames)
# \code
# traj = loos.pyloos.Trajectory('foo.dcd', model, iterator=range(19,40))
# \endcode
#
# An alternative way of only iterating over a subset...
# \code
# model = loos.createSystem('foo.pdb')
# traj = loos.pyloos.Trajectory('foo.dcd', model, subset='name == "CA"')
# \endcode
#
# Decouple the model stored in the trajectory,
# \code
# traj = loos.pyloos.Trajectory('foo.dcd', model.copy())
# \endcode
#

class Trajectory(object):
    """
    Python-based wrapper for LOOS Trajectories
    >>> model = loos.createSystem('foo.pdb')
    >>> traj = loos.pyloos.Trajectory('foo.dcd', model)

    keyword args:
        skip = # of frames to skip from start
      stride = # of frames to step through
    iterator = Python iterator used to pick frame (overrides skip and stride)
      subset = Selection used to pick subset for each frame

    See the Doxygen documentation for more details.
    """

    ## Instantiate a Trajectory object.
    def __init__(self, fname, model, **kwargs):

        self._skip = 0
        self._stride = 1
        self._iterator = None

        if 'skip' in kwargs:
            self._skip = kwargs['skip']
        if 'stride' in kwargs:
            self._stride = kwargs['stride']
        if 'iterator' in kwargs:
            self._iterator = kwargs['iterator']
        if 'subset' in kwargs:
            self._subset = loos.selectAtoms(model, kwargs['subset'])
        else:
            self._subset = model

        self._model = model
        self._fname = fname
        self._traj = loos.createTrajectory(fname, model)

        self._stale = 1
        self._initFrameList()


    def _initFrameList(self):
        self._framelist = []
        if self._iterator is None:
            it = range(self._skip, self._traj.nframes(), self._stride)
        else:
            it = iter(self._iterator)

        for i in it:
            self._framelist.append(i)

        self._index = 0
        self._stale = 0

    def stride(self, n):
        """
        Step through the trajectory by this number of frames
        """
        self._stride = n

    def skip(self, n):
        """
        Skip this number of frames at the start of the trajectory
        """
        self._skip = n

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
        if self._stale:
            self._initFrameList()
        self._index = 0
        return(self)

    def __len__(self):
        """
        Number of frames in the trajectory
        """
        if self._stale:
            self._initFrameList()
        return(len(self._framelist))


    def reset(self):
        """Reset the iterator"""
        self._index = 0

    def __next__(self):
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
        if self._stale:
            self._initFrameList()
        if (i < 0 or i >= len(self._framelist)):
            raise IndexError
        self._traj.readFrame(self._framelist[i])
        self._traj.updateGroupCoords(self._model)
        return(self._subset)

    def frame(self):
        """Return the current frame (subset)"""
        return(self._subset)

    def realIndex(self):
        """The 'real' frame in the trajectory for this index"""
        if self._stale:
            self._initFrameList()
        return(self._framelist[self._index-1])

    def index(self):
        """The state of the iterator"""
        return(self._index-1)


    def frameNumber(self, i):
        """
        Returns the real frame numbers corresponding to the passed indices.  Can accept
        either an integer or a list of integers.

        For example:
        >>> t = loos.pyloos.Trajectory('foo.dcd', model, skip=50)
        >>> t.frameNumber(0)
        50
        >>> t.frameNumber(range(0,2))
        [50, 51]
        """

        if self._stale:
            self._initFrameList()
        if type(i) is int:
            if (i < 0):
                i += len(self._framelist)
            return(self._framelist[i])

        indices = [x if x >=0 else len(self._framelist)+x for x in i]
        framenos = [self._framelist[x] for x in indices]
        return(framenos)


    def _getSlice(self, s):
        if self._stale:
            self._initFrameList()
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
        if self._stale:
            self._initFrameList()
        if isinstance(i, slice):
            return(self._getSlice(i))

        if (i < 0):
            i += len(self._framelist)
        if (i >= len(self._framelist) or i < 0):
            raise IndexError
        self._traj.readFrame(self._framelist[i])
        self._traj.updateGroupCoords(self._model)
        return(self._subset)



## Virtual trajectory composed of multiple Trajectory objects
# This class can combine multiple loos.pyloos.Trajectory objects
# into one big "virtual" trajectory.  Any skips or strides set in
# the contained trajectories will be honored.  In addition, a skip
# and a stride for the whole meta-trajectory are available.  These
# can be set via keyword arguments when creating a VirtualTrajectory,
#
# Keyword    | Description
# -----------|------------------------------------------------------------------------------
# skip=n     | Skip the first n-frames of the virtual trajectory
# stride=n   | Step through the virtual trajectory n frames at a time
# iterator=i | Use the python iterator object i to select frames from the virtual trajectory
#
# There is no requirement that the subsets used for all trajectories
# must be the same.  Ideally, the frame (subset) that is returned
# should be compatible (e.g. same atoms in the same order), but the
# models used for each trajectory (and the corresponding subset
# selection) can vary.  Why would you want to do this?  Imagine
# combining three different GPCRs where the subsets are the common
# trans-membrane C-alphas.  This makes processing all of the
# ensembles together easier.
#
# <h2>WARNING</h2>
# Since each contained trajectory can have a different set of shared
# atoms it updates, care must be taken when pre-selecting atoms.
#
#
# Examples:
# \code
# model = loos.createSystem('foo.pdb')
# Takes a filename (\a fname) and an AtomicGroup model (\a model).
# Takes a filename (\a fname) and an AtomicGroup model (\a model).
# Takes a filename (\a fname) and an AtomicGroup model (\a model).
# traj1 = loos.pyloos.Trajectory('foo-1.dcd', model)
# traj2 = loos.pyloos.Trajectory('foo-2.dcd', model)
# vtraj = loos.pyloos.VirtualTrajectory(traj1, traj2)
#
# for frame in vtraj:
#   print frame.centroid()
# \endcode
#
# Adding another trajectory in with its own stride and skip
# \code
# traj3 = loos.pyloos.Trajectory('foo-3.dcd', skip=50, stride=2)
# vtraj.append(traj3)
# \endcode
#
# Same as above but stride through the combined trajectory
# \code
# vtraj10 = loos.pyloos.VirtualTrajectory(traj1, traj2, stride=10)
# \endcode

class VirtualTrajectory(object):
    """
    Combines multiple loos.pyloos.Trajectory objects into one big virtual trajectory
    >>> vtraj = loos.pyloos.VirtualTrajectory(traj1, traj2, traj3)
        Keyword arguments:
            skip = # of frames to skip at start of composite traj
          stride = # of frames to step through in the composite traj
        iterator = Python iterator used to pick frames from the composite traj

    See the Doxygen documentation for more details.
    """


    def __init__(self, *trajs, **kwargs):
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

        # make sure the trajectories are trajectories
        for t in self._trajectories:
            if not isinstance(t, loos.pyloos.trajectories.Trajectory):
                raise TypeError("Inputs to VirtualTrajectory must be pyloos Trajectory objects")


    def append(self, *traj):
        """
        Add a trajectory to the end of the virtual trajectory.  Resets
        the iterator state
        """

        # make sure we're appending only Trajectory objects
        if isinstance(traj, tuple):
            for t in traj:
                if not isinstance(t, loos.pyloos.trajectories.Trajectory):
                    raise TypeError("Inputs to VirtualTrajectory must be pyloos Trajectory objects")
        else:
            if not isinstance(traj, loos.pyloos.trajectories.Trajectory):
                    raise TypeError("Inputs to VirtualTrajectory must be pyloos Trajectory objects")

        self._trajectories.extend(traj)
        self._stale = 1

    def stride(self, n):
        """
        Set the stride of the combined trajectory
        """
        self._stride = n
        self._stale = 1

    def skip(self, n):
        """
        Set the skip of the combined trajectory
        """
        self._skip = n
        self._stale = 1

    def allStride(self, n):
        """
        Sets the stride of all contained trajectories
        """
        self._stale = 1
        for t in self._trajectories:
            t.stride(n)

    def allSkip(self, n):
        """
        Sets the skip of all contained trajectories
        """
        self._stale = 1
        for t in self._trajectories:
            t.skip(n)

    def setSubset(self, selection):
        """
        Set the subset selection for all managed trajectories
        """
        for t in self._trajectories:
            t.setSubset(selection)

    def frame(self):
        """
        Return the current frame/model.  If the iterator is past the
        end of the trajectory list, return the last valid frame.
        """
        if self._stale:
            self._initFrameList()

        if self._index >= len(self._framelist):
            i = len(self._framelist) - 1
        else:
            i = self._index

        return(self._trajectories[self._trajlist[i]].frame())

    def index(self):
        """
        Return index into composite trajectory for current frame
        """
        return(self._index-1)



    ## Returns information about the ith frame in the VirtualTrajectory
    # The tuple returned has the following format:
    # \code
    # (frame-index, traj-index, trajectory, real-frame-within-trajectory)
    # \endcode
    #
    # Consider the following,
    # \code
    # t1 = loos.pyloos.Trajectory('foo.dcd', model)   # 50 frames
    # t2 = loos.pyloos.Trajectory('bar.dcd', model)   # 25 frames
    # vt = loos.pyloos.VirtualTrajectory(t1, t2)
    # \endcode
    #
    # * \c frame-index is the index into the corresponding trajectory object.  For example,
    # frameLocation(50) would have a frame-index of 0 because vt[50] would return the first
    # frame from t2.
    #
    # * \c traj-index is the index into the list of managed trajectories for the frame.
    # In the above example, the traj-index will be 1.
    # * \c trajectory is the actual loos.pyloos.Trajectory object that contains the frame.
    # * \c real-frame-within-trajectory is the same as calling trajectory.frameNumber(frame-index).
    #
    # Instead of the \c t1 above, imagine it was setup this way,
    # \code
    # t1 = loos.pyloos.Trajectory('foo.dcd', model, skip=25)
    # \endcode
    # Now, <tt>vt.frameLocation(0)</tt> will return <tt>(0, 0, t1, 25)</tt>,
    # and <tt>vt.frameLocation(25)</tt> will return <tt>(25, 1, t2, 0)</tt>
    #
    # Python documentation:

    def frameLocation(self, i):
        """
        Return info about where a frame comes from.
        >>> (frame-index, traj-index, trajectory, real-frame-within-trajectory) = vtraj.frameLocation(i)
        """
        if (self._stale):
            self._initFrameList()

        if (i < 0):
            i += len(self._framelist)

        t = self._trajectories[self._trajlist[i]]
        return( self._framelist[i], self._trajlist[i], t, t.frameNumber(self._framelist[i]))

    def frameBoundaries(self):
        """
        Return a list containing the index of the first frame associated with
        each traj
        >>> b = vt.frameBoundaries()
        len(b) will be the number of trajectories in vt.
        -> can slice the data from the nth traj from an array with b[n]:b[n+1]
        """
        from numpy import searchsorted
        if (self._stale):
            self._initFrameList()
        boundaries = [0]
        for i in range(1, len(self._trajectories)):
            loc = searchsorted(_trajlist, i)
            boundaries.append(loc)
        return boundaries

    def _initFrameList(self):
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
            self._initFrameList()
        return(len(self._framelist))


    def __getitem__(self, i):
        """
        Return the ith frame in the composite trajectory.  Supports
        Python slicing.  Negative indices are relative to the end of
        the composite trajectory.
        """
        if self._stale:
            self._initFrameList()

        if isinstance(i, slice):
            return(self._getSlice(i))

        if (i < 0):
            i += len(self)
        if (i >= len(self)):
            raise IndexError

        return(self._trajectories[self._trajlist[i]][self._framelist[i]])


    def __iter__(self):
        if self._stale:
            self._initFrameList()
        self._index = 0
        return(self)

    def reset(self):
        self._index = 0

    def __next__(self):
        if self._stale:
            self._initFrameList()
        if (self._index >= len(self._framelist)):
            raise StopIteration
        frame = self.__getitem__(self._index)
        self._index += 1
        return(frame)

    def _getSlice(self, s):
        indices = list(range(*s.indices(self.__len__())))
        ensemble = []
        for i in indices:
            frame = self._trajectories[self._trajlist[i]][self._framelist[i]].copy()
            ensemble.append(frame)
        return(ensemble)



## A virtual trajectory that supports iterative alignment.
# Only the
# transformation needed to align each frame is stored.  When a frame
# is accessed, it is automatically transformed into the aligned
# orientation.  All keywords from VirtualTrajectory are supported,
# along with the following new ones,
#
# Keyword     | Description
# ------------|------------------------------------------------------------------------------
# alignwith=s | Use 's' to select what part of the model is used for aligning.
# reference=g | Use the AtomicGroup g as a reference structure.  All frames will be aligned to it.
#
# There are two ways that a trajectory can be aligned.  The first
# uses in iterative alignment method (the same used in LOOS).  This
# is the default method.  In order to do the alignment, the
# alignwith subset must be read into memory and temporarily stored.
# This can potentially use a lot of memory and create delays in
# execution.  Once the alignment is complete, however, those cached
# frames are released and subsequent frame accesses will be quick.
#
# The second method is to align each frame to a reference
# structure.  This method is selected when a reference structure is
# passed to the constructor (with the 'reference' keyword), or when
# setReference() is called.  Note that you can pass None to
# setReference() which will return the AlignedVirtualTrajectory to
# the iterative method.  Also note that the reference structure is
# copied into the AVT object as a deep copy (i.e. it does not share
# any atoms).
#
# See VirtualTrajectory for some basic examples in addition to
# below...
#
# Align using only C-alphas (the default)
# \code
# vtraj = loos.pyloos.AlignedVirtualTrajectory(traj1, traj2)
# \endcode
#
# Align using only backbone atoms
# \code
# vtraj = loos.pyloos.AlignedVirtualTrajectory(traj1, traj2, alignwith='name =~ "^(C|N|O|CA)$"')
# \endcode
#
# Add another trajectory
# \code
# vtraj.append(traj3)
# \endcode
#
# Align using only C-alphas and a reference structure
# \code
# refmodel = loos.createSystem('foo-ref.pdb')
# refsubset = loos.selectAtoms(refmodel, 'name == "CA"')
# vtraj = loos.pyloos.AlignedVirtualTrajectory(traj1, traj2, reference = refsubset)
# \endcode

class AlignedVirtualTrajectory(VirtualTrajectory):
    """
    Combine loos.pyloos.Trajectory objects and align them

    >>> aligned = loos.pyloos.AlignedVirtualTrajectory(traj1, traj2, alignwith = 'backbone')
    Supports the same keywords as VirtualTrajectory.
    New keywords:
      alignwith = Selection used for alignment (default is all C-alphas)
      reference = AtomicGroup that all frames are aligned to (disables iterative alignment)

    See the Doxygen documentation for more details.
    """

    def __init__(self, *trajs, **kwargs):
        super(AlignedVirtualTrajectory, self).__init__(*trajs, **kwargs)
        self._aligned = False
        self._xformlist = []
        self._rmsd = -1
        self._iters = -1
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
        # make sure we're appending only Trajectory objects
        if isinstance(traj, tuple):
            for t in traj:
                if not isinstance(t, loos.pyloos.trajectories.Trajectory):
                    raise TypeError("Inputs to VirtualTrajectory must be pyloos Trajectory objects")
        else:
            if not isinstance(traj, loos.pyloos.trajectories.Trajectory):
                    raise TypeError("Inputs to VirtualTrajectory must be pyloos Trajectory objects")

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
        self._align()
        self._index = 0
        return(self)


    def setReference(self, reference):
        self._reference = copy.deepcopy(reference)
        self._aligned = False

    def _align(self):
        """
        Align the frames (called implicitly on iterator or array access)
        """
        current_traj = None
        current_subset = None
        ensemble = []

        if self._stale:
            self._initFrameList()

        if self._reference:       # Align to a reference structure
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

        else:                      # Iterative alignment

            ensemble = loos.DoubleVectorMatrix()

            for i in range(len(self._framelist)):
                t = self._trajectories[self._trajlist[i]]
                if t != current_traj:
                    current_traj = t
                    current_subset = loos.selectAtoms(t.model(), self._alignwith)
                t.readFrame(self._framelist[i])
                ensemble.push_back(current_subset.coordsAsVector())

            result = loos.iterativeAlignmentPy(ensemble)
            (self._xformlist, self._rmsd, self._iters) = (loos.xformVectorToList(result.transforms), result.rmsd, result.iterations)

        self._aligned = True


    def rmsd(self):
        return(self._rmsd)

    def iters(self):
        return(self._iters)


    def _getSlice(self, s):
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
            self._align()

        if isinstance(i, slice):
            return(self._getSlice(i))

        if (i < 0):
            i += len(self._framelist)
        if (i >= len(self._framelist)):
            raise IndexError

        frame = self._trajectories[self._trajlist[i]][self._framelist[i]]
        frame.applyTransform(self._xformlist[i])
        return(frame)
