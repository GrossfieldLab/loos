"""
PyLOOS interface to the LOOS iterative alignment routines
"""


import loos



##    Iteratively align a loos.pyloos.Trajectory object (or a list of AtomicGroups)

def iterativeAlignment(ensemble, threshold=1e-8, maxiter=1000):
    """ 
    Iteratively align a loos.pyloos.Trajectory object (or a list of AtomicGroups).
    Returns the transformations needed to perform the alignment.  Note that it does
    not currently transform what's passed into it.
    
    Args

      ensemble: loos.pyloos.Trajectory, or list of AtomicGroups
      threshold (opt): change in average less than this ends alignment
      maxiter(opt): maximum number of iterations allowed


    Returns

      (list-of-xforms, final-rmsd, number-of-iterations)


    Examples

      model = loos.createSystem('foo.pdb')
      traj = loos.pyloos.Trajectory('foo.dcd', model, subset = 'name == "CA"')
      (xforms, rmsd, iters) = loos.pyloos.iterativeAlignment(traj)

      ensemble = [frame.copy() for frame in traj]
      (xforms, rmsd, iters) = loos.pyloos.iterativeAlignment(ensemble)
    """
    enlist = loos.DoubleVectorMatrix()
    for e in ensemble:
       enlist.push_back(e.coordsAsVector())
    result = loos.iterativeAlignmentPy(enlist, threshold, maxiter)
    return(loos.xformVectorToList(result.transforms), result.rmsd, result.iterations)



# Optional 'framelist' argument specifies indices of frames to use
def iterativeAlignTrajectory(model, traj, threshold=1e-8, maxiter=1000, **kwargs):
    """ 
    Interface to the standard LOOS iterative alignment routines.
    Returns the transformations needed to perform the alignment.  Note that it does
    not currently transform what's passed into it.

    Args

      model: AtomicGroup (subset of trajectory model) to use for aligning
      traj: trajectory to align
      threshold (opt): change in avg less than this ends alignment
      maxiter (opt): max number of iterations allowed
      framelist (opt): a python list (or loos.UIntVector) of frame indices
          from the trajectory to use


    Returns

      (list-of-xforms, final-rmsd, number-of-iterations)


    Examples

      model = loos.createSystem('foo.pdb')
      traj = loos.createTrajectory('foo.dcd', model)
      subset = loos.selectAtoms(model, 'backbone')

      (xforms, rmsd, iters) = loos.pyloos.iterativeAlignTrajectory(subset, traj)
    """
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
        result = loos.iterativeAlignmentPy(model, traj, flist, threshold, maxiter)

    else:
        result = loos.iterativeAlignmentPy(model, traj, threshold, maxiter)

    return(loos.xformVectorToList(result.transforms), result.rmsd, result.iterations)

