"""
PyLOOS functions applicable to ensembles (trajectories or trajectory-like collections)

"""
import loos
import numpy
import numpy.linalg


## Python version of averageStructure (using loos.pyloos.Trajectory-like objects)
# The subset defined in the trajectory controls what is averaged over.  The trajectory
# may actually be a VirtualTrajectory or an AlignedVirtualTrajectory.

def averageStructure(traj):
    """Returns the average structure for a trajectory.
    >>> avg = loos.pyloos.averageStructure(traj)
    """
    avg = numpy.zeros((len(traj.frame()), 3))
    for frame in traj:
        coords = frame.getCoords()
        avg += coords
    avg /= len(traj)

    structure = traj.frame().copy()
    for i in range(len(structure)):
        structure[i].coords(loos.GCoord(avg[i][0], avg[i][1], avg[i][2]))
        
    return(structure)

## Returns the coordinates for an entire trajectory as an MxN numpy matrix
# where M is 3*natoms and N is the length of the trajectory.  The subset
# in the trajectory controls what is extracted.  The trajectory
# may actually be a VirtualTrajectory or an AlignedVirtualTrajectory.

def extractCoords(traj):
    """
    Extracts coords from a trajectory as a NumPy matrix
    >>> A = loos.pyloos.extractCoords(traj)
    """
    m = len(traj.frame()) * 3
    n = len(traj)

    A = numpy.zeros((m, n))
    for i in range(n):
        coords = traj[i].getCoords()
        A[:,i] = numpy.reshape(coords, (1,m))
    return(A)



## Returns a tuple containing the SVD result for a trajectory, along with the average structure.
# The tuple is,
# <tt>(left-singular-vectors, singular-values, right-singular-vectors, average)</tt>
#
# The subset set in the trajectory controls what is used for the SVD.  The
# trajectory may be a VirtualTrajectory or an AlignedVirtualTrajectory.
# The following example computes an SVD using backbone atoms for a trajectory
# that has been iteratively aligned using C-alphas
# \code
#traj = loos.pyloos.Trajectory('foo.dcd', model, subset = 'backbone')
#aligned = loos.pyloos.AlignedVirtualTrajectory(traj, alignwith = 'name == "CA"')
#(L,S,V,avg) = loos.pyloos.svd(aligned)
# \endcode


def svd(traj):
    """
    Returns a tuple containing SVD results along with the average structure for a trajectory
    >>> (L,S,V,avg) = loos.pyloos.svd(traj)
    """
    A = extractCoords(traj)
    avg = numpy.average(A, 1)
    avg = avg[:,numpy.newaxis]
    A -= avg
    (U,s,V) = numpy.linalg.svd(A)

    structure = traj.frame().copy()
    for i in range(len(structure)):
        structure[i].coords(loos.GCoord(avg[i*3], avg[i*3+1], avg[i*3+2]))

    return(U,s,V,structure)
