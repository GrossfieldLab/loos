import numpy
import numpy.linalg

def averageStructure(traj):
    avg = numpy.zeros((len(traj.currentFrame()), 3))
    for frame in traj:
        coords = frame.getCoords()
        avg += coords
    avg /= len(traj)
    return(avg)


def extractCoords(traj):
    m = len(traj.currentFrame()) * 3
    n = len(traj)

    A = numpy.zeros((m, n))
    for i in range(n):
        coords = traj[i].getCoords()
        A[:,i] = numpy.reshape(coords, (1,m))
    return(A)



def svd(traj):
    A = extractCoords(traj)
    avg = numpy.average(A, 1)
    avg = avg[:,numpy.newaxis]
    A -= avg
    (U,s,V) = numpy.linalg.svd(A)
    return(U,s,V,avg)
