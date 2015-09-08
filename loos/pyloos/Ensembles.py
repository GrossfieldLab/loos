import numpy

def averageStructure(traj):
    avg = numpy.zeros((len(traj.currentFrame()), 3))
    for frame in traj:
        coords = frame.getCoords()
        avg += coords
    avg /= len(traj)
    return(avg)
