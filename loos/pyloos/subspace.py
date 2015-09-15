"""
Functions related to subspace calculations...
"""


import numpy
import math


def subspaceOverlap(A, B, nmodes=0):
    (am, an) = A.shape
    (bm, bn) = B.shape
    if am != bm:
        raise RuntimeError('Eigenvectors must have the same number of rows')

    if nmodes == 0:
        nmodes = int(0.25 * min(an, bn))
    
    U = A[:,0:nmodes]
    V = B[:,0:nmodes]
    D = numpy.dot(numpy.transpose(U), V)
    D = numpy.multiply(D,D)
    return(numpy.sum(D) / nmodes)



def covarianceOverlap(ls, lU, rs, rU, nmodes = 0):

    (lm, ln) = lU.shape
    (rm, rn) = rU.shape

    if lm != rm :
        raise RuntimeError('Eigenvectors must have the same number of rows')
    if nmodes == 0:
        nmodes = min(ln, rn)
    
    X = numpy.absolute(numpy.dot(numpy.transpose(rU[:, 0:nmodes]), lU[:, 0:nmodes]))

    m = ls.shape
    if len(m) == 1:
        ls=numpy.reshape(ls, (m[0], 1))
        rs=numpy.reshape(rs, (m[0], 1))

    L = numpy.dot(rs[0:nmodes, :], numpy.transpose(ls[0:nmodes, :]))
    X2 = numpy.multiply(X, X)
    y = numpy.sum(numpy.multiply(numpy.sqrt(L), X2))

    e = numpy.sum(numpy.add(ls[0:nmodes], rs[0:nmodes]))
    
    num = e - 2.0 * y
    co = 1.0 - math.sqrt( math.fabs(num) / e)

    return(co)

    
