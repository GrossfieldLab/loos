#!/usr/bin/env python3

class SymmMatrix:
    """
    1D storage of a symmetric matrix (eg a contact matrix)
    Useful as an intermediate to be sent to PCA, etc
    """

    import numpy as np

    def __init__(self, dimension, doStorage=False):
        self.dimension = int(dimension)
        self.doStorage = doStorage
        if self.doStorage:
            length = int(self.dimension * (self.dimension-1) / 2)
            self._array = np.zeros([length], np.float64)
        else:
            self._array = None

    def toFlat(self, i, j):
        """
        Given indices i,j to the symmetric matrix, generate the 1-D
        index
        """
        # canonicalize order
        if i > j:
            i, j = j, i

        index = 0
        for k in range(1, i+1):
            index += self.dimension - k
        index += j - i - 1
        return index

    def toSymm(self, index):
        """
        Given an index into the 1D representation, return the indices for the
        symmetric matrix (0-based)
        """
        i = 0
        j = 1
        val = 0
        next = 0
        while (val + next <= index):
            val += next
            i += 1
            next = self.dimension - i
        j = index - val + i
        return i-1, j

    def set(self, i, j, val):
        if not self.doStorage:
            return None
        index = self.toFlat(i, j)
        self._array[index] = val

    def incr(self, i, j, val):
        if not self.doStorage:
            return None
        index = self.toFlat(i, j)
        self._array[index] += val

    def get(self, i, j):
        if not self.doStorage:
            return None
        index = self.toFlat(i, j)
        return self._array[index]
