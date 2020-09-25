#!/usr/bin/env python3
import sys
import json
import numpy
from math import sqrt


def read_frame_ranges(filename):
    filenames = []
    ends = []
    with open(filename) as f:
        # skip the first 2 lines, stop when a line doesn't start with "#"
        for line in f.readlines()[2:]:
            if line.startswith("#"):
                try:
                    hash, index, start, end, filename = line.split()
                    filenames.append(filename)
                    ends.append(int(end))
                except ValueError:
                    break
            else:
                break
    return (filenames, numpy.array(ends))


if __name__ == "__main__":
    json_filename = sys.argv[1]
    rmsds_filename = sys.argv[2]
    split = int(sys.argv[3])

    filenames, ends = read_frame_ranges(rmsds_filename)

    with open(json_filename) as f:
        cluster_data = json.load(f)

    clusters = cluster_data['clusters']

    data = numpy.zeros((len(clusters), len(filenames)))

    for i in range(len(clusters)):
        arr = numpy.array(clusters[i])
        # TODO: check this
        indices = numpy.searchsorted(ends, arr, side='right')
        indices -= 1
        for j in indices:
            data[i, j] += 1

    # normalize for each trajectory
    data /= numpy.add.reduce(data)

    print("# Clust Ave Dev StdErr")
    first_mean = numpy.mean(data[:,:split], axis=1)
    first_dev = numpy.std(data[:,:split], axis=1)
    first_err = first_dev / sqrt(split) # assumes same number of trajs in both
    sec_mean = numpy.mean(data[:,split:], axis=1)
    sec_dev = numpy.std(data[:, split:], axis=1)
    sec_err = sec_dev / sqrt(split)
    for i in range(len(first_mean)):
        print(i, first_mean[i], first_dev[i], first_err[i],
              sec_mean[i], sec_dev[i], sec_err[i]
             )
