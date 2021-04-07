#!/usr/bin/env python3

import argparse

# Set up the parser
parser = argparse.ArgumentParser()


# Set up some default arguments
def modelSelectionOptions(parser):
    parser.add_argument('--model',
                        help="Model file describing system contents")
    parser.add_argument('--sel',
                        help='Use this selection for computation',
                        default='all')
    return parser


def trajOptions(parser):
    parser.add_argument('--traj',
                        help='Filename of trajectory or trajectories',
                        nargs='+')
    parser.add_argument('--skip',
                        help='Skip frames from the start of the trajectory',
                        type=int,
                        default=0)
    parser.add_argument('--stride',
                        help='Step through the trajectory by this many frames',
                        type=int,
                        default=1)
    return parser
