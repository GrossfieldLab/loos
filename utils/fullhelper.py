#!/usr/bin/env python3
from sys import argv
import re


def periodic_insert(period, insert, target):
    l = len(target)
    if l < period:
        return target
    else:
        return target[:period] + insert + periodic_insert(period, insert, target[period:])

# IOJUNK
msg = ''
with open(argv[1], 'r') as f:
    msg = f.readlines()
wrapcount = int(argv[2])
targetfn = argv[3]

# a regex that matches nonwhitespace characters
nws = re.compile('\S')
for line in msg:
    if line.match(nws):
