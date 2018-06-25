#!/usr/bin/env python3
from sys import argv

# int period string head string tail string target
# prints head + target % period + tail, calls itself on remainder of target
# returns void
def periodic_insert(period, head, tail, target):
    l = len(target)
    if l < period:
        print(head + target + tail)
    else:
        print(head + target[:period] + tail) 
        periodic_insert(period, head, tail, target[period:])

# IOJUNK
msg = []
with open(argv[1], 'r') as f:
    msg = f.readlines()
wrapcount = int(argv[2])
signal = argv[4]
head = '\"'
tail = '\\n\"'


with open(argv[3], 'r') as loos_source_file:
    for source_line in loos_source_file:
        if signal in source_line:
            for line in msg:
                periodic_insert(wrapcount, head, tail, line)
        else:
            print(source_line)