#!/usr/bin/env python3

import sys
import loos

atomic_number = int(sys.argv[1])

ff = loos.FormFactor(atomic_number)

for i in range(40):
    q = 0.05 * i
    val = ff.compute(q)
    print(q, val)
