#!/usr/bin/env python3

import sys
import loos

filename = sys.argv[1]

system = loos.createSystem(filename)
print (system.totalCharge())
