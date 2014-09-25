#!/usr/bin/env pypenv
# -*- Mode: Python; tab-width: 4; py-indent-offset: 4; -*-

import sys
import pylab as p
import pypedata as pd
from tank2hdf5 import Block, gettanks

if len(sys.argv) < 2:
    sys.stderr.write('usage: %s [^]pypefile\n' % sys.argv[0])
    sys.exit(1)
    
tankdir, blocklist = gettanks(sys.argv[1])
for b in blocklist:
    block = Block(tankdir, b)
    block.getall()
    block.plot()
