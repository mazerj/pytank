#!/usr/bin/env pypenv
# -*- Mode: Python; tab-width: 4; py-indent-offset: 4; -*-

import time
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
#block.refindsnips(1, 1)

if 0:
    print 'writing hdf5...'
    ti = time.time()
    block.savehdf5('foo')
    ti = time.time() - ti
    print 'done:',ti,'secs'

if 0:
    print 'writing mat...'
    ti = time.time()
    block.savemat('foo')
    ti = time.time() - ti
    print 'done:',ti,'secs'

p.show()

