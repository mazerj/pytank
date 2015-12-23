#!/usr/bin/env pypenv
# -*- Mode: Python; tab-width: 4; py-indent-offset: 4; -*-
#
# Convert tdt tanks associated with specified pypefile into HDF5 data
# stores.  This is FAST!!!
#
# Note -- this really converts tank blocks, but when run as a program
# it takes a pype datafile and then tracks down all the tanks and
# blocks associated with the file and converts them all.
#
# Each block is split into chunks of 50,000 segments (~10 mins) to
# avoid running out of memory on typical desktops when extracting the
# raw signal traces. Each chunk gets it's own sequentially numbered
# HDF5 file to allow resequencing on load. Time base is continues
# across chunks. Header information is the same in all chunk files
# just to make things easy to load.
#
# This application is not at all CPU intensive, but it is memory
# intensive, so breaking it up into managable chunks prevents the
# machine from locking up when converting big blocks/tanks/runs
#
#
# HDF5 structure:
#   /hdr                                general header info
#   /hdr/dacq_fs_hz
#   /hdr/src
#   /hdr/tr_starts
#   /hdr/tr_stops
#   /continuous
#   /continuous/lfp                     voltages for continues LFP waveform
#     attr fs_hz: sampling rate
#     attr units: 'V' for volts
#     attr filters: (lp,hp)
#     attr tstart,tend: timestamps
#   /continuous/spk                     voltages for continuous spike waveform
#     attr fs_hz: sampling rate
#     attr units: 'V' for volts
#     attr filters: (lp,hp)
#     attr tstart,tend: timestamps
#   /snip                               tdt generated snips
#   /snip/ch                            dacq channel
#   /snip/sc                            sort code
#   /snip/t                             time (s)
#   /snip/v                             voltage (v)
#

import os
import struct
import string
import types
import time
import sys

MAXSEGS_PER_FILE = 50000
H5DUMP = '/auto/th5'
        

import numpy as np
import pylab as p
from scipy.signal import iirdesign, filtfilt

def icode(x):
	"""Convert 4byte event code (Snip, eNeu etc) into string and vice versa.
	"""
	if isinstance(x, types.StringType):
		# string -> code
		if x is '@@@@':
			return 0
		else:
			return np.sum(map(ord, x) * np.power(2, 8*np.arange(4)))
	else:
		# code -> string (this could be faster, but doesn't happen much)
		if x == 0:
			# this is internal to TDT, I think..
			return '@@@@'
		else:
			y = ''
			for n in range(4):
				y = y + chr(((0xff<<(8*n)) & x)>>(8*n))
			return y

class Block():
    lfpfilt = None
    spkfilt1 = None                       #hipass
    spkfilt2 = None                       #lowpass
    # from TDTbin2mat.c:
    formats = {
        0: (4, 'f'),                      # DFORM_FLOAT
		1: (4, 'i'),                      # DFORM_LONG
		2: (2, 'h'),                      # DFORM_SHORT
		3: (1, 'b'),                      # DFORM_BYTE
		4: (8, 'd'),                      # DFORM_DOUBLE
		5: (8, 'Q'),                      # DFORM_QWORD
        }
    
	def __init__(self, tank, block,
                 lfpcut=200.0, spikecut=400.0, spiketop=8000.0):
		"""TDT block (one file inside tank). A pypefile could reference
        multiple blocks if it has been appended to! The block index (.TSQ)
        is read into memory and parsed on init. Index is used by get*
        methods to retreive raw data streams from .SEV and .TEV files.
		"""
		if tank[-1] == os.sep:
			tank = tank[:-1]
		tankdir, tankname = os.path.split(tank)
        self.src = tank+os.sep+block
        
        if not os.path.exists(self.src):
            sys.stderr.write('%s missing\n' % self.src)
            sys.exit(1)

        h = '%s_%s' % (tankname, block)
        base = string.join((tankdir, tankname, block, h), os.sep)
		self.tsqfile = base + '.tsq';
		self.tevfile = base + '.tev';
		self.fs = None
        self.lfpcut = lfpcut
        self.spikecut = spikecut
        self.spiketop = spiketop

		# TSQ record format is:
		#	0 size		  int32
		#	1 type		  int32
		#	2 icode		  int32
		#	3 channel	  uint16
		#	4 sortcode	  uint16
		#	5 timestamp	  double
		#	6 offset	  uint64  <--+ note: these are the same field in record,
		#  *6 strobe	  double  <--|		  but how to parse depends on icode
		#	7 format	  int32
		#	8 frequency	  single float
		
		recfmt = 'iiiHHdQif'
		s = struct.Struct(recfmt)
		recfmt = 'iiiHHddif'
		s2 = struct.Struct(recfmt)

		buf = open(self.tsqfile, 'rb').read()
		nrec = (len(buf)/s.size) - 2
		
		self.size = [0] * nrec
		self.type = [0] * nrec
		self.icode = [0] * nrec
		self.channel = [0] * nrec
		self.sortcode = [0] * nrec
		self.timestamp = [0] * nrec
		self.offset = [0] * nrec
		self.format = [0] * nrec
		self.frequency = [0] * nrec
		self.strobe = [0] * nrec

		# skip records 0 & 1 (tdt-internal info/state)
		for n in range(2, nrec):
            ix = n - 2
            
			# parsing each record twice is not very efficient.. but
			# this is not the speed bottle neck	 -- it's the data
			# loading part.
			a = s.unpack_from(buf, n*s.size)
			b = s2.unpack_from(buf, n*s2.size)

			self.size[ix] =      a[0]
			self.type[ix] =      a[1]
			self.icode[ix] =     a[2]
			self.channel[ix] =   a[3]
			self.sortcode[ix] =  a[4]
			self.timestamp[ix] = a[5]
			self.offset[ix] =    a[6]     # alternative 1
			self.strobe[ix] =    b[6]     # alternative 2
			self.format[ix] =    a[7]
			self.frequency[ix] = a[8]
			
        # convert them all to numpy arrays for speed
		self.size = np.array(self.size)
		self.type = np.array(self.type)
		self.icode = np.array(self.icode)
		self.channel = np.array(self.channel)	# 1,2..nchan (not zero based!)
		self.sortcode = np.array(self.sortcode)
		self.timestamp = np.array(self.timestamp)
		self.offset = np.array(self.offset)
		self.format = np.array(self.format)
		self.frequency = np.array(self.frequency)
		self.strobe = np.array(self.strobe)
		self.recno = np.arange(len(self.size))   # record number
        self.nrec = nrec

	def printsummary(self):
		for c in np.unique(self.icode):
			n = np.sum(self.icode == c)
			print '%4s %d' % (icode(c), n,)
		print 'channels: ', self.getchns()

    def splits(self):
        """Split block up into chunks of up to MAXSEGS_PER_FILE records.
        This is to avoid running out of memory during extraction or extraction.
        """
        s = []
        a = 0
        n = 0
        while a < self.nrec:
            b = min(a + MAXSEGS_PER_FILE, self.nrec)
            s.append((('_%03d'%n),a,b,))
            a = b
            n = n + 1
        return s
            
			
	def _getsegments(self, seglist, snips=False):
		"""Retrieve indicate segment numbers the .TEV or .SEV files
		based on data loaded from the index loaded from the .TSQ file.
		Note: all segements should be of same icode/type and length
		Returns: time, voltage (2d arrays - nsegments x segmentlen)
		"""
		
		# v: voltage (volts)
		# t: time (secs)
		v, t = [], []
		if snips:
			sortcodes = []
		f, fname = None, None
		for n in seglist:
			sevfile = \
			  string.replace(self.tevfile, '.tev',
							 '_%s_Ch%d.sev' % (icode(self.icode[n]),
											   self.channel[n]))
			if os.path.exists(sevfile):
				if fname is not sevfile:
					if f:
                        f.close()
					f = open(sevfile, 'r')
					fname = sevfile
			else:
				if fname is not self.tevfile:
					if f:
                        f.close()
					f = open(self.tevfile, 'r')
					fname = self.tevfile
				
			f.seek(self.offset[n], 0)                    # from BOF
			nlongs = self.size[n] - 10					 # 4byte units
			buf = f.read(4 * nlongs)
            sz, code = self.formats[self.format[n]]
            nsamp = nlongs * 4 / sz
			v.append(np.array(struct.unpack(code * nsamp, buf)))
            tt = np.arange(nsamp) / self.frequency[n] + self.timestamp[n]
			if snips:
				# Snip timestamps indicate time of 1st thresh crossing,
				# which is 1/4 way into the trace, subtract this out
				# this offset so snip time is veridical.
				tt = tt - ((nsamp + 2) / 4) / self.frequency[n];
				sortcodes.append(self.sortcode[n])
			t.append(tt)
            
		if f:
            f.close()
            
        # convert list of arrays into 2d array
		v = np.array(v)
		t = np.array(t)
		if snips:
            return t, v, np.array(sortcodes)
        else:
            return t, v

	def getchns(self):
		"""Get a list of all channels active in this block
		Note: channels start with 1!!!
		"""
		mask = 0
		for m in np.unique(self.strobe[np.where(self.icode==icode('CHNS'))[0]]):
			mask |= int(m)
		channels = []
		for n in range(32):
			if mask & (1<<n):
				channels.append(n+1)
		return channels

	def getsnips(self, channel, first=0, last=+np.inf):
		"""Get snip data for specified channel.
		Returns: time, voltage (2d arrays - nsnips x sniplen)
		"""
		ix = np.where(((self.icode == icode('eNeu')) | \
                       (self.icode == icode('Snip'))) &
                       (self.channel == channel) &
                       (self.recno >= first) &
                       (self.recno < last))[0]
		return self._getsegments(ix, snips=True)

	def getraw(self, channel, first=-np.inf, last=+np.inf):
		"""Get raw 16bit voltage trace for channel over entire block.
		Returns: time, voltage (1d vectors)
		"""
		ix = np.where((self.icode==icode('RAW0')) &
                      (self.channel==channel) &
                      (self.recno >= first) &
                      (self.recno < last))[0]
		t, v = self._getsegments(ix)
		return t.flatten(), v.flatten()

	def gettrials(self):
		"""Get start and stop time of each trial in this block.
		"""
		trl1 = self.timestamp[np.where(self.icode==icode('TRL1'))[0]]
		trl2 = self.timestamp[np.where(self.icode==icode('TRL2'))[0]]
		if trl1.size > trl2.size:
			# discard last trl1 event if no matching trl2 event
			trl1 = trl1[:-1]
		return trl1, trl2

	def W(self, hz):
		"""Normalized frequency (1=nyquist).
		"""
		return hz / (self.fs / 2.0)

	def getall(self, first=-np.inf, last=+np.inf):
		"""Pull all key from the tank into memory.
		This includes generating spike and lfp waveforms by filtering
		the RAW0 signal.
		"""
		
		self.channellist = self.getchns()
		self.start, self.stop = self.gettrials()
		if len(self.channellist) < 1:
			return False

		
		self.raw = {}
		for c in self.channellist:
			self.raw[c] = self.getraw(c, first=first, last=last)
			if self.fs is None:
				self.fs = (1.0/np.diff(self.raw[c][0][0:2]))[0]

		if self.lfpfilt is None:
			self.lfpfilt = iirdesign(self.W(self.lfpcut), \
                                     self.W(1.2*self.lfpcut), \
									 0.1, 45.0)
			self.spkfilt1 = iirdesign(self.W(self.spikecut), \
                                      self.W(0.8*self.spikecut), \
                                      0.1, 45.0)
			self.spkfilt2 = iirdesign(self.W(self.spiketop), \
                                      self.W(1.2*self.spiketop), \
                                      0.1, 45.0)
									 
		self.lfp = {}
		self.spk = {}
		for c in self.channellist:
			t, v = self.raw[c]
			# lowpass and decimate for lfp
			l = filtfilt(self.lfpfilt[0], self.lfpfilt[1], v)
			ds = np.ceil(self.fs/self.lfpcut/4.0)
			l = l[::ds]
			self.lfp[c] = (np.linspace(t[0], t[-1], len(l)), l)
            s = filtfilt(self.spkfilt1[0], self.spkfilt1[1], v)
            s = filtfilt(self.spkfilt2[0], self.spkfilt2[1], s)
			self.spk[c] = (t, s)

		self.snips = {}
		for c in self.channellist:
			self.snips[c] = self.getsnips(c, first=first, last=last)

		return True

	def showfilters(self):
		"""Plot frequency response for lfp and spikefilters.
		"""
		p.figure()
		w, h = freqz(self.lfpfilt[0],self.lfpfilt[1])
		p.plot(w, 20*np.log10(abs(h)), 'r-')
		w, h = freqz(self.spkfilt[0],self.spkfilt[1])
		p.plot(w, 20*np.log10(abs(h)), 'b-')
		p.xscale('log')
		p.title('lfp/spike freq response')
			
    def meansnip(self, channel, sortcode):
        t, v, sc = self.snips[channel]
        ix = np.where(sc==sortcode)[0]
        snip = np.mean(v[ix,:], 0)
        return snip

    def savehdf5(self, fname, force=False):
        """hdf5 structure:

        /hdr/src
        /hdr/dacq_fs_hz
        /hdr/channellist
        /hdr/tr_starts
        /hdr/tr_stops
        /hdr/tr_stops
        /continuous/RAW0  <- not included to save space!
        /continuous/spk
        /continuous/lfp
        /snip/t   time base
        /snip/v   voltage traces
        /snip/ch  tdt dacq channel
        /snip/sc  tdt sortcodes

        """
        
        import h5py
        
        # don't overwrite existing files unless force==True
        if os.path.exists(fname) and not force:
            sys.stderr.write('%s already exists\n' % fname)
            return
        f = h5py.File(fname, 'w')
        h = f.create_dataset('hdr/src', data=self.src)
        h = f.create_dataset('hdr/dacq_fs_hz', data=(self.fs,))
        h = f.create_dataset('hdr/channellist', data=self.channellist)
        
        h = f.create_dataset('hdr/tr_starts', data=self.start,
                             compression='gzip')
        h.attrs['units'] = 's'
        h = f.create_dataset('hdr/tr_stops', data=self.stop,
                             compression='gzip')
        h.attrs['units'] = 's'

        for k in range(1,3):
            if k == 0:
                d,name,filt = self.raw, 'RAW0', (-1,-1)
            elif k == 1:
                d,name,filt = self.spk, 'spk', (self.spikecut, self.spiketop)
            elif k == 2:
                d,name,filt = self.lfp, 'lfp', (-1, self.lfpcut)

            # save matrix c1,c2...cN x #samples:
            # need to save timebase for each signal since there's not
            # reason to assume they're all at same sampling rate..
            # times can be calculated from the 'tzero' and 'tend' attributes
            m = np.zeros([d[self.channellist[0]][0].shape[0],
                          len(self.channellist)])
            t = d[self.channellist[0]][0]   # sample times
            fs = 1.0/(t[1]-t[0])
            n = 0
            for c in self.channellist:
                m[:,n] = d[c][1]
                n = n + 1
            h = f.create_dataset('continuous/%s' % name, data=m,
                                 compression='gzip')
            h.attrs['fs_hz'] = fs
            h.attrs['units'] = 'V'
            h.attrs['filters'] = '%s' % (filt,)
            h.attrs['tstart'] = t[0]
            h.attrs['tend'] = t[-1]
                    
        for c in self.channellist:
            t, v, sc = self.snips[c]
            ch = c + (0 * sc)
            if c == self.channellist[0]:
                ts, vs, chs, scs = t, v, ch, sc
            else:
                ts = np.concatenate((ts, t))
                vs = np.concatenate((vs, v))
                chs = np.concatenate((chs, ch))
                scs = np.concatenate((svs, sv))
                
        h = f.create_dataset('snip/t', data=ts, compression='gzip')
        h.attrs['units'] = 's'
        h = f.create_dataset('snip/v', data=vs, compression='gzip')
        h.attrs['fs_hz'] = self.fs
        h.attrs['units'] = 'V'
        h = f.create_dataset('snip/ch', data=chs, compression='gzip')
        h.attrs['units'] = '#1-N'
        h = f.create_dataset('snip/sc', data=scs, compression='gzip')
        h.attrs['units'] = '#0-N'

        f.close()
                
	def plot(self):
        import pylab as p
        
        for c in self.channellist:
            colors = 'krgbymrgbymrgbym'
            
            p.figure(figsize=(8, 10))
            p.clf()
            p.subplot(2,2,3)
			t, v, sc = self.snips[c]
			for n in range(v.shape[0]):
				p.plot(1000*(t[n,:]-t[n,0]), 1e6*v[n,:],
                       '%s-' % colors[sc[n]])
            p.xlabel('time (ms)')
            p.ylabel('voltage (uV)')
            p.title('channel=%d' % c)

            p.subplot(2,2,4)
			t, v, sc = self.snips[c]
            for code in np.unique(sc):
                ix = np.where(sc==code)[0]
                x = 1000*(t[n,:]-t[n,0])
                y = 1e6 * np.mean(v[ix,:],0)
                e = 1e6 * np.std(v[ix,:],0)
                #e = np.std(v[ix,:],0) / sqrt(len(ix))
                p.plot(x, y, '%s-' % colors[code])
                p.fill_between(x, y-e, y+e, alpha=0.2, facecolor=colors[code])
            p.xlabel('time (ms)')
            p.autoscale(axis='x', tight=True)

            p.subplot(6,1,1)
            t, v = self.raw[c]
            p.plot(t-t[0], v*1e6, 'k-')
            p.ylabel('raw (uV)')
            p.autoscale(axis='x', tight=True)
            p.title(self.src)
            
            p.subplot(6,1,2)
            t, v = self.lfp[c]
            p.plot(t-t[0], v*1e6, 'r-')
            p.autoscale(axis='x', tight=True)
            p.ylabel('lfp (uV)')
            
            p.subplot(6,1,3)
            t, v = self.spk[c]
            p.plot(t-t[0], v*1e6, 'b-')
            p.ylabel('spikes (uV)')
            p.xlabel('time (s)')
            p.autoscale(axis='x', tight=True)
            
def getblocks(pypefile):
    """Get tankdir and list of all blocks references in pypefile.
    """
    import pypedata as pd

    pf = pd.PypeFile(pypefile)
    n = 0
    blocks = {}
    while 1:
        rec = pf.nth(n)
        if rec is None:
            break
        if n == 0:
            tankdir = string.replace(rec.params['tdt_tank'], \
                                     'C:\\DataTanks\\', \
                                     '/auto/data/critters/DataTanks/')
        blocks[rec.params['tdt_block']] = 1
        n += 1
    return tankdir, blocks.keys()
        
if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.stderr.write('usage: %s [-force] ..pypefiles..\n' % sys.argv[0])
        sys.exit(1)

    force = False
    files = []
    for f in sys.argv[1:]:
        if f == '-force':
            force = True
        else:
            files.append(f)

    for f in files:
        if f.find('0000') > 0:
            # always skip 'training' files with no physiology data
            continue
        tankdir, blocklist = getblocks(f)
        for b in blocklist:
            print '%s -->' % f
            block = Block(tankdir, b)
            for (k, a, b) in block.splits():
                if H5DUMP:
                    base = string.join(block.src.split('/')[-2:], '-')
                    outfile = '%s/%s%s.th5' % (H5DUMP, base, k,)
                else:
                    outfile = '%s%s.th5' % (block.src, k,)
                if block.getall(first=a, last=b):
				    print '  %s' % (outfile,)
					try:
						block.savehdf5('%s' % (outfile,), force=force)
					except IOError:
						sys.stderr.write('Can''t write: %s\n' % outfile)
						sys.exit(1)
				else:
				    print '  no trials.'
