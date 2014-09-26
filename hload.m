function [times, volts, trials] = hload(pf, seg, channel)
%function [times, volts, trials] = hload(pf, seg, channel)
%
% Load hdf5 files containing tdt tank data matching specified
% single pype file. This doesn't attempt to generate the hdf5 files
% automatically -- use tank2hdf5 (shell + python) to extract
% tanks to hdf5 files first.
%
% INPUTS
%   pf - pypefile struct from p2mLoad()
%   seg - seg num to load
%
% OUTPUTS
%   times - vector of times (secs)
%   volts - vector of continuous spike voltages
%   trials - Nx2 matrix of trial start (col1) and stop (col2) times
%
%   if isempty(times), then you asked for a segment past end of file.
%
%  NOTE: this is hard coded to pull data from CHANNEL 1!!!
%    see data(1,:) line below..


H5DUMP='/auto/th5';

if ~exist('channel', 'var')
  channel=1;
end

blocks = tdtblocks(pf);

hfiles = {};
for n = 1:length(blocks)
  b = strsplit(blocks{n}, '/');
  files = jls(sprintf('%s/%s-%s_???.th5', H5DUMP, b{end-1}, b{end}));
  for k = 1:length(files)
    hfiles{length(hfiles)+1} = files{k};
  end
end

times = [];
volts = [];
trials = [];
if seg > length(hfiles)
  return
end

hf = hfiles{seg};
data = h5read(hf, '/continuous/spk');
tstart = h5readatt(hf, '/continuous/spk', 'tstart');
tend = h5readatt(hf, '/continuous/spk', 'tend');
t = linspace(tstart, tend, size(data, 2));
times = [times t];
volts = [volts data(channel,:)];
trials = [h5read(hf, '/hdr/tr_starts') h5read(hf, '/hdr/tr_stops')];

fprintf('%s (%.1fs)\n', hf, t(end)-t(1));

function blocks = tdtblocks(pf)
%function blocks = tdtblocks(pf)
%
% get list of tdt blocks (including name mangling from local
% storage target to raid location) containing analog data
% traces for specified pypefile.
%
% INPUTS
%  pf - pypefile struct from p2mLoad
%
% OUTPUTS
%  blocks - cell array list of block names
%

blocks = {};
for n = 1:length(pf.rec)
  t = pf.rec(n).params.tdt_tank;
  t = strrep(t, '\', '/');
  t = strrep(t, 'C:', '/auto/data/critters');
  t = [t '/' pf.rec(n).params.tdt_block];
  blocks{n} = t;
end
blocks = unique(blocks);
