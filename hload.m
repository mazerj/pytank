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
%   channel - defaults to 1
%
% OUTPUTS
%   times - vector of times (secs)
%   volts - vector of continuous spike voltages
%   trials - Nx2 matrix of trial start (col1) and stop (col2) times
%
%   if isempty(times), then you asked for a segment past end of file.
%
% NOTES
%   if the trial start/stop boundaries span a segment border, you
%   need to be careful re-assembling the continuous data record!
%

if ~exist('channel', 'var')
  channel=1;
end

if ~exist('seg', 'var') || seg == 0
  times = [];
  volts = [];
  trials = [];
  seg = 1;
  while 1
    [t, v, tr] = hload(pf, seg, channel);
    if isempty(t)
      return
    else
      times = [times t];
      volts = [volts v];
      trials = [trials tr];
      seg = seg + 1;
    end
  end
end

blocks = list_blocks(pf);

hfiles = {};
for n = 1:length(blocks)
  b = strsplit(blocks{n}, '/');
  files = jls(sprintf('%s/%s-%s_???.th5', h5dump, b{end-1}, b{end}));
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
% this is a problem -- what happens when the segment splits
% a trial:
%   trials = [h5read(hf, '/hdr/tr_starts') h5read(hf, '/hdr/tr_stops')];

starts = h5read(hf, '/hdr/tr_starts');
stops = h5read(hf, '/hdr/tr_stops');
if length(starts) > length(stops)
  % stop event is in next segment:
  stops = [stops; NaN];
elseif length(starts) < length(stops)
  % start event is in previous segment:
  starts = [NaN; starts];
end
trials = [starts stops];

fprintf('%s (%.1fs)\n', hf, t(end)-t(1));

