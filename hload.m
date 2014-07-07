function [times, volts, trials] = hload(pf, maxsegs)
%function [times, volts, trials] = hload(pf, maxsegs)
%
% Load hdf5 files containing tdt tank data matching specified
% single pype file. This doesn't attempt to generate the hdf5 files
% automatically -- use tank2hdf5 (shell + python) to extract
% tanks to hdf5 files.
%
% INPUTS
%   pf - pypefile struct from p2mLoad()
%   maxsegs - maximum segments to load
%
% OUTPUTS
%   times - vector of times (secs)
%   volts - vector of continuous spike voltages
%   trials - Nx2 matrix of trial start (col1) and stop (col2) times
%   pf - Nx2 matrix of trial start (col1) and stop (col2) times
%

%  NOTE: this is hard coded to pull data from CHANNEL 1!!!
%    see data(1,:) line below..

CHANNEL=1;

if ~exist('maxsegs', 'var')
  maxsegs = +Inf;
end

blocks = tdtblocks(pf);
fprintf('%d Block(s)\n', length(blocks));

times = [];
volts = [];
for n = 1:length(blocks)
  hfiles = jls([blocks{n} '?.hdf5']);
  % try to skip files randomly to get down to maxsegs
  skip = 1;
  while 1
    hix = 1:skip:length(hfiles);
    if length(hix) <= maxsegs
      break;
    end
    skip = skip + 1;
  end
  for k = hix
    hf = hfiles{k};
    data = h5read(hf, '/continuous/spk');
    tstart = h5readatt(hf, '/continuous/spk', 'tstart');
    tend = h5readatt(hf, '/continuous/spk', 'tend');
    t = linspace(tstart, tend, size(data, 2));
    times = [times t];
    volts = [volts data(CHANNEL,:)];
    if skip > 1
      times = [times NaN];
      volts = [volts NaN];
    end
      
    if n == 1 && k == 1
      t0 = times(1);
    end
    fprintf('seg %d %c %f-%f\n', n, char('a'+k-1), t(1)-t0, t(end)-t0);
  end
  if length(hix) < length(hfiles)
    fprintf('** only loaded %d/%d segments\n', length(hix), length(hfiles));
  end
end

trials = [h5read(hf, '/hdr/tr_starts') h5read(hf, '/hdr/tr_stops')];


