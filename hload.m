function [times, volts, trials] = hload(pypefile)
%function [times, volts, trials] = hload(pypefile)
%
% Load hdf5 files containing tdt tank data matching specified
% pype file. This doesn't attempt to generate the hdf5 files
% automatically -- use tank2hdf5 (shell + python) to extract
% tanks to hdf5 files.
%
% INPUTS
%   pypefile - name of pypefile or pf struct from p2mLoad()
%
% OUTPUTS
%   times - vector of times (secs)
%   volts - vector of continuous spike voltages
%   trials - Nx2 matrix of trial start (col1) and stop (col2) times
%


if ischar(pypefile)
  pf = p2mLoad2(pypefile);
elseif isstruct(pypefile)
  pf = pypefile
else
  error('must provide pype filename or pypefile struct');
end
blocks = {};
for n = 1:length(pf.rec)
  t = pf.rec(n).params.tdt_tank;
  t = strrep(t, '\', '/');
  t = strrep(t, 'C:', '/auto/data/critters');
  t = [t '/' pf.rec(n).params.tdt_block];
  blocks{n} = t;
end
blocks = unique(blocks);
fprintf('%d block(s)\n', length(blocks));
times = [];
volts = [];
for n = 1:length(blocks)
  hfiles = jls([blocks{n} '?.hdf5']);
  for k = 1:length(hfiles)
    hf = hfiles{k}
    data = h5read(hf, '/continuous/spk');
    %h5disp(hf);
    if n == 1 && k == 1
      t0 = data(1,1);
    end
    times = [times data(1,:)];
    volts = [volts data(2,:)];
    fprintf('seg %d %c %f-%f\n', n, char('a'+k-1), data(1,1)-t0, data(1,end)-t0);
  end
end

trials = [h5read(hf, '/hdr/tr_starts') h5read(hf, '/hdr/tr_starts')];

