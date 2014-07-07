function [times, volts, trials] = hloadall(unit, validate, maxsegs)
%function [times, volts, trials] = hloadall(unit, validate, maxsegs)
%
% use dbfind & hload to find and load all the data associated with
% a single recording site/unit to use for global template making.
%
% results are much like hload for a single pypefile, but NaNs are
% used to mark file (not block!) boundaries.
%
% INPUTS
%   unit - string unit descriptor (eg, 'bert0321'; dbfind-friendly)
%   validate - if 1, then generate validation plots for timing etc
%   maxsegs - maximum number of segs from each block to load
%    (use this to control memory consumption)
%
% OUTPUTS
%    times - time vector (s)
%    volts - voltage vector (v)
%    trial - ntrials x 2 matrix of trial start/stop times (s)
%

if ~exist('validate', 'var')
  validate = 0;
end

if ~exist('maxsegs', 'var')
  maxsegs = +Inf;
end

files = dbfind(unit, 'all', 'list');
times = [];
volts = [];
trials = [];

tt = [];
et = 0;
fnames = {};
for n = 1:length(files)
  fprintf('%s\n', files{n});
  fnames{n}= basename(files{n});
  pf = p2mLoad2(files{n});
  [t, v, tr] = hload(pf, maxsegs);
  times = [times t NaN];
  volts = [volts v NaN];
  trials = [trials; tr; [NaN NaN]];
  if validate
    for k = 1:length(pf.rec)
      % pypefile's idea of trial start times
      tt = [tt sscanf(pf.rec(k).ev_e{2}, 'TOD_START %f')];
    end
    tt = [tt NaN];
  end
end

if validate
  subplot(3,3,[1 2]);
  skip = 50;
  plot(times(1:skip:end)-times(1), 1e6*volts(1:skip:end), '-');
  x = find(isnan(times));
  for n = 1:length(x)
    vline(times(x(n)-1)-times(1));
  end
  axis tight;
  title('tank time');
  ylabel('uv');

  subplot(3,3,[4 5]);
  plot(trials(:,1)-times(1), -0.5, 'go', trials(:,2)-times(1), 0.5, 'ro');
  axis tight;
  yrange(-1,1);
  title('tank trial bounds');
  
  subplot(3,3,[7 8]);
  plot(tt-tt(1), 0, 'k.');
  x = find(isnan(tt));
  for n = 1:length(x)
    vline(tt(x(n)-1)-tt(1));
  end
  axis tight;
  title('pype trial starts');
  
  subplot(3,3,[3 6 9]);
  cla;
  set(text(0, 1 , fnames), ...
      'HorizontalAlignment', 'left', ...
      'VerticalAlignment', 'top');
  axis off;
end


