function [times, volts, trials] = hloadall(unit)
%function [times, volts, trials] = hloadall(unit)
%
% use dbfind & hload to find and load all the data associated with
% a single recording site/unit at once for template making.
%
% results are much like hload for a single pypefile, but NaNs are
% used to mark file boundaries
%

files = dbfind(unit, 'all', 'list');
times = [];
volts = [];
trials = [];

tt = [];
et = 0;
for n = 1:length(files)
  fprintf('%s\n', files{n});
  [t, v, tr, pf] = hload(files{n});
  times = [times t NaN];
  volts = [volts v NaN];
  trials = [trials; tr; [NaN NaN]];
  for k = 1:length(pf.rec)
    % pypefile's idea of trial start times
    tt = [tt sscanf(pf.rec(k).ev_e{2}, 'TOD_START %f')];
  end
  tt = [tt NaN];
end

subplot(3,1,1);
skip = 50;
plot(times(1:skip:end)-times(1), volts(1:skip:end), '.');
x = find(isnan(times));
for n = 1:length(x)
  vline(times(x(n)-1)-times(1));
end
axis tight;
title('tank time');

subplot(3,1,2);
plot(trials(:,1)-times(1), 0, 'go', trials(:,2)-times(1), 1, 'ro');
axis tight;

subplot(3,1,3);
plot(tt-tt(1), 0, 'k.');
x = find(isnan(tt));
for n = 1:length(x)
  vline(tt(x(n)-1)-tt(1));
end
axis tight;
title('pype time');


