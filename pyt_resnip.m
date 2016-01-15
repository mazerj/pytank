function t = pyt_resnip(t)
%function t = pyt_resnip(t)
%
%  Goes back to 'raw' highpass spike trace and extract each snip
%  from the original minimally filtered highpass signal. The result
%  is a new pyt tank structure with the new snip waveforms. Note
%  that there are occasional edge cases where a snip goes off the
%  end of the available highpass data, in which case the snip is
%  all NaNs.
%
%INPUT
%  t - pyt tank structure from pyt_load()
%
%OUTPUT
%  t - new tank structure
%

chsc = unique([t.snips_ch t.snips_sc], 'rows');

% find index into t.highpass and t.highpass_t for each snip time
% this is the main optimization step -- this generates a lookup
% table that contains the mapping from snip times to indices in
% the highpass datastream but serial search of the highpass_t
% vector.
[sniptimes, rawixs] = findsnips_(t);

newsnips = NaN * t.snips_v;

OFF=5;                                  % offset -- depends on snipsize/circuit
offv = (1:size(t.snips_v, 2))-OFF;      % offset vector (for speed)

for ns = 1:size(chsc,1)
  ch = chsc(ns, 1);
  sc = chsc(ns, 2);
  ixs = find(t.snips_ch == ch & t.snips_sc == sc);
  for n = 1:length(ixs)
    t0 = t.snips_t(ixs(n), 1);
    vix = rawixs(sniptimes == t0);
    try
      % don't bother trying to find edge cases (for speed)
      newsnips(ixs(n), :) = t.highpass(vix+offv, ch);
    end
  end
end
if any(isnan(newsnips(:,1)))
  warning('resnip: %.1f%% edge snips', ...
          100 * sum(isnan(newsnips(:,1))) / size(newsnips, 1));
end
t.snips_v = newsnips;

function [times, ixs] = findsnips_(t)
% generate fast lookup table for mapping between sniptimes and
% highpass indices.

times = unique(t.snips_t(:,1));
ixs = zeros(size(times));

n = 1;
k = 1;
while k <= length(times) && n <= length(t.highpass_t)
  if t.highpass_t(n) >= times(k)
    ixs(k) = n;
    k = k + 1;
  end
  n = n + 1;
end
