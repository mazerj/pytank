function [snips, events] = xsnips(y, a, b, spiketh, artifactth, randfrac)
%function [snips, events] = xsnips(y, a, b, spiketh, artifactth, randfrac)
% Extract snips from voltage trace.
% Find all snips that cross +/- spiketh, wihout exceeding the
%   +/- artifactth -- this is intended to reject motor noise from
%   off the bat.
% Snips are peak-aligned automatically
%
% INPUTS
%   y - 'voltage' trace (could be zscored etc.. doesn't matter)
%   a,b - pre and post cross windows (# samples) to snip out
%   spiketh - threshold voltage for spike
%   artifactth - threshold voltage for artifact
%   randfrac - random fraction (in addition to real snips)
%
% OUTPUTS
%   snips - Nsnips x sniplen matrix of snip voltage traces
%   events - index into y for the threshold crossing for each snip
%

if ~exist('randfrac', 'var')
  randfrac = 0;
end

ALIGN = 1;

events = find(diff(y > spiketh | y < -spiketh) > 0);

if randfrac > 0
  % also include a random sample of 'events' by dropping threshold
  r = find(diff(y > (spiketh/10) | y < -(spiketh/10)) > 0);
  %r = randperm(length(y));
  r = r(1:round(randfrac*length(events)));
  events = sort([events r]);
end

events = events((events > a) & (events < (length(y)-b)));
  
if ~isempty(events)
  snips = zeros([length(events) a+b+1]);
  for n = 1:length(events)
    snips(n, :) = y((events(n)-a):(events(n)+b));
  end
  
  % exclude any snips exceeding the artifact limit (motor noise..)
  ix = all(abs(snips) < artifactth,2);
  snips = snips(ix, :);
  events = events(ix);

  % exclude any snips with a Nan -- edge effect or it spans
  % a data file boundary.
  ix = all(~isnan(snips),2);
  snips = snips(ix, :);
  events = events(ix);

  if ALIGN
    m = max(abs(snips), [], 2);
    for n = 1:size(snips,1)
      k = median(find(abs(snips(n,:)) == m(n)));
      snips(n,:) = circshift(snips(n,:), a-k+1, 2);
      events(n) = events(n)- (a-k+1);
    end
  end
else
  snips = [];
end

