function sortit(pypefile, force)

LSE=1;                                  % use template-LSE as sort criterion

if ~exist('force', 'var')
  force = 0;
end

pf = dbfind(pypefile);
[templates, times, volts] = mktemplates(pf, force);
if isempty(times)
  [times, volts] = hload(pf.src);
end
t0 = times(1);

% test sort.. use a reduced the threshold to see what we
% misclassify.. sometimes it's good to be conservative when
% generating the templates to avoid getting overwhelmed by
% noise, then more permissive during the sort.
nsig = 10;

nunits = length(templates.units);
[snips, events] = xsnips(volts, templates.a, templates.b, ...
                         std(volts)*nsig, std(volts)*nsig*3);
snipt = times(events);

if LSE
  % minimal LSE
  scores = zeros([size(snips, 1) nunits]);
  for nu = 1:nunits
    for ns = 1:size(snips, 1)
      scores(ns, nu) = sum((templates.v(nu,:) - snips(ns, :)).^2);
    end
  end
  [~, sortcodes] = find(scores == repmat(min(scores, [], 2), ...
                                         [1 nunits]));
else
  scores = snips * templates.v';
  % negative projections are bad..
  scores(scores < 0) = -Inf;
  [~, sortcodes] = find(scores == repmat(max(scores, [], 2), [1 nunits]));
end

colors = 'rgbmyrgbmy';

for n = 1:nunits
  subplot(3, nunits, n);
  hist(scores(:,n));
  xrange(min(scores(:)),max(scores(:)));
  set(title(sprintf('%d: n=%d', templates.units(n), sum(sortcodes==n))), ...
            'Color', colors(n));
  if n == 1, ylabel('count'); end
  xlabel('score');

  subplot(3, nunits, n+nunits);
  eplot(templates.t, ...
        1e6*templates.v(n,:), 1e6*templates.ve(n,:));
  yrange(1e6*min(templates.v(:) - templates.ve(:)), ...
         1e6*max(templates.v(:) + templates.ve(:)));
  hline(0, 'linestyle', '-');
  vline(0, 'linestyle', '-');
  if n == 1, ylabel('uV'); end
  xlabel('time (s)');
end

subplot(3,1,3);
plot(times-t0, 1e6*volts, 'k');
xlabel('block time (s)');
ylabel('uV');
axis tight;
hold on;
t = templates.t;
for n = 1:length(events)
  set(plot(t + times(events(n)) - t0, 1e6*snips(n, :), ...
           [colors(sortcodes(n)) '-']), ...
      'linewidth', 2);
end
hold off;


