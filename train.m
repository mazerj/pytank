function templates = train(times, volts)
%function templates = train(times, volts)
%
% generate templates interactively from spike data
%
%  left-click raw snip plots to select/unselect cluster(s) in lower
%  left snip plot window
%    'm' to merge all selected clusters into one
%    'r' to reset sorts (ie, throw away all merges)
%    'x' exit without saving
%    'esc' save templates and terminate interactive mode
%
% INPUTS
%   times - time stamps in s
%   volts - voltage trace in volts
%
% OUTPUTS
%   templates - structure containing templates
%     - templates.units - vector of unit ids (small ints)
%     - templates.t - 1 x snipsize matrix of time bases for templatex
%     - templates.v - nunits x snipsize matrix of voltage trace for templates
%     - templates.ve - nunits x snipsize voltage trace stdevs
%       NOTE: Only .v is really needed for sorting. The other stuff
%             is just useful to have around.
%

a = 5;                                  % nsamps before
b = 15;                                 % nsamps after for snip
nsig = 20;                              % default snip threshold
nart = 3;                               % artifact th = nart * nsig
npcs = 3;                               % use first n PCs for clustering
nc = 0;                                 % number of clusters

fs = 1.0 / (times(2)-times(1));

tvs = (-a:b)/fs;                        % snip time vector in secs
tvms = 1000*tvs;                        % snip time vector in ms
vmean = nanmean(volts);
vstd = nanstd(volts);
vnorm = (volts - vmean) ./ vstd;

oldfig = get(0, 'CurrentFigure');
h = figure;
set(gcf, 'MenuBar', 'none', 'ToolBar', 'none');

while 1
  oldpointer = get(gcf, 'pointer');
  set(gcf, 'pointer', 'watch'); drawnow;

  artifact = nsig * nart;               % anything exceeding this is artifact!
  x = [(NaN * zeros([1 2*a])) vnorm (NaN * zeros([1 2*b]))];
  t = [(NaN * zeros([1 2*a])) times (NaN * zeros([1 2*b]))];
  t = t - min(t);
  [snips, events] = xsnips(x, a, b, nsig, artifact, 0.20);
  
  clf;
  subplot(2,1,1);
  skip = 10;
  plot(t(1:skip:end), x(1:skip:end), '-');
  if ~isempty(events)
    hold on;
    plot(t(events+1), x(events+1), 'r.');
    hold off;
  end
  hmarks(nsig, artifact);
  title(sprintf('nsig=%.1f nsnips=%d', nsig, size(snips,1)));
  
  subplot(2,3,4);
  cla;
  if ~isempty(snips)
    plot(tvms, snips');
    hmarks(nsig, artifact);
    title(sprintf('n=%d', size(snips,1)));
  end
    
  subplot(2,3,5);
  cla;
  if ~isempty(snips)
    eshade(tvms, nanmean(snips,1), nanstd(snips,1));
    hold on;
    plot(tvms, nanmean(snips,1), 'r-');
    hold off;
    hmarks(nsig, artifact);
  end
    
  subplot(2,3,6);
  cla;
  if ~isempty(snips)
    dt = diff(t(events)) * 1000;
    hist(dt(dt < 25), 0:25);
    ylabel('freq');
    xlabel('isi (ms)');
    axis tight;
  end
  
  set(gcf, 'pointer', oldpointer); drawnow;
  
  subplot(2,3,5);
  [xx, yy, bb] = ginput(1);
  if bb == 27, break; end
  if bb == 1, nsig = abs(yy); end
end

metrics = [];

% pcs are "features" ((a+b) x (a+b), assuming (a+b) > nsnips..)
% scores are the projections of each snip onto each feature
% latent is the eigenvalues
[pcs, scores, latent] = pca(snips);
metrics = [metrics scores(:,1:npcs)];

% add a combination of hand-picked spike metrics (spike
% heights, widths etc) and first n PCs as feature metrics
% for clustering.

if 0
% note: this doesn't quite work right for monophasic spikes...
% the tmin (or tmax) is not really correct in that case,
% really need to do something like t > 1sigma..
vmax = max(snips, [], 2);
[~, tmax] = find(snips == repmat(vmax, [1 size(snips,2)]));
vmin = min(snips, [], 2);
[~, tmin] = find(snips == repmat(vmin, [1 size(snips,2)]));
metrics = [metrics vmin vmax tmin tmax tmax-tmin];
end

metrics = [metrics snips(:,a)];                 % voltage at 'peak'
metrics = [metrics mean(snips,2)];              % zscore is global,
                                                % this is local..
metrics = [metrics nansum(snips,2)];            % total area under snip
metrics = [metrics nansum(max(0, snips), 2)];   % positive area
metrics = [metrics nansum(min(0, snips), 2)];   % negative area

% calculate best estimate for number of clusters
eva = evalclusters(metrics, 'kmeans', 'CalinskiHarabasz', ...
                   'KList', [1 2 3 4 5 10 20]);
subplot(4,2,1);
plot(eva);

if nc == 0
  nc = eva.OptimalK;
end

c = [];
while 1
  if isempty(c)
    % now actually cluster based on 'optimal' number
    [c, ~, within, ~] = kmeans(metrics, nc, 'replicates', 10);
    
    % sort clusters from tight to loose
    c2 = zeros(size(c));
    [~, ix] = sort(within);
    for n = 1:length(within)
      c2(find(c == ix(n))) = 1+length(within)-n;
    end
    c = c2;
  end
    
  cs = unique(c);
  sp = [];
  clear templates
  templates.units = cs;
  for cn = 1:length(cs)
    subplot(4, length(cs), length(cs)+cn);
    plot(tvms, snips(c == cs(cn), :)');
    sp = [sp gca];
    
    axis tight;
    yrange(min(snips(:)), max(snips(:)));
    if cn == 1, ylabel('zscore'); end
    vline(0, 'linestyle', '-');
    title(sprintf('#%d n=%d', cs(cn), sum(c==cs(cn))));
    
    subplot(4, length(cs), 2*length(cs)+cn);
    templates.t{cs(cn)} = tvs;
    templates.v{cs(cn)} = (nanmean(snips(c == cs(cn), :)) * vstd) + vmean;
    templates.ve{cs(cn)} = (nanstd(snips(c == cs(cn), :)) * vstd);
    eplot(tvms, 1e6*templates.v{cs(cn)}, 1e6*templates.ve{cs(cn)});
    axis tight;
    yrange(1e6*vstd*min(snips(:)), 1e6*vstd*max(snips(:)));
    if cn == 1, ylabel('uV'); end
    vline(0, 'linestyle', '-');
    
    subplot(4, length(cs), 3*length(cs)+cn);
    dt = diff(t(events(c == cs(cn)))) * 1000;
    hist(dt(dt <= 25), 0:25);
    title(sprintf('%d >25ms', sum(dt > 25)));
    if cn == 1, ylabel('freq'); end
    xlabel('isi (ms)');
    axis tight;
  end
  selbox(sp);
  
  set(gcf, 'KeyPressFcn', {@keypress, h});
  set(gcf, 'UserData', []);
  waitfor(gcf, 'UserData');
  ch = get(gcf, 'UserData');
  switch ch
    case {char(27), 'q', 'x'}           % done/exit
      break
    case 'm'                            % merge
      v = get(sp(1), 'UserData');
      merge = find(v==1);
      for n = 2:length(merge)
        c(c == merge(n)) = merge(1);
        fprintf('merge %d->%d\n', merge(n), merge(1));
      end
      newcs = unique(c);
      for n=1:length(newcs)
        if newcs(n) ~= n
          c(c == newcs(n)) = n;
          fprintf('changed %d->%d\n', newcs(n), n);
        end
      end
    case 'r'
      c = [];
    otherwise
      ;
  end
end
close(gcf);
if ~isempty(oldfig)
  figure(oldfig);
end

%templates.site = pf.rec(1).params.cell;
templates.t = templates.t{1};
templates.v = cell2mat(templates.v');
templates.ve = cell2mat(templates.ve');
templates.a = a;
templates.b = b;

templates.nsig = nsig;                  % spike threshold (units of sigma)
templates.art = nsig * nart;            % artifact thresh (units of sigma)

if 0
  templates.v(size(templates.v,1)+1,:) = nanmean(volts);
  templates.ve(size(templates.ve,1)+1,:) = nanstd(volts);
  templates.units = [templates.units; length(templates.units)+1]; 
end

if 1
  for n = 1:size(metrics, 2)
    for k = 1:n
      subplot(size(metrics, 2), size(metrics, 2), (n-1)*size(metrics,2)+k);
      for cn = templates.units'
        ix = find(c == cn);
        plot(metrics(ix,n), metrics(ix,k), [pcolors(cn) '.']);
        hold on;
        axis equal
        axis square
        if n == size(metrics,2)
          xlabel(k);
        end
        if k == 1
          ylabel(n);
        end
      end
      hold off;
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selbox(axs)

set(axs(1), 'UserData', zeros(size(axs)));

n = 1;
for ax = axs
  invisax = axes('Position',get(ax,'OuterPosition'), ...
                 'PlotBoxAspectRatioMode',get(ax,'PlotBoxAspectRatioMode'),...
                 'DataAspectRatioMode',get(ax,'DataAspectRatioMode'),...
                 'NextPlot','add');
  axis off;
  p = patch([0 0 1 1], [0 1 1 0], [0 0 0 0], 'w');
  set(p,'FaceAlpha',0);
  set(p,'EdgeAlpha',0);
  
  setappdata(p, 'PlotNumber', n); n = n + 1;
  set(p,'ButtonDownFcn', {@toggleplot, axs(1)});
end

function toggleplot(src, eventdata, ax)
pn = getappdata(src,'PlotNumber');
a = get(src,'FaceAlpha');                % 0 when selected..
set(src,'FaceAlpha',0.25 - a);

v = get(ax, 'UserData');
v(pn) = ~a;
set(ax, 'UserData', v);

function keypress(src, eventdata, h)
set(h, 'UserData', eventdata.Character)


function hmarks(nsig, artifact)
for v = [-nsig nsig 0 -artifact artifact]
  hline(v, 'linestyle', '-', 'color', 'k');
end
