function pyt_showsnips_te(tank, style, mode)
%function pyt_showsnips_te(tank, style, mode)
%
% Show on-line sorted snip shapes over time (temporal evolution)
% 
%INPUT
%  tank - data tank structure from pyt_load() -- only snip data are required!
%  style - 'h' or 'v' for emphasis on files vs sorts, respectively
%  mode - 'raw' or 'mean'
%
%  for 'h' style, all sorts in each data file are plotted by
%  channel for 'v' style all sorts across files are plotted on the
%  same axis to emphasize changes in spike shape over time.
%
%OUTPUT
%  none
%

if ~exist('style', 'var')
  style = 'v';
end

if ~exist('mode', 'var')
  mode = 'mean';
end

switch style(1)
  case 'h'
    plot_h(tank, mode);
  otherwise
    plot_v(tank, mode);
end

boxtitle(tank.exper);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_h(tank, mode);

MAXSHOW = 1000;

nsrcs = length(tank.srcs);
chs = unique(tank.snips_ch);

vrange = 10 * 1e6 * mean(std(tank.snips_v,[],2));
fs = 1.0 / mean(diff(tank.snips_t(1,:)));
tbase = 1000 * (1/fs) * ((0:size(tank.snips_v,2)-1));

% the time of threshold cross is hardcoded here -- 6 samples in:
ix = 6;
t0 = mean(tbase(ix));
tbase = tbase-t0;

clf;
np = 1;
for srcn = 1:nsrcs
  for ch = chs
    subplot(length(chs), nsrcs, np);
    scs = unique(tank.snips_sc(tank.snips_srcn == srcn & ...
                               tank.snips_ch == ch))';
    for sc = scs
      ix = find(tank.snips_srcn == srcn & ...
                tank.snips_ch == ch & tank.snips_sc == sc);
      switch mode
        case 'raw'
          ix = ix(randperm(length(ix)));
          rix = ix(1:min(MAXSHOW, length(ix)));
          plot(tbase, ...
               1e6 * tank.snips_v(rix, :)', ...
               [tdtsnipcolors(sc) '-']);
          set(gca, 'Color', [0.5 0.5 0.5]);
        case 'mean'
          ls = eplot(tbase, ...
                     1e6*mean(tank.snips_v(ix, :), 1), ...
                     1e6*std(tank.snips_v(ix, :), 1));
          set(ls(2), 'Color', tdtsnipcolors(sc));
          set(ls(1), 'FaceColor', tdtsnipcolors(sc), 'FaceAlpha', 0.3);
          set(gca, 'Color', [0.5 0.5 0.5]);
        otherwise
          error('mode must be raw or mean');
      end
      hold on;
    end
    title(sprintf('srcn=%d ch=%d', srcn, ch));
    axis square;
    axis tight;
    yrange(-vrange, +vrange);
    
    hline(0, 'LineStyle', '-', 'Color', 'w');
    vline(0, 'LineStyle', '-', 'Color', 'w');
    
  end
  
  grid on;
  
  xt = get(gca, 'XTick'); xt = xt(end-1);
  dx = get(gca, 'XTick'); dx = diff(dx(end-1:end));
  yt = get(gca, 'YTick'); yt = yt(end);
  for n = 0:5
    set(plot(xt+(n/6*dx), yt, 'ko'), ...
        'MarkerFaceColor', tdtsnipcolors(n));
  end
  
  hold off;
  
  if np == 1
    xlabel('time (ms)');
    ylabel('voltage (uv)');
  end
  np = np + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_v(tank, mode)

MAXSHOW = 100;

nsrcs = length(tank.srcs);
chs = unique(tank.snips_ch);

vrange = 10 * 1e6 * mean(std(tank.snips_v,[],2));
fs = 1.0 / mean(diff(tank.snips_t(1,:)));
tbase = 1000 * (1/fs) * ((0:size(tank.snips_v,2)-1));

% the time of threshold cross is hardcoded here -- 6 samples in:
ix = 6;
t0 = mean(tbase(ix));
tbase = tbase-t0;

ysep = .4*vrange;

clf;
np = 1;
for ch = chs
  scs = unique(tank.snips_sc(tank.snips_ch == ch))';
  for sc = scs
    subplot(length(chs), length(scs), np);

    for srcn = 1:nsrcs
      offset = -ysep * (srcn - 1);
      scs = unique(tank.snips_sc(tank.snips_srcn == srcn & ...
                                 tank.snips_ch == ch))';
      ix = find(tank.snips_srcn == srcn & ...
                tank.snips_ch == ch & tank.snips_sc == sc);
      switch mode
        case 'raw'
          ix = ix(randperm(length(ix)));
          rix = ix(1:min(MAXSHOW, length(ix)));
          plot(tbase, ...
               offset + 1e6 * tank.snips_v(rix, :)', ...
               [tdtsnipcolors(sc) '-']);
          set(gca, 'Color', [0.5 0.5 0.5]);
        case 'mean'
          ls = eplot(tbase, ...
                     offset + 1e6*mean(tank.snips_v(ix, :), 1), ...
                     1e6*std(tank.snips_v(ix, :), 1));
          set(ls(2), 'Color', tdtsnipcolors(sc));
          set(ls(1), 'FaceColor', tdtsnipcolors(sc), 'FaceAlpha', 0.3);
          set(gca, 'Color', [0.5 0.5 0.5]);
        otherwise
          error('mode must be raw or mean');
      end
      hold on;
      axis tight;
      hline(offset, 'LineStyle', '-', 'Color', 'w');
    end
    yrange(offset-ysep, ysep);
    
    title(sprintf('ch=%d sc=%d', ch, sc));
    grid on;
    vline(0, 'LineStyle', '-', 'Color', 'w');

    if np == 1
      xlabel('time (ms)');
      ylabel('[<---#999]   voltage (uv)    [#000-->]');
    end
    
    np = np + 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ls = eplot(x, y, ye)
%function ls = eplot(x, y, ye)
%
% Uses eshade to generate nice looking simple x,y,error plots in
% one function call.
%

if ~ishold
  cla;                                  % otherwise old eshade's persist..
end
h = ishold;
hold on;
ls = [eshade(x, y, ye); plot(x, y, 'r-')];
if ~h, hold off, end;




