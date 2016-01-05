function pyt_showsnips_te(tank, mode)
%function pyt_showsnips_te(tank, mode)
%
% Show on-line sorted snip shapes over time (temporal evolution)
% 
%INPUT
%  tank - data tank structure from pyt_load() -- only snip data are required!
%  mode - 'raw' or 'mean'
%
%OUTPUT
%  none
%

if ~exist('mode', 'var')
  mode = 'mean';
end

MAXSHOW = 1000;

nsrcs = length(tank.srcs);
chs = unique(tank.snips_ch);

vrange = 10 * 1e6 * mean(std(tank.snips_v,[],2));
fs = 1.0 / mean(diff(tank.snips_t(1,:)));
tbase = 1000 * (1/fs) * ((0:size(tank.snips_v,2)-1));

v = std(tank.snips_v,[],1);
ix = find(v == min(v));
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
        otherwise
          error('mode must be raw or mean');
      end
      hold on;
    end
    title(sprintf('srcn=%d ch=%d', srcn, ch));
    axis square;
    axis tight;
    yrange(-vrange, +vrange);
    
    hline(0, 'LineStyle', '-', 'Color', 'm');
    vline(0, 'LineStyle', '-', 'Color', 'm');
    
    if 0
      xt = get(gca, 'XTick'); xt = mean(xt(end-1:end-1));
      yt = get(gca, 'YTick'); yt = mean(yt(end-1:end));
      for n = 0:4
        set(plot(xt+2*n, yt, 'ko'), ...
            'MarkerFaceColor', tdtsnipcolors(n));
      end
    end
  end
  hold off;
  
  if np == 1
    xlabel('time (ms)');
    ylabel('voltage (uv)');
  end
  np = np + 1;
end

boxtitle(tank.exper);

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

