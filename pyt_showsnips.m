function pyt_showsnips(tank)
%function pyt_showsnips(tank)
%
% Show on-line sorted snips tdt datasets. This generates
% per-channel plots of average snip shapes, isi's etc..
% 
%INPUT
%  tank - data tank structure from pyt_load() -- only snip data are required!
%
%OUTPUT
%  none
%

MAXSHOW = 1000;

np = 1;
chs = unique(tank.snips_ch);
for ch = chs
  subplot(length(chs), 3, 3*np-2);
  scs = unique(tank.snips_sc(tank.snips_ch == ch))';
  for sc = scs
    ix = find(tank.snips_ch == ch & tank.snips_sc == sc);
    ix = ix(randperm(length(ix)));
    rix = ix(1:min(MAXSHOW, length(ix)));
    plot(1:size(tank.snips_v, 2), ...
         1e6 * tank.snips_v(rix, :)', ...
         [tdtsnipcolors(sc) '-']);
    hold on;
  end
  title(sprintf('raw ch=%d', ch));
  hline(0, 'LineStyle', '-', 'Color', 'w');
  hold off;
  ax = axis;
  xlabel('sample #');
  ylabel('volt (uv)');
  set(gca, 'Color', [0.5 0.5 0.5]);
  axis square;
  
  subplot(length(chs), 3, 3*np-1);
  scs = unique(tank.snips_sc(tank.snips_ch == ch))';
  for sc = scs
    ix = find(tank.snips_ch == ch & tank.snips_sc == sc);
    ls = eplot(1:size(tank.snips_v, 2), ....
               1e6*nanmean(tank.snips_v(ix, :), 1), ...
               1e6*nanstd(tank.snips_v(ix, :), 1));
    set(ls(2), 'Color', tdtsnipcolors(sc));
    set(ls(1), 'FaceColor', tdtsnipcolors(sc), 'FaceAlpha', 0.3);
    hold on;
  end
  title(sprintf('mean{\\pm}std ch=%d', ch));
  xlabel('sample #');
  axis(ax);
  hline(0, 'LineStyle', '-', 'Color', 'w');
  hold off;  
  set(gca, 'Color', [0.5 0.5 0.5]);
  axis square;
  
  subplot(length(chs), 3, 3*np-0);
  scs = unique(tank.snips_sc(tank.snips_ch == ch))';
  bins = (0:50)/1000;
  centers = bins + (bins(2)-bins(1))/2.0;
  for sc = scs
    ix = find(tank.snips_ch == ch & tank.snips_sc == sc);
    counts = hist(diff(tank.snips_t(ix, 1)), bins);
    counts = counts ./ length(diff(tank.snips_t(ix, 1)));
    plot(centers(1:end-1)*1000, counts(1:end-1), [tdtsnipcolors(sc) 'o-']);
    hold on;
  end
  xlabel('isi (ms)');
  ylabel('proportion');
  title(sprintf('isi hist ch=%d', ch));
  hold off;  
  set(gca, 'Color', [0.5 0.5 0.5]);
  axis square;
  
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

