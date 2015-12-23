function showsnips(tank)

SHOW = 1000;

x = unique([tank.snips_ch tank.snips_sc], 'rows');

np = 1;
chs = unique(tank.snips_ch);
for ch = chs
  subplot(length(chs), 3, 3*np-2);
  scs = unique(tank.snips_sc(tank.snips_ch == ch))';
  for sc = scs
    ix = find(tank.snips_ch == ch & tank.snips_sc == sc);
    ix = ix(randperm(length(ix)));
    rix = ix(1:min(SHOW, length(ix)));
    plot(1:size(tank.snips_v, 2), ...
         1e6 * tank.snips_v(rix, :)', ...
         [pcolors(sc) '-']);
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
               1e6*mean(tank.snips_v(ix, :), 1), ...
               1e6*std(tank.snips_v(ix, :), 1));
    set(ls(2), 'Color', pcolors(sc));
    set(ls(1), 'FaceColor', pcolors(sc), 'FaceAlpha', 0.3);
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
    plot(centers(1:end-1)*1000, counts(1:end-1), [pcolors(sc) 'o-']);
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

