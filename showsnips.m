function showsnips(tank)

SHOW = 1000;

x = unique([tank.snips_ch tank.snips_sc], 'rows');

np = 1;
chs = unique(tank.snips_ch);
for ch = chs
  subplot(length(chs), 1, np); np = np + 1;
  ns = 1;
  scs = unique(tank.snips_sc(tank.snips_ch == ch))'
  for sc = scs
    ix = find(tank.snips_ch == ch & tank.snips_sc == sc);
    rix = ix(randperm(min(SHOW, length(ix))));
    plot(1:size(tank.snips_v, 2), tank.snips_v(rix, :)', ...
         [pcolors(ns) '-']);
    ns = ns + 1;
    hold on;
    title(sprintf('ch=%d sc=%d', ch, sc));
  end
  hold off;
end


