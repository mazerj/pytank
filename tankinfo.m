function tankinfo(pattern)
%function tankinfo(pattern)
%
%  show list of tanks associcated with making pypefiles
%

files = dbfind(pattern, 'list', 'all');

for k = 1:length(files)
  pf = p2mLoad2(files{k});
  lasttank = '';
  lastblack = '';
  fprintf('%s\n', basename(pf.src))
  for n = 1:length(pf.rec)
    if strcmp(pf.rec(n).params.tdt_tank, lasttank) == 0 || ...
          strcmp(pf.rec(n).params.tdt_block, lastblock) == 0
      t = sprintf('%s\\%s', ...
                  pf.rec(n).params.tdt_tank, ...
                  pf.rec(n).params.tdt_block);
      t = strrep(t, 'C:\', '/auto/data/critters/');
      t = strrep(t, '\', '/');
      fprintf(' --> %s\n', t)
      lasttank = pf.rec(n).params.tdt_tank;
      lastblock = pf.rec(n).params.tdt_block;
    end
  end
end
