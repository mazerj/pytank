function blocks = list_tdtblocks(pf)
%function blocks = list_tdtblocks(pf)
%
% get list of tdt blocks (including name mangling from local
% storage target to raid location) containing analog data
% traces for specified pypefile.
%
% INPUTS
%  pf - pypefile struct from p2mLoad
%
% OUTPUTS
%  blocks - cell array list of block names
%

blocks = {};
for n = 1:length(pf.rec)
  t = pf.rec(n).params.tdt_tank;
  t = strrep(t, '\', '/');
  t = strrep(t, 'C:', '/auto/data/critters');
  t = [t '/' pf.rec(n).params.tdt_block];
  blocks{n} = t;
end
blocks = unique(blocks);
