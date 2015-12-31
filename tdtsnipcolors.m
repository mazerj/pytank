function cname = tdtsnipcolors(no)
%function cname = tdtsnipcolors(no)
%
% convert snip/sort number to a plotable color
%
%INPUT
%  no - snip/box/sort number
%
%OUTPUT
%  cname - color name as string
%

colors = 'kyrgmbyrgmbyrgmb';
cname = colors(no+1);
