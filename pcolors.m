function c = pcolors(n)
colors = 'kyrckyrc';
%c = colors(mod(n-1, length(colors)) + 1);
c = colors(n+1);
