function c = pcolors(n)
colors = 'rgbymck';
c = colors(mod(n-1, length(colors)) + 1);
