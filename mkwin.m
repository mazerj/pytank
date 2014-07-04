function mkwin(new, x, y, w, h)
%function mkwin(new, x, y, w, h)
%
% create or move window at location specified (xy = upper left corner),
% where values are fractions [0-1] of the total screen size.
%

s = get(0, 'ScreenSize')
sw = s(3);
sh = s(4);

if new
  f = figure(new);
else
  f = gcf;
end

w = round(sw * w);
h = round(sh * h);
x = round(x * sw);
y = round(y * sh);
y = sh - y;

set(f, 'Position', [x y w h]);


