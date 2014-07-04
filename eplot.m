function eplot(x, y, ye)

if ~ishold
  cla;                                  % otherwise old eshade's persist..
end
h = ishold;
hold on;
eshade(x, y, ye);
plot(x, y, 'r-');
if ~h, hold off, end;

