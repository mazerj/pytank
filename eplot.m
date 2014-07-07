function eplot(x, y, ye)
%function eplot(x, y, ye)
%
% Uses eshade to generate nice looking simple x,y,error plots in
% one function call.
%

if ~ishold
  cla;                                  % otherwise old eshade's persist..
end
h = ishold;
hold on;
eshade(x, y, ye);
plot(x, y, 'r-');
if ~h, hold off, end;

