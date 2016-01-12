function ls = pyt_eplot(x, y, ye)
%function ls = pyt_eplot(x, y, ye)
%
% Uses eshade to generate nice looking simple x,y,error plots in
% one function call.
%

if ~ishold
  cla;                                  % otherwise old eshade's persist..
end
h = ishold;
hold on;
if length(ye) > 1
  ls = [eshade(x, y, ye); plot(x, y, 'r-')];
else
  ls = [plot(x, y, 'r-')];
end
if ~h, hold off, end;

