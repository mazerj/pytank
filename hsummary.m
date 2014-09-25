function hsummary(pattern)

clf;
dt = 40;

% train with data from first segment from first datafile
flist = dbfind(pattern, 'list', 'all');
maxrms = -Inf;
globaloffset = 0;
c = 'kr';
for n = 1:length(flist)
  pf = p2mLoad2(flist{n});
  seg = 1;
  times = []; volts = [];
  while 1
    [t_, v_, ~] = hload(pf, seg);
    if isempty(t_), break; end
    v_ = v_ * 1e6;
    t_ = t_ / 60;
    maxrms = max([maxrms std(v_)]);
    
    times = t_(1:dt:end);
    volts = v_(1:dt:end);
    if seg == 1
      t0 = times(1);
    end
    seg = seg + 1;
    
    subplot(1+length(flist), 1, n+1);
    plot(times-t0, volts);
    hold on;
    ylabel(basename(pf.src));
    set(get(gca, 'ylabel'), 'Rotation', 70);
    xlabel('time (min)');

    subplot(1+length(flist), 1, 1);
    plot(globaloffset+(times-t0), volts, c(rem(n,2)+1));
    hold on;
    ylabel('voltage (uv)');
    set(get(gca, 'ylabel'), 'Rotation', 70);
    xlabel('time (min)');
    title(pattern);
    
    drawnow;
  end
  globaloffset = globaloffset + (times(end)-t0) + .2;
end

for n = 0:length(flist)
  subplot(1+length(flist), 1, n+1);
  yrange(-20*maxrms, 20*maxrms);
end
