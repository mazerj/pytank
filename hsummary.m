function hsummary(pattern, hide)

if ~exist('hide', 'var')
  hide = 0;
end

fig = figure;
if hide
  set(fig, 'visible', 'off');
end
pos = get(fig, 'position');
pos(3) = 700;
pos(4) = 900;
set(gcf, 'position', pos);
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
    if isempty(t_)
      if seg == 1
        error('no hdf5: %s; run exper2hdf5', basename(pf.src))
      else
        break
      end
    end
    v_ = v_ * 1e6;
    t_ = t_ / 60;
    ve = nanstd(v_);
    maxrms = max([maxrms ve]);
    
    times = t_(1:dt:end);
    volts = v_(1:dt:end);
    if seg == 1
      t0 = times(1);
    end
    
    subplot(1+length(flist), 1, 1);
    plot(globaloffset+(times-t0), volts, c(rem(n,2)+1));
    hold on;
    ylabel('uV');
    xlabel('time (min)');
    title(sprintf('%s [lines indicate ~ +-5sig]', pattern));
    
    subplot(1+length(flist), 1, n+1);
    plot(times-t0, volts);
    hold on;
    ylabel(strsplit(basename(pf.src), '.'));
    xlabel('time (min)');
    drawnow;
    
    seg = seg + 1;
  end
  hline(5*ve, 'color', 'r', 'linestyle', '-');
  hline(-5*ve, 'color', 'r', 'linestyle', '-');
  globaloffset = globaloffset + (times(end)-t0) + 0;
end

for n = 0:length(flist)
  subplot(1+length(flist), 1, n+1);
  yrange(-20*maxrms, 20*maxrms);
end

if hide
  set(fig, 'visible', 'on');
end
