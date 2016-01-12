function stab_met(t, ch, sc, flat)
%function stab_met(t, ch, sc, flat)
%
%  Assess sort stability over time. This uses heuristics
%  instead of PCs (see stab_pca): peak and trough voltages
%  and the spike width (time between peak and trough). Metrics
%  and computed coarsely for speed -- no effort to be exact.
%
%  Note: this looks only at the on-line TDT sorts!
%
%INPUTS
%  t - tank data from pyt_load()
%  ch - channel number (starting with 1)
%  sc - sort code (starting with 0, for unsorted)
%  flat - 2d vs 3d plots
%
%OUPUTS
%  none
%


if ~exist('flat', 'var')
  flat = 1;
end

MAXPTS = 1000;

ix = find(t.snips_ch == ch & t.snips_sc == sc);
snips = t.snips_v(ix, :);
fs = 1.0 / mean(diff(t.snips_t(1,:)));


vmax = max(snips, [], 2);
vmin = min(snips, [], 2);
wid = zeros(size(vmax));
for n = 1:length(wid)
  wid(n) = 1e6 * abs(mean(find(snips(n,:) == vmax(n))) - ...
                     mean(find(snips(n,:) == vmin(n)))) / fs;
end
scores = [1e6*vmax 1e6*vmin wid];

names = {'vmax (uv)', 'vmin (uv)', 'width (us)'};

N=50;
bn = round(size(scores,1)/N);
cm = icolormap([1 0 0; 0 1 0; 0 0 1], N/2);

if flat
  NMET = size(scores, 2);
  np = 1;
  for met1 = 1:NMET
    for met2 = 1:NMET
      if met2 >= met1
        mets = [met1 met2];
        subplot(NMET,NMET,np);
        mx = [];
        my = [];
        for n = 1:N
          ix = ((n-1)*bn)+(1:bn);
          ix = ix(ix <= size(scores,1));
          mx = [mx mean(scores(ix, mets(1)))];
          my = [my mean(scores(ix, mets(2)))];
          while length(ix) > MAXPTS
            ix = ix(1:2:end);
          end
          set(plot(scores(ix,mets(1)), scores(ix,mets(2)), '.'), ...
              'Color', cm(n,:));
          hold on;
        end
        
        set(plot(mx(1), my(1), 'go'), ...
            'markerfacecolor', 'g', 'markeredgecolor', 'k');
        set(plot(mx(end), my(end), 'ro'), ...
            'markerfacecolor', 'r', 'markeredgecolor', 'k');
        arrow([mx(1) my(1)], [mx(end) my(end)]);
        hold off;
        grid on;
        xlabel(names{met1});
        ylabel(names{met2});
        if np == 1
          title(sprintf('ch=%d sc=%d (time: R->G->B)', ch, sc));
        end
      end
      np = np + 1;
    end
  end
else
  mx = [];
  my = [];
  mz = [];
  mets = 1:3;
  for n = 1:N
    ix = ((n-1)*bn)+(1:bn);
    ix = ix(ix <= size(scores,1));
    mx = [mx mean(scores(ix, mets(1)))];
    my = [my mean(scores(ix, mets(2)))];
    mz = [mz mean(scores(ix, mets(3)))];
    while length(ix) > MAXPTS
      ix = ix(1:2:end);
    end
    set(plot3(scores(ix,mets(1)), scores(ix,mets(2)), ...
              scores(ix,mets(3)), '.'), 'Color', cm(n,:));
    hold on;
  end
  set(plot3(mx(1), my(1), mz(1), 'go'), ...
      'markerfacecolor', 'g', 'markeredgecolor', 'k');
  set(plot3(mx(end), my(end), mz(end), 'ro'), ...
      'markerfacecolor', 'r', 'markeredgecolor', 'k');
  arrow([mx(1) my(1) mz(1)], [mx(end) my(end) mz(end)]);

  axis square;
  hold off;
  xlabel(names{1});
  ylabel(names{2});
  zlabel(names{3});
  grid on;
  title(sprintf('ch=%d sc=%d (time: R->G->B)', ch, sc));
end  

