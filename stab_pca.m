function stab_pca(t, ch, sc, flat)
%function stab_pca(t, ch, sc, flat)
%
%  Assess sort stability over time. This computes the principle
%  components of the entire set of TDT snips for the specified channel
%  and sort and then plots the temporal evolution of the
%  loadings across the first 3 PCs. If flat is specified, this is
%  done as a matrix of scatter plots, otherwise a 3d cloud.
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
[~, scores, ~] = pca(snips);

N=50;
bn = round(size(scores,1)/N);
cm = icolormap([1 0 0; 0 1 0; 0 0 1], N/2);

if flat
  NPC = 3;
  np = 1;
  for pc1 = 1:NPC
    for pc2 = 1:NPC
      if pc2 >= pc1
        pcs = [pc1 pc2];
        subplot(NPC,NPC,np);
        mx = [];
        my = [];
        for n = 1:N
          ix = ((n-1)*bn)+(1:bn);
          ix = ix(ix <= size(scores,1));
          mx = [mx mean(scores(ix, pcs(1)))];
          my = [my mean(scores(ix, pcs(2)))];
          while length(ix) > MAXPTS
            ix = ix(1:2:end);
          end
          set(plot(scores(ix,pcs(1)), scores(ix,pcs(2)), '.'), ...
              'Color', cm(n,:));
          hold on;
        end
        
        set(plot(mx(1), my(1), 'go'), ...
            'markerfacecolor', 'g', 'markeredgecolor', 'k');
        set(plot(mx(end), my(end), 'ro'), ...
            'markerfacecolor', 'r', 'markeredgecolor', 'k');
        arrow([mx(1) my(1)], [mx(end) my(end)]);
        
        axis square;
        hold off;
        xlabel(sprintf('pc%d', pcs(1)));
        ylabel(sprintf('pc%d', pcs(2)));
      end
      np = np + 1;
    end
  end
else
  mx = [];
  my = [];
  mz = [];
  pcs = 1:3;
  for n = 1:N
    ix = ((n-1)*bn)+(1:bn);
    ix = ix(ix <= size(scores,1));
    mx = [mx mean(scores(ix, pcs(1)))];
    my = [my mean(scores(ix, pcs(2)))];
    mz = [mz mean(scores(ix, pcs(3)))];
    while length(ix) > MAXPTS
      ix = ix(1:2:end);
    end
    set(plot3(scores(ix,pcs(1)), scores(ix,pcs(2)), ...
              scores(ix,pcs(3)), '.'), 'Color', cm(n,:));
    hold on;
  end
  
  set(plot3(mx(1), my(1), mz(1), 'go'), ...
      'markerfacecolor', 'g', 'markeredgecolor', 'k');
  set(plot3(mx(end), my(end), mz(end), 'ro'), ...
      'markerfacecolor', 'r', 'markeredgecolor', 'k');
  arrow([mx(1) my(1) mz(1)], [mx(end) my(end) mz(end)]);
  
  axis square;
  hold off;
  xlabel(sprintf('pc%d', pcs(1)));
  ylabel(sprintf('pc%d', pcs(2)));
  zlabel(sprintf('pc%d', pcs(3)));
  grid on;
end  

title(sprintf('ch=%d sc=%d (time: R->G->B)', ch, sc));
