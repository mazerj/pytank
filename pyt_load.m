function tank = pyt_load(exper, what)
%function pyt_load(exper, what)
%
% load all data associated with experiment (multiple pypefiles)
%
% what is string with:
%   's' for snips
%   'h' for high-pass continuous spike signal
%   'l' for low-pass continuous lfp signal

try
  DUMP = HDF5DUMP;
catch
  DUMP = '/auto/th5';
  warning('loading from /auto/th5');
end

if ~exist('what', 'var') || isempty('what')
  what = 'shl';                         % load everything
end

pflist = dbfind(exper, 'list', 'all')';

tank.exper = exper;
tank.srcs = pflist;                % p2m/pf #'d source files
tank.snips_v = [];                      % snip voltage trace
tank.snips_t = [];                      % snip time base (file relative)
tank.snips_ch = [];                     % snip channel
tank.snips_sc = [];                     % snip sort code (0=unsorted)
tank.snips_srcn = [];                   % snip source file #
tank.highpass = [];                     % continuous HP voltage trace
tank.highpass_t = [];                   % continuous HP time trace
tank.highpass_srcn = [];                % HP source file #
tank.lowpass = [];                      % continuous LP voltage trace
tank.lowpass_t = [];                    % continuous LP time trace
tank.lowpass_srcn = [];                 % LP source file #
tank.t0 = 0;


for pfn = 1:length(tank.srcs)
  fprintf('%s\n', tank.srcs{pfn});
  pf = p2mLoad2(tank.srcs{pfn});

  % get list of all tank blocks referenced in pypefile
  blocks = list_tdtblocks(pf);

  hfiles = {};
  for bn = 1:length(blocks)
    blk = strsplit(blocks{bn}, '/');
    % for each block, find all segment files
    
    segfiles = jls(sprintf('%s/%s-%s_???.th5', ...
                           DUMP, blk{end-1}, blk{end}));
    for n = 1:length(segfiles)
      hfiles{length(hfiles)+1} = segfiles{n};
    end
  end
  
  if 0 & isnan(tank.t0)
    % time for very first data sample in entire dataset
    tank.t0 = h5readatt(hfiles{1}, '/continuous/spk', 'tstart');
  end
  
  if any('s' == what)
    for hn = 1:length(hfiles)
      hf = hfiles{n};
      tank.snips_v = [tank.snips_v; h5read(hf, '/snip/v')'];
      tank.snips_t = [tank.snips_t; h5read(hf, '/snip/t')'-tank.t0];
      tank.snips_ch = [tank.snips_ch; h5read(hf, '/snip/ch')];
      sc = h5read(hf, '/snip/sc');
      tank.snips_sc = [tank.snips_sc; sc];
      tank.snips_srcn = [tank.snips_srcn; uint16(pfn + zeros(size(sc)))];
    end
  end

  if any('h' == what)
    for hn = 1:length(hfiles)
      hf = hfiles{n};
      v = h5read(hf, '/continuous/spk')';
      t = linspace(h5readatt(hf, '/continuous/spk', 'tstart'), ...
                   h5readatt(hf, '/continuous/spk', 'tend'), ...
                   length(v))';
      tank.highpass = [tank.highpass; v; NaN];
      tank.highpass_t = [tank.highpass_t; t-tank.t0; NaN];
      tank.highpass_srcn = [tank.highpass_srcn; ...
                          uint16(pfn + zeros(size(v))); 0];
    end
  end

  if any('l' == what)
    for hn = 1:length(hfiles)
      hf = hfiles{n};
      v = h5read(hf, '/continuous/lfp')';
      t = linspace(h5readatt(hf, '/continuous/lfp', 'tstart'), ...
                   h5readatt(hf, '/continuous/lfp', 'tend'), ...
                   length(v))';
      tank.lowpass = [tank.lowpass; v; NaN];
      tank.lowpass_t = [tank.lowpass_t; t-tank.t0; NaN];
      tank.lowpass_srcn = [tank.lowpass_srcn; ...
                          uint16(pfn + zeros(size(v))); NaN];
    end
  end
end

