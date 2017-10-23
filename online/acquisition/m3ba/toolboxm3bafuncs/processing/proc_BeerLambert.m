function [dat,state] = proc_BeerLambert(dat, state, varargin)
% PROC_BEERLAMBERT - apply the BeerLambert function, which is useful to analyze
% NIRS data:  relative concentration is calculated as a function of total photon path length
%
% Synopsis:
%   DAT = nirs_LB(DAT, 'Property1',Value1, ...)
%
% Arguments:
%   DAT: data struct with NIRS cnt data split in the x field
%
%   OPT - struct or property/value list of optional properties:
%   'Citation'   - if set, the epsilon (extinction values) is taken from
%                  the specified citation number (see GetExtinctions.m for
%                  details). (default 1 = Gratzer et al)
%   'Epsilon'    - sets extinction coefficients manually (wl1 wl2 vs deoxy
%                  and oxy). In this case, citation is ignored. State in
%                  millimol/liter(?).
%   'Opdist'     - optode (source-detector) distance in cm (default 2.5)
%   'Ival'       - either 'all' (default) or a vector [start end] in samples
%                  specifying the baseline for the LB transformation
%   'DPF'        - differential pathlength factor: probabilistic average
%                  distance travelled by photons, default [5.98 7.15]
%   'Baseline'   - individual baseline of dimension: nchannels x 1
%                  if set, Ival will be ignored
%
% Returns:
%   DAT - updated data struct with fields specifying absorption
%         values in mmol/l.
%         dat.x = timepoints x oxy/deoxy x channel
%
% See also: nirs_* nirsfile_* GetExtinctions
%
% Note: Based on the nirX Nilab toolbox functions u_LBG and u_popLBG.

% matthias.treder@tu-berlin.de 2011, mail@aewald.net 2013
% AE: changed channels label from wl1&wl2 to oxa and deoxy.
% Markus Wenzel 2013 (adapted it to the new toolbox)
% Jan Mehnert February 2014 (ready for public BBCI toolbox) (jan@mehnert.org)
% Stephanie Brandl (stephanie.brandl@tu-berlin.de):
% bug fixing, online functionality and individual baseline
% results now coincide with homer2 toolbox


props={'Citation'   1             'INT'
    'Opdist'    2.5           'DOUBLE'
    'Ival'      'all'         'CHAR|DOUBLE'
    'DPF'       [5.98 7.15]   'DOUBLE'
    'Epsilon'   []            'DOUBLE'
    'Verbose'   0             'BOOL'
    'Baseline'  []            'DOUBLE'
    'idx'       []            'DOUBLE[2-]'
    'clab'      []             'CELL'};

if nargin==0,
    dat= props; return
end

opt= opt_proplistToStruct(varargin{:});
[opt, ~]= opt_setDefaults(opt, props);

opt_checkProplist(opt, props);
misc_checkType(dat, 'STRUCT');

if strcmp(opt.Ival,'all')
    opt.Ival = [1 size(dat.x,1)];
end

if isempty(opt.idx)
    [wl1,idx1] = proc_selectChannels(dat,'*low*');
    idx1 = find(sum(idx1,2)==1);
    [wl2, idx2] = proc_selectChannels(dat,'*high*');
    idx2 = find(sum(idx2,2)==1);
    out.idx = cat(1,idx1,idx2);
else
    idx1 = idx(1,:);
    idx2 = idx(2,:);
end

[s1,s2] =size(wl1.x);   %number of timepoints %number of channels

wl1.x = wl1.x+eps;
wl2.x = wl2.x+eps;

if ~isempty(opt.Baseline)
    if size(opt.Baseline,2)==1
        opt.Baseline=opt.Baseline';
    end
    opt.Baseline = opt.Baseline+eps;
end

%% Get epsilon
if isempty(opt.Epsilon)
    if ~isfield(dat,'wavelengths')
        error('Wavelengths should be given in the .wavelengths field.')
    end
    [ext,nfo] = procutil_getExtinctions(dat.wavelengths,opt.Citation);
    if opt.Verbose, fprintf('Citation: %s\n',nfo), end
    epsilon = ext(:,1:2);
    
    % Divide by 1000 to obtain the required unit
    epsilon = epsilon/1000;
    out.Epsilon = epsilon;
    
else
    epsilon = opt.Epsilon;
end

%% Arrange epsilon so that higher wavelength is on top
[~,idx] = max(dat.wavelengths);
if idx==1 % higher wavelength is on bottom
    epsilon = flipud(epsilon);
    if opt.Verbose, fprintf('Epsilon matrix was rearranged so that lower WL is on top\n'), end
end

%% Apply Lambert-Beer law
if isempty(opt.Baseline)
    
    Att_lowWL= -log(abs( wl1.x)./ ...
        ( repmat(mean(wl1.x(opt.Ival(1):opt.Ival(2),:),1), [s1,1]))   );
    
    Att_highWL= -log(abs( wl2.x) ./ ...
        ( repmat(mean(wl2.x(opt.Ival(1):opt.Ival(2),:),1), [s1,1]))   );
    
else
    
    Att_lowWL= -log(abs( wl1.x)./ ...
        ( repmat(opt.Baseline(idx1), [s1,1]))   );
    
    Att_highWL= -log(abs( wl2.x )./ ...
        ( repmat(opt.Baseline(idx2), [s1,1]))   );
    
    out.Baseline = opt.Baseline;
    
end
%
A=[];
A(:,1)=reshape(Att_lowWL,s1*s2,1);
A(:,2)=reshape(Att_highWL,s1*s2,1);


%----------------------------------
%       3.cc
%----------------------------------
% e=...looks like this
%               oxy-Hb         deoxy-Hb
% higherWL: 830 | e: 0.974       0.693
% lowerWL : 690 | e: 0.35         2.1

%not sure why this was in here (maybe because of line 42 in hmr2ODConc)
% e= epsilon/10;

e = epsilon;

c= (e\(A./(ones(s2*s1,1)*opt.Opdist*opt.DPF))');

dat.x = [];
% dat.x(:,1,:) = reshape(c(1,:),s1,s2); %in mmol/l
% dat.x(:,2,:) = reshape(c(2,:),s1,s2); %in mmol/l
dat.x= reshape(c(1,:),s1,s2); %in mmol/l
dat.x= [dat.x reshape(c(2,:),s1,s2)];
%% Change Channel labels wl1 & wl2 to 'oxy' and 'deoxy'

if isempty(opt.clab)
    dat.clab = strrep(wl1.clab, 'low', '');
    dat.clab = strrep(dat.clab, 'WL', '');
    out.clab = dat.clab;
else
    dat.clab = opt.clab;
end

dat.signal = 'NIRS (oxy, deoxy)';

dat.yUnit = 'mmol/l';