function [dat, varargout]= proc_cssdp(cnt,mrk,freqs, varargin)
%PROC_CSSDP - Common Spatio-Frequency Decomposition Pattern (CSP) Analysis
%
%Synopsis:
% [DAT, CSP_W, CSP_A, SCORE]= proc_csp(CNT,MRK,FREQS, <OPT>);
%
%Arguments:
% DAT - data structure of epoched data
% OPT - struct or property/value list of optional properties:
%  'CovFcn' - Handle of a function to estimate a covariance matrix.
%             Can also be a CELL {@FCN, PARAM1, PARAM2, ...} that
%             provides further arguments to that function.
%             Default: @cov
%             Other options, e.g., @procutil_covShrinkage
%  'ScoreFcn' - Handle of a function that calculates a score of
%               the extracted components.
%               Can also be a CELL {@FCN, PARAM1, PARAM2, ...}.
%               Default: @procutil_score_eigenvalues
%               Other options, e.g., @procutil_score_ratioOfMedians
%  'SelectFcn' - Handle of a function that selects a subset of the
%                extracted components.
%                Can also be a CELL {@FCN, PARAM1, PARAM2, ...}.
%                Default: {@procutil_selectMinMax, 3}
%                Other options, e.g., {@procutil_selectMaxAbs, 6}
%  'Verbose' - Print warnings and other output if larger than 0. Default 1
% 
%Returns:
% DAT   - updated data structure
% CSP_W - CSP projection matrix (spatial filters, in the columns)
% CSP_A - estimated mixing matrix (activation patterns, in the columns)
% SCORE - score of each CSP component
%
%Description:
% calculate common spatial patterns (CSP).
% please note that this preprocessing uses label information. so for
% estimating a generalization error based on this processing, you must
% not apply csp to your whole data set and then do cross-validation
% (csp would have used labels of samples in future test sets).
% you should use the .proc (train/apply) feature in xvalidation, see
% demos/demo_validation_csp
%
%See also demos/demo_validate_csp

props= {'CovFcn'      {@cov}                            '!FUNC|CELL'
        'ScoreFcn'    {@score_eigenvalues}              '!FUNC|CELL'
        'SelectFcn'   {@cspselect_equalPerClass, 3}     '!FUNC|CELL'
        'Verbose'     1                                 'INT'
        'filterOrder'   3                               'INT'
        'ival'  []                               'DOUBLE[- -2]'
       };

if nargin==0,
  dat= props; return
end

misc_checkType(cnt, 'STRUCT(x)'); 
misc_checkType(mrk.y, 'DOUBLE[2 -]', 'mrk.y');    % two classes only
opt= opt_proplistToStruct(varargin{:});
opt= opt_setDefaults(opt, props);
opt_checkProplist(opt, props);
cnt= misc_history(cnt);

% make sure FREQS has the correct dimensions
if not( size(freqs,1)==3 && size(freqs,2)==2 )
  error('freq must be a 3 by 2 matrix, i.e. three bands must be specified!');
end

% check the given frequency bands
signal_band = freqs(1,:); % signal bandpass band
noise_bp_band = freqs(2,:); % noise bandpass band
noise_bs_band = freqs(3,:); % noise bandstop band
if not( noise_bs_band(1) < signal_band(1) && ...
        noise_bp_band(1) < noise_bs_band(1) && ...
        signal_band(2) < noise_bs_band(2) && ...
        noise_bs_band(2) < noise_bp_band(2) )
  error('Wrongly specified frequency bands!\nThe first band (signal band-pass) must be within the third band (noise band-stop) and the third within the second (noise band-pass)!');
end


nChans= size(cnt.x, 2);

% Creating filters
[b,a]=butter(opt.filterOrder, signal_band/(cnt.fs/2));
[b_f,a_f]=butter(opt.filterOrder, noise_bp_band/(cnt.fs/2));
[b_s,a_s]=butter(opt.filterOrder, noise_bs_band/(cnt.fs/2),'stop');


% Covariance matrix for the center frequencies (signal)
X_s = proc_filtfilt(cnt,b,a);
epo = proc_segmentation(X_s, mrk, opt.ival);
% Calculate classwise covariance matrices
[covFcn, covPar]= misc_getFuncParam(opt.CovFcn);
C= zeros(nChans, nChans, 2);
for k= 1:2,
  idx= find(epo.y(k,:));
  X= permute(epo.x(:,:,idx), [1 3 2]);
  X= reshape(X, [], nChans);
  C(:,:,k)= covFcn(X, covPar{:});
end


% Covariance matrix for the flanking frequencies (noise)
X_tmp = proc_filtfilt(cnt,b_f,a_f);
X_tmp = proc_filtfilt(X_tmp,b_s,a_s);
X_tmp = proc_segmentation(X_tmp, mrk, opt.ival);
X_tmp= permute(X_tmp.x, [1 3 2]);
X_tmp= reshape(X_tmp,[], nChans);
C_n = cov(X_tmp,1);
clear X_tmp



% get the whitening matrix
M = procutil_whiteningMatrix([], 'C', mean(C,3));
if (opt.Verbose > 0) && (size(M,2) < nChans)
    warning('Due to dimensionality reduction a maximum of only %d CSP components can be computed', size(M,2))
end

% ORIGINAL CODE FOR COMPUTING CSP IN CHANNEL SPACE
% % Do actual CSP computation as generalized eigenvalue decomposition
[W, D]= eig( C(:,:,1)-C(:,:,2), C_n );

% Calculate score for each CSP channel
[scoreFcn, scorePar]= misc_getFuncParam(opt.ScoreFcn);
score= scoreFcn(cnt, W, D, scorePar{:});

% Select desired CSP filters
[selectFcn, selectPar]= misc_getFuncParam(opt.SelectFcn);
idx= selectFcn(score, W, selectPar{:});
W= W(:,idx);
score= score(idx);

% Save old channel labels
if isfield(cnt, 'clab'),
  dat.origClab= cnt.clab;
end

% Apply CSP filters to time series
cnt= proc_linearDerivation(cnt, W, 'prependix','csp');

% Determine patterns according to [Haufe et al, Neuroimage, 87:96-110, 2014]
% http://dx.doi.org/10.1016/j.neuroimage.2013.10.067
C_avg = mean(C,3);
A= C_avg * W / (W'*C_avg*W);

varargout= {W, A, score};
