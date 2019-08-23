function [dat2, varargout]= proc_csssp_onlineNoiseupdate(dat, W,score,C, varargin)
%PROC_CSSDP - Common Spatio-Frequency Decomposition Pattern (CSP) Analysis
%
%Synopsis:
% [DAT, CSSDP_W, CSSDP_A, SCORE]= proc_cssdp_prefilt(DAT, <OPT>);
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
        'lambda'      0.001                              'DOUBLE'
        'alpha'       0.1                             'DOUBLE'
        'Verbose'     1                                 'INT'
       };

if nargin==0,
  dat= props; return
end

misc_checkType(dat, 'STRUCT(x y)'); 
misc_checkType(dat.y, 'DOUBLE[2 -]', 'dat.y');    % two classes only
opt= opt_proplistToStruct(varargin{:});
opt= opt_setDefaults(opt, props);
opt_checkProplist(opt, props);
dat= misc_history(dat);

% Calculate classwise covariance matrices
[covFcn, covPar]= misc_getFuncParam(opt.CovFcn);

epo_noise=proc_selectChannels(dat,'*noise*');
dat=proc_selectChannels(dat,'not','*noise*');
nChans= size(W, 2);
nEpo= size(epo_noise.x, 3);
dat2= proc_linearDerivation(dat, W, 'prependix','csssp');
D=diag(score);

for ii=1:nEpo
    X_n=epo_noise.x(:,:,ii);
    C_temp= covFcn(dat.x(:,:,ii), covPar{:});
    C=(1-opt.lambda)*C+opt.lambda*C_temp;
    [~,D_k]=eig(C_temp-C);
    C_k=abs((C_temp-C)*sign(D_k));
    C_n= (1-opt.lambda)*eye(nChans)+opt.lambda*W'*(covFcn(X_n, covPar{:})+opt.alpha*C_k)*W;
    % ORIGINAL CODE FOR COMPUTING CSSDP IN CHANNEL SPACE
    [V, D]= eig( D, C_n );
    %resort (eig mixes them up) THIS IS RATHER A HOTFIX
    [~,imax]=max(abs(V));
    [~,inds]=sort(imax);
    V=V(:,inds).*repmat(sign(diag(V(:,inds))),1,nChans)';
    D=diag(diag(D(inds,inds)));
    W=W*V;
    dat2.x(:,:,ii)=dat.x(:,:,ii)*W;
end

% Calculate score for each CSP channel
% [scoreFcn, scorePar]= misc_getFuncParam(opt.ScoreFcn);
% score= scoreFcn(dat2, W, D, scorePar{:});
    

% Save old channel labels
if isfield(dat2, 'clab'),
  dat2.origClab= dat2.clab;
end

% Apply CSP filters to time series
% dat= proc_linearDerivation(dat, W, 'prependix','cssdp');

% Determine patterns according to [Haufe et al, Neuroimage, 87:96-110, 2014]
% http://dx.doi.org/10.1016/j.neuroimage.2013.10.067
% C_avg = mean(C,3);
% A= C_avg * W / (W'*C_avg*W);

varargout= {W , C};
