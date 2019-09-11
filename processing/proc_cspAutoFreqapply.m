function [dat, varargout]= proc_cspAutoFreqapply(dat, freqs, varargin)
%PROC_CSPAUTO - Common Spatial Pattern Analysis with Auto Filter Selection
%
%Synopsis:
% [DAT, CSP_W, CSP_A, CSP_EIG]= proc_cspAuto(DAT, <OPT>);
% [DAT, CSP_W, CSP_A, CSP_EIG]= proc_cspAuto(DAT, NPATS);
%
%Arguments:
% DAT    - data structure of epoched data
% NPATS  - number of patterns per class to be calculated, 
%          default nChans/nClasses.
% OPT - struct or property/value list of optional properties:
%  .patterns - 'all' or matrix of given filters or number of filters. 
%      Selection depends on opt.selectPolicy.
%  .selectPolicy - determines how the patterns are selected. Valid choices are 
%      'all', 'equalperclass', 
%      'maxvalues' (attention: possibly not every class is represented!)
%      'maxvalueswithcut' (like maxvalues, but components with less than
%      half the maximal score are omitted),
%      'directorscut' (default, heuristic, similar to maxvalueswithcut, 
%      but every class will be represented), 
%      'matchfilters' (in that case opt.patterns must be a matrix),
%      'matchpatterns' (not implemented yet)
%  .covPolicy - 'normal', 'average' (default) or manually defined cov matrices
%      of size [nchans, nchans, 2] . 'average' calculates the average of the 
%      single-trial covariance matrices.
%  .score - 'eigenvalues', 'medianvar' (default) or 'auc' to
%       determine the components by the respective score
%
%Returns:
% DAT    - updated data structure
% CSP_W  - CSP projection matrix (spatial filters, in the columns)
% CSP_A  - estimated mixing matrix (activation patterns, in the columns)
% CSP_EIG- eigenvalue score of CSP projections 
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

% Author(s): Benjamin Blankertz
props= { 'patterns'     3           'INT|CHAR'
         'score'        'medianvar' '!CHAR(eigenvalues medianvar auc)'
         'covPolicy'    'average'   'CHAR|DOUBLE[- - 2]'
         'scaling'      'none'      'CHAR'
         'normalize'    0           'BOOL'
         'selectPolicy' 'directorscut'  'CHAR'
         'weight'       ones(1,size(dat.y,2))   'DOUBLE'
         'weightExp'    1           'BOOL'};

if nargin==0,
  dat = props; return
end

dat = misc_history(dat);
misc_checkType(dat, 'STRUCT(x clab y)'); 
if length(varargin)==1 && isnumeric(varargin{1}),
  opt= struct('patterns', varargin{1});
else
  opt= opt_proplistToStruct(varargin{:});
end
[opt, isdefault]= opt_setDefaults(opt, props);
opt_checkProplist(opt, props);

[T, nChans, nEpochs]= size(dat.x);

if size(dat.y,1)~=2,
  error('this function works only for 2-class data.');
end

[filt_b, filt_a] = butters(4,freqs/dat.fs*2);

dat=proc_filterbank(dat,filt_b,filt_a);

[dat, W, A, score]=proc_multiBandSpatialFilter(dat,...
    {@proc_cspAuto,opt});
varargout={dat,W,A,score,freqs};   
