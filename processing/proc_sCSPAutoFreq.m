function [fv, varargout]= proc_sCSPAutoFreq(dat, freqs, varargin)
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
props= {'patterns'     3           'INT|CHAR'
        'ival'          []          'DOUBLE|CHAR'
        'maxIval'       [250 5000]  'DOUBLE[- 2]'
        'startIval'     [750 3500]  'DOUBLE[- 2]'
        'score'        'medianvar' '!CHAR(eigenvalues medianvar auc)'
        'covPolicy'    'average'   'CHAR|DOUBLE[- - 2]'
        'scaling'      'none'      'CHAR'
        'normalize'    0           'BOOL'
        'selectPolicy' 'directorscut'  'CHAR'
        'weight'       []   'DOUBLE'
        'weightExp'    1           'BOOL'
        'alpha'      1                              'DOUBLE'
        'chunksize'      10                              'DOUBLE'
        };

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
if isempty(opt.weight)
    opt.weight=ones(1,size(dat.y,2)) ;
end

[T, nChans, nEpochs]= size(dat.x);

if size(dat.y,1)~=2,
    error('this function works only for 2-class data.');
end
clab_load=dat.clab;

dat_lap=proc_laplacian(dat);

%% first apply a broad-band filt before using the heuristic
[filt_b, filt_a]= butter(5, [min(freqs(:,1)) max(freqs(:,2))]/dat.fs*2);
dat_flt= proc_filt(dat_lap, filt_b, filt_a);

%% select best time interval on broad band filtered data
if ischar(opt.ival)&&strcmp(opt.ival,'auto')
    ival= select_timeivalEpo(dat_flt,'maxIval',opt.maxIval,...
    'startIval',opt.startIval);
    % if diff(ival)<1000
    %     ival= [mean(ival)-500 mean(ival)+500];
    % end
elseif isempty(opt.ival)
    ival=[dat.t(1) dat.t(end)];
else
    ival=opt.ival;
end
dat_lap=proc_selectIval(dat_lap,ival);
%%
bands=[];
for iFreq=1:size(freqs,1)
%band= select_bandnarrow_epo(dat_lap, 'band',freqs(iFreq,:),...
%    'bandTopscore',freqs(iFreq,:),'areas',{});
band= select_bandbroad_epo(dat_lap, 'band',freqs(iFreq,:));
if band(1)==band(2)
    band=freqs(iFreq,:);
end
bands=[bands;band];
%     fprintf('freqs %f %f selected\n',band(1),band(2))
end

[filt_b, filt_a] = butters(4,bands/dat.fs*2);
dat=proc_filterbank(dat,filt_b,filt_a);

if ischar(opt.ival)&&strcmp(opt.ival,'auto')
    ivals=nan(size(bands,1),2);
    for iFreq=1:size(bands,1)
        dat_tmp=proc_selectChannels(dat,sprintf('*flt%i*',iFreq));
        dat_tmp.clab=clab_load;
        ival= select_timeivalEpo(dat_tmp,'maxIval',opt.maxIval,...
            'startIval',opt.startIval);
        %     if diff(ival)<1000
        %         ival= [mean(ival)-500 mean(ival)+500];
        %     end
        ivals(iFreq,:)=ival;
    end
else
    ivals=[ival;ival];
end
% dat=proc_selectIval(dat,ival);
% [dat, W, A, score]=proc_multiBandSpatialFilter(dat,...
%     {@proc_cspAuto,opt});

W={};
A={};
score={};
fv=[];
for ifreq=1:size(freqs,1)
    epo= proc_selectIval(proc_selectChannels(dat,sprintf('*flt%d',ifreq)),ivals(ifreq,:));
    [fv_i, csp_w_i,csp_a_i,csp_score_i]=proc_sCSPAuto(epo,opt_pickProps(opt, proc_sCSPAuto()));
    fv_i= proc_variance(fv_i);
    fv_i= proc_logarithm(fv_i);
    fv=proc_appendChannels(fv,fv_i);
    W{end+1}=csp_w_i;
    A{end+1}=csp_a_i;
    score{end+1}=csp_score_i;
end

varargout={W,A,score,bands, ivals};
