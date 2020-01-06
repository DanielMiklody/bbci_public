function [fv, varargout]= proc_csssp_orig_auto(dat,bands, varargin)
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

props= {
    'CovFcn'      {@cov}                            '!FUNC|CELL'
    'ScoreFcn'    {@score_eigenvalues}              '!FUNC|CELL'
    'SelectFcn'   {@cspselect_equalPerClass, 3}     '!FUNC|CELL'
    'Verbose'     1                                 'INT'
    'filterOrder'   3                               'INT'
    'ival'      'auto'                              'DOUBLE|CHAR'
    'maxIval'    [250 5000]                         'DOUBLE[- 2]'
    'noiseband'   0                                 'INT'
    'startIval'  [750 3500]                         'DOUBLE[- 2]'
    'alpha'      1                                  'DOUBLE'
    'maxFreqIval' [4 45]                             'DOUBLE[- 2]'
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


dat_lap=proc_laplacian(dat);

%% first apply a broad-band filt before using the heuristic
[filt_b, filt_a]= butter(5, [min(bands(:,1)) max(bands(:,2))]/dat.fs*2);
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
end
dat_lap=proc_selectIval(dat_lap,ival);
%% Do banpassfiltering
freqs=cell(1,3);
for iFreq=1:size(bands ,1)
band= select_bandbroad(dat_lap, 'band',bands(iFreq,:));
%band= select_bandnarrow_epo(dat_lap, 'band',bands(iFreq,:),...
    %'bandTopscore',bands(iFreq,:),'areas',{});
if band(1)==band(2)
    band=bands(iFreq,:);
end
freqs{1}=[freqs{1};band];
if opt.noiseband==3
    searchbandsL=[band(1)*0.7 band(1)-0.5];
    searchbandsH=[band(2)+0.5 band(2)*1.3];
    if searchbandsL(1)<opt.maxFreqIval(1)
        searchbandsL(1)=opt.maxFreqIval(1);
    end
%     if serachbandsL(2)>opt.maxival(1)
%         serachbandsL(1)=opt.maxival(1);
%     end
    if searchbandsH(2)>opt.maxFreqIval(2)
        searchbandsH(2)=opt.maxFreqIval(2);
    end
    band_n1=select_bandbroad(dat_lap,...
        'scoreProc',@proc_rSquareInv,...
        'band',searchbandsL);
    band_n2=select_bandbroad(dat_lap,...
        'scoreProc',@proc_rSquareInv,...
        'band',searchbandsH);
    freqs{2}=[freqs{2};band_n1(1) band_n2(2)];
    freqs{3}=[freqs{3};band_n1(2) band_n2(1)];
else
    freqs{2}=[freqs{2};[band(1)*0.75 band(2)*1.25]];
    freqs{3}=[freqs{3};[band(1)*0.95 band(2)*1.05]];
end
% fprintf('Freqs %i: noise1 %d -%d signal %d %d noise2 %d %d\n',iFreq,...
%     band_n1,band,band_n2)
end

[filt_b, filt_a] = butters(4,freqs{1}/dat.fs*2);
[filt_b_pb, filt_a_pb] = butters(4,freqs{2}/dat.fs*2);
[filt_b_sb, filt_a_sb] = butters(4,freqs{3}/dat.fs*2,'stop');
%%
origclab_v=dat.clab;

epo_noise1=proc_filterbank(dat,filt_b_pb,filt_a_pb);
epo_noise=[];
for iFreq=1:size(bands,1)
        epo_noise_tmp=proc_selectChannels(epo_noise1,sprintf('*flt%i*',iFreq));
        epo_noise_tmp.clab=strcat(origclab_v,sprintf('noise_flt%i',iFreq));
        epo_noise_tmp=proc_filt(epo_noise_tmp,filt_b_sb{iFreq},filt_a_sb{iFreq});
        epo_noise=proc_appendChannels(epo_noise,epo_noise_tmp);
end
clear epo_noise1 epo_noise_tmp

dat=proc_filterbank(dat,filt_b,filt_a);
if ischar(opt.ival)&&strcmp(opt.ival,'auto')
    ivals=nan(size(bands,1),2);
    for iFreq=1:size(bands,1)
        dat_tmp=proc_selectChannels(dat,sprintf('*flt%i*',iFreq));
        dat_tmp.clab=origclab_v;
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
%fprintf('ival: %f -%f and %f - %f',ivals(1,1),ivals(1,2),ivals(2,1),ivals(2,2))
dat=proc_appendChannels(dat,epo_noise);

% [dat, W, A, score,Ctr]=proc_multiBandSpatialFilter(dat,...
%     {@proc_csssp_prefilt,'SelectFcn',opt.SelectFcn,'alpha',opt.alpha});
W={};
A={};
score={};
Ctr={};
fv=[];
for ifreq=1:size(freqs{1},1)
    epo= proc_selectIval(proc_selectChannels(dat,...
        sprintf('*flt%d',ifreq)),ivals(ifreq,:));
    [fv_i, csssp_w_i,csssp_a_i,csssp_score_i,Ctr_i]=proc_csssp_prefilt_orig(...
        epo,'SelectFcn',opt.SelectFcn,'alpha',opt.alpha);
    fv_i= proc_variance(fv_i);
    fv_i= proc_logarithm(fv_i);
    fv=proc_appendChannels(fv,fv_i);
    W{end+1}=csssp_w_i;
    A{end+1}=csssp_a_i;
    score{end+1}=csssp_score_i;
    Ctr{end+1}=Ctr_i;
end
varargout= {W, A, score, Ctr, freqs, ivals};
