function [fv, varargout]= proc_csssp_auto(dat,bands, varargin)
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
    'SelectFcn'   {@cspselect_equalPerClass, 3}     '!FUNC|CELL'
    'Verbose'     1                                 'INT'
    'filterOrder'   3                               'INT'
    'ival'  []                               'DOUBLE[- -2]'
    'alpha'      1                              'DOUBLE'
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
ival= select_timeivalEpo(dat_flt);
% if diff(ival)<1000
%     ival= [mean(ival)-500 mean(ival)+500];
% end           
dat_lap=proc_selectIval(dat_lap,ival);
%% Do banpassfiltering
band1= select_bandnarrow_epo(dat_lap, 'band',bands(1,:),...
    'bandTopscore',bands(1,:),'areas',{});
if band1(1)==band1(2)
    band1=bands(1,:);
end
if size(bands,1)>1
    band2= select_bandnarrow_epo(dat_lap, 'band',bands(2,:),...
        'bandTopscore',bands(2,:),'areas',{});
    if band2(1)==band2(2)
        band2=bands(2,:);
    end
    freqs={[band1;band2],...
        [band1(1)*0.75 band1(2)*1.25;band2(1)*0.75 band2(2)*1.25],...
        [band1(1)*0.95 band1(2)*1.05;band2(1)*0.95 band2(2)*1.05]};
else
    freqs={[band1],...
        [band1(1)*0.75 band1(2)*1.25],...
        [band1(1)*0.95 band1(2)*1.05]};
end
[filt_b, filt_a] = butters(4,freqs{1}/dat.fs*2);
[filt_b_pb, filt_a_pb] = butters(4,freqs{2}/dat.fs*2);
[filt_b_sb, filt_a_sb] = butters(4,freqs{3}/dat.fs*2,'stop');
%%
origclab_v=dat.clab;

epo_noise=proc_filterbank(dat,filt_b_pb,filt_a_pb);

epo_noise=proc_filterbank(epo_noise,filt_b_sb,filt_a_sb);

epo_noise=proc_selectChannels(epo_noise,'*flt1_flt1*','*flt2_flt2*');

epo_noise.clab(1:numel(origclab_v))=strcat(origclab_v,'noise_flt1');
if size(bands,1)>1
    epo_noise.clab(numel(origclab_v)+1:end)=strcat(origclab_v,'noise_flt2');
end
dat=proc_filterbank(dat,filt_b,filt_a);

ivals=nan(size(bands,1),2);
for iFreq=1:size(bands,1)
    dat_tmp=proc_selectChannels(dat,sprintf('*flt%i*',iFreq));
    dat_tmp.clab=origclab_v;
    ival= select_timeivalEpo(dat_tmp);
%     if diff(ival)<1000
%         ival= [mean(ival)-500 mean(ival)+500];        
%     end
    ivals(iFreq,:)=ival;
end
%fprintf('ival: %f -%f and %f - %f',ivals(1,1),ivals(1,2),ivals(2,1),ivals(2,2))
dat=proc_appendChannels(dat,epo_noise);

% [dat, W, A, score,Ctr]=proc_multiBandSpatialFilter(dat,...
%     {@proc_csssp_prefilt,'SelectFcn',opt.SelectFcn,'alpha',opt.alpha});
W={};
A={};
score={};
Ctr={};
for ifreq=1:size(freqs{1},1)
    epo= proc_selectIval(proc_selectChannels(dat,sprintf('*flt%d',ifreq)),ivals(ifreq,:));
    [fv_i, csssp_w_i,csssp_a_i,csssp_score_i,Ctr_i]=proc_csssp_prefilt(epo,'SelectFcn',opt.SelectFcn,'alpha',opt.alpha);
    fv_i= proc_variance(fv_i);
    fv_i= proc_logarithm(fv_i);
    if ifreq==1
        fv=fv_i;
    else
        fv=proc_appendChannels(fv,fv_i);
    end
    W{end+1}=csssp_w_i;
    A{end+1}=csssp_a_i;
    score{end+1}=csssp_score_i;
    Ctr{end+1}=Ctr_i;
end
varargout= {W, A, score, Ctr, freqs, ivals};
