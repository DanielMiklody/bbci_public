function [fv]= proc_csssp_orig_onlineNoiseupdateauto(dat, W,score,Ctr,freqs, ivals, varargin)
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
        'updateStationarity' 1                           'INT'
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

[filt_b, filt_a] = butters(4,freqs{1}/dat.fs*2);
[filt_b_pb, filt_a_pb] = butters(4,freqs{2}/dat.fs*2);
[filt_b_sb, filt_a_sb] = butters(4,freqs{3}/dat.fs*2,'stop');

origclab_v=dat.clab;

epo_noise=proc_filterbank(dat,filt_b_pb,filt_a_pb);

epo_noise=proc_filterbank(epo_noise,filt_b_sb,filt_a_sb);

epo_noise=proc_selectChannels(epo_noise,'*flt1_flt1*','*flt2_flt2*');

epo_noise.clab(1:numel(origclab_v))=strcat(origclab_v,'noise_flt1');
if size(freqs{1},1)>1
    epo_noise.clab(numel(origclab_v)+1:end)=strcat(origclab_v,'noise_flt2');
end
dat=proc_filterbank(dat,filt_b,filt_a);
dat=proc_appendChannels(dat,epo_noise);
% 
% [dat]=proc_multiBandFunction(dat,{@proc_csssp_onlineNoiseupdate,...
%     'lambda',opt.lambda,'alpha',opt.alpha},W,score,Ctr);

for ifreq=1:size(freqs{1},1)        
        fv_i= proc_selectIval(proc_selectChannels(dat,sprintf('*flt%d',ifreq)), ivals(ifreq,:));
        fv_i=proc_csssp_orig_onlineNoiseupdate(fv_i,W{ifreq},score{ifreq},Ctr{ifreq},'lambda',opt.lambda,'alpha',opt.alpha);
        fv_i= proc_variance(fv_i);
        fv_i= proc_logarithm(fv_i);
        if ifreq==1
            fv=fv_i;
        else
            fv=proc_appendChannels(fv,fv_i);
        end
    end

