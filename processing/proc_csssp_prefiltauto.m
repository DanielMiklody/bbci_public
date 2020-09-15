function [dat, varargout]= proc_csssp_prefiltauto(dat, varargin)
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

%% Do banpassfiltering
dat_lap=proc_laplacian(dat);
band1= select_bandnarrow_epo(dat_lap, 'band',[5 7],...
    'bandTopscore',[5 7],'areas',{});
band2= select_bandnarrow_epo(dat_lap, 'band',[7.5 14],...
    'bandTopscore',[7.5 14],'areas',{});
if band1(1)==band1(2)
    band1=[5 7];
end
if band2(1)==band2(2)
    band2=[7.5 14];
end
freqs={[band1;band2],...
    [band1(1)*0.75 band1(2)*1.25;band2(1)*0.75 band2(2)*1.25],...
    [band1(1)*0.95 band1(2)*1.05;band2(1)*0.95 band2(2)*1.05]};
[filt_b, filt_a] = butters(4,freqs{1}/dat.fs*2);
[filt_b_pb, filt_a_pb] = butters(4,freqs{2}/dat.fs*2);
[filt_b_sb, filt_a_sb] = butters(4,freqs{3}/dat.fs*2,'stop');

origclab_v=dat.clab;

epo_noise=proc_filterbank(dat,filt_b_pb,filt_a_pb);

epo_noise=proc_filterbank(epo_noise,filt_b_sb,filt_a_sb);

epo_noise=proc_selectChannels(epo_noise,'*flt1_flt1*','*flt2_flt2*','*flt3_flt3*');

epo_noise.clab(1:numel(origclab_v))=strcat(origclab_v,'noise_flt1');
epo_noise.clab(numel(origclab_v)+1:end)=strcat(origclab_v,'noise_flt2');

dat=proc_filterbank(dat,filt_b,filt_a);
dat=proc_appendChannels(dat,epo_noise);

[dat, W, A, score,Ctr]=proc_multiBandSpatialFilter(dat,...
    {@proc_csssp_prefilt,'SelectFcn',opt.SelectFcn,'alpha',opt.alpha});
        
varargout= {W, A, score, Ctr, freqs};