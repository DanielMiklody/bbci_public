function C= proc_updatePmean(epo, C, varargin)
%PROC_AVERAGE - Update the LDA bias using the pmean rule
%
%Synopsis:
% EPO= proc_updatePmean(EPO,C, <OPT>)
%
%Arguments:
% EPO -      data structure of epoched data
%            (can handle more than 3-dimensional data, the average is
%            calculated across the last dimension)
% OPT struct or property/value list of optional arguments:
% 'lambda' - the update ratio
%
%Returns:
% C  - updated classifier struct
%
% 08-2019 Daniel Miklody

props= {  'lambda' []   'DOUBLE'};

opt= opt_proplistToStruct(varargin{:});

misc_checkType(epo, 'STRUCT(x clab y)'); 

[opt, isdefault]= opt_setDefaults(opt, props);
opt_checkProplist(opt, props);        
epo = misc_history(epo);

pmean=mean(epo.x,3);
C.b=(1-opt.lambda)*C.b+opt.lambda*pmean;

