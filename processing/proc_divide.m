function epo= proc_divide(epo1, epo2, varargin)
%PROC_SUBTRACTREFERENCECLASS - subtracts a reference class in epo_ref from the data in epo.
%
%Synopsis:
% EPO= proc_subtractReferenceClass(EPO, EPO_REF, <OPT>)
%
%Arguments:
% EPO -      data structure of epoched data
% EPO_REF -  data structure of epoched data. Must have the same number of
%               epochs and data per epoch as EPO
%
% OPT struct or property/value list of optional arguments:
%  .SubtractFrom - '*' (default), a string matching one of the classes of
%                   EPO_REF.


% props= {};
% 
% opt = opt_proplistToStruct(varargin{:});
% [opt, isdefault] = opt_setDefaults(opt, props);
% opt_checkProplist(opt, props);        

misc_checkType(epo1, 'STRUCT(x className y)|DOUBLE[1]'); 
misc_checkType(epo2, 'STRUCT(x className y)|DOUBLE[1]'); 

if isstruct(epo1)
    epo=epo1;
elseif isstruct(epo2)
    epo=epo2;
else
    error('Incorrect datatype: One argument has to be a data struct')
end

if isstruct(epo1)&&isstruct(epo2)
    epo.x=epo1.x./epo2.x;
elseif isstruct(epo1)
    epo.x=epo1.x./epo2;
elseif isstruct(epo2)
    epo.x=epo1./epo2.x;
end