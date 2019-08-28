function idx= cssdpselect_greater(score, ~, value,varargin)
%CSSDPSELECT_LARGER1 - Select components with an absolute eigenvalue
%greater than 1
%
%Synopsis:
% CI = select_fixedNumber(SCORE, ~, NCOMP)
%
%Arguments:
% SCORE  - score of components
% NCOMPS - number of components (total number depends on MODE)
% 
%Returns:
% CI     - index of components
%
%See also processing/proc_csp
props= {  'altfunc'  {@cspselect_fixedNumber, 0,'absolutemax'}     '!FUNC|CELL'
          'Nmin'     0                                           'INT'
          'Verbose'  1                                           'INT'
            };

if nargin==0,
  out = props; return
end
opt= opt_proplistToStruct(varargin{:});
[opt, isdefault]= opt_setDefaults(opt, props);
opt_checkProplist(opt, props);        


idx= find(abs(score)>value);
if numel(idx)<opt.Nmin
    [selectFcn, selectPar]= misc_getFuncParam(opt.altfunc);
    idx2= selectFcn(score, 0, selectPar{:});   
    idx=unique([idx; idx2]);
end
if opt.Verbose
    fprintf('%i compponents selected\n',numel(idx))
end
