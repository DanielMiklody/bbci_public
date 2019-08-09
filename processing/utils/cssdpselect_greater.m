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

if nargin>3
    N_min=varargin{1};
else
    N_min=0;
end

idx= find(abs(score)>value);
if numel(idx)<N_min
    idx2=cspselect_equalPerClass(score,0,N_min/2);
    idx=unique([idx; idx2]);
end
fprintf('%i compponents selected\n',numel(idx))
