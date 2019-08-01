function idx= cssdpselect_greater(score, ~, value)
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


idx= abs(score)>value;
fprintf('%i compponents selected\n',sum(idx))
