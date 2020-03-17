function [out,C] = apply_Schlauch(C, x)
% APPLY_RLDAshrinkPmean - apply RLDA with Pmean update of the bias
%Synopsis:
%  C   [STRUCT with fields 'subC' and 'sublab_unique'] - RSLDA Classifier
%  X   [DOUBLE [ndim nsamples]] - Data
%
%Returns:
%  OUT  - Classifier output
%
%Description:
%  Apply RSLDA Classifier
%
%See also:
%  APPLY_SEPARATINGHYPERPLANE, TRAIN_RSLDASHRINK

out= real( C.w'*x + C.b );

end

