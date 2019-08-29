function [out,C] = apply_RLDAshrinkPmean(C, x)
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

out = nan(1,size(x,2));

for xi=1:size(x,2)
C.b=(1-C.lambda)*C.b-C.lambda*C.w'*x(:,xi);
out(xi)= real( C.w'*x(:,xi) + C.b );
end
end

