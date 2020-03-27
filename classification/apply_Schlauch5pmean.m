function [out,C] = apply_Schlauch5pmean(C, x)
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
    out(xi)=C.w'*x(:,xi) +C.b;
%     out(xi)=geomean(abs(C.w(C.w~=0)*(out1)/norm(C.w(C.w~=0)).^2),1).*sign(out1);
    %out(xi)=out1;
end

end


