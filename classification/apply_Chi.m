function [out,C] = apply_Chi(C, x)
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
% x=bsxfun(@plus,x,C.b./C.w/size(x,1));
k1=C.k(:,2);
k2=C.k(:,1);
PSD1=prod(ChiSquarePSD(bsxfun(@times,x,k1./C.D)',k1),2)...
    *prod(k1./C.D);
PSD2=prod(ChiSquarePSD(bsxfun(@times,x,k2./(1-C.D))',k2),2)...
    *prod(k2./(1-C.D));




%%
out=squeeze(PSD1-PSD2)';
end

