function [out,C] = apply_RLDAshrinkPmean(C, x, varargin)
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
props= {  'lambda' 1e-5   'DOUBLE'};
opt= opt_proplistToStruct(varargin{:});
[opt, isdefault]= opt_setDefaults(opt, props);
opt_checkProplist(opt, props); 

out = nan(1,size(x,2));

for xi=1:size(x,2)
pmean=mean(x(:,xi));
C.b=(1-opt.lambda)*C.b+opt.lambda*pmean;
out(xi)= real( C.w'*x(:,xi) + C.b );
end
end

