function out = apply_NCC(C, x)
% APPLY_RSLDA - apply fcn for each subclass
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

out = C.w'*x+C.b;
%check if there is a NaN... If so, there had been a data point with a
%subclass label which was not specified in the RSLDA classifier
if any(isnan(out))
    warning('at least one data point is nan!')
end

end

