function [out,C] = apply_Schlauch3(C, x)
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
out1=C.w'*x+C.b;
% h=x-C.w*(out1)/norm(C.w).^2;
% out=geomean(abs(x-h),1).*sign(out1);
out=geomean(abs(C.w*(out1)/norm(C.w).^2),1).*sign(out1);


%Schlauch3000
% %out1=C.w'*x+C.b;
% h=x-C.w*(C.w'*x+C.b)/norm(C.w).^2;
% out=sum(log(x)-log(h),1);

%out=sum(bsxfun(@times,sign(C.w),log(bsxfun(@times,abs(C.w),log(x_)))));
%hdir=sign(C.w'*x);
%h=Chi2ToGauss(h',C.k);
%x=Chi2ToGauss(x',C.k);
%out= sum(real( x-h )',1);
%B=h(1,:).*x(1,:)-sum(h(2:end,:).*x(2:end,:),1);
%d=acosh(B);
%out=abs(d);%.*;
% h=log(abs(h)).*sign(h);
% x=log(x);
% out= sum(real( x-h ),1);
%out= real( sum(log(bsxfun(@times,C.w,x) + C.b )));
%out= real( sum(Chi2ToGauss(bsxfun(@times,C.w,x) + C.b )));

end

