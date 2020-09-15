function [out,C] = apply_Wishart(C, x)
% APPLY_CHI2 - apply Chi2 with Pmean update of the bias
%Synopsis:
%  C   [STRUCT with fields 'subC' and 'sublab_unique'] - RSLDA Classifier
%  X   [DOUBLE [ndim nsamples]] - Data
%
%Returns:
%  OUT  - Classifier output
%
%Description:
%  Apply chi2  Classifier
%
%See also:
x=reshape(x,sqrt(size(x,1)),sqrt(size(x,1)),size(x,2));
X1=tprod(pinv(C.C_k(:,:,1)),[1 2],x,[1 2 3]);
X2=tprod(pinv(C.C_k(:,:,2)),[1 2],x,[1 2 3]);

out=C.n/2*(log(det(C.C_k(:,:,1)))-log(det(C.C_k(:,:,2))))+0.5*(trace(X2)-trace(X1));
% if C.lambda==0
%     %out=b'*log(x)+c'*x+A;
%     out=b1'*log(X1)-b2'*log(X2)+c1'*X1-c2'*X2+A;
% else
%     out = nan(1,size(x,2));
%     mu_x=0.5;
%     var_train=2*[C.D1 C.D2].^2./C.k;
%     var_train(:,2)=var_train(:,2)+(C.D2-0.5).^2;
%     var_train(:,1)=var_train(:,1)+(-C.D2+0.5).^2;
%     var_train=mean(var_train,2);
%     %g_train=1./var_train;
%     %g=1;
%     vari=var_train;
%     for xi=1:size(x,2)
%         mu_x=(1-C.lambda)*mu_x+C.lambda*x(:,xi);
%         %g=(1-C.lambda)*g+C.lambda*1./(x(:,xi)-mu_x).^2.*149/150./g_train;
%         vari=(1-C.lambda)*vari+C.lambda*(x(:,xi)-mu_x).^2;
%         g=vari./var_train;
%         out(xi)=sum(b'*log(x(:,xi)./g)+c'*x(:,xi)./g+A);
%     end    
% end
end

