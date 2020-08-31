function [out,C] = apply_Chi2(C, x)
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

b1=C.k(:,2)/2;
b2=C.k(:,1)/2;
b=b1-b2;
c1=-0.5.*(C.k(:,2)./C.D);
c2=-0.5.*C.k(:,1)./(1-C.D);
c=c1-c2;
A1=-C.k(:,2)*(log(2))/2-gammaln(C.k(:,2)/2)+...
    (C.k(:,2)/2).*log(C.k(:,2)./C.D);
A2=-C.k(:,1)*(log(2))/2-gammaln(C.k(:,1)/2)+...
    (C.k(:,1)/2).*log(C.k(:,1)./(1-C.D));
A=sum(A1-A2);

if C.lambda==0
    out=b'*log(x)+c'*x+A;
else
    out = nan(1,size(x,2));
    mu_x=0.5;
    var_train=2*[1-C.D C.D].^2./C.k;
    var_train(:,2)=var_train(:,2)+(C.D-0.5).^2;
    var_train(:,1)=var_train(:,1)+(-C.D+0.5).^2;
    var_train=mean(var_train,2);
    %g_train=1./var_train;
    %g=1;
    vari=var_train;
    for xi=1:size(x,2)
        mu_x=(1-C.lambda)*mu_x+C.lambda*x(:,xi);
        %g=(1-C.lambda)*g+C.lambda*1./(x(:,xi)-mu_x).^2.*149/150./g_train;
        vari=(1-C.lambda)*vari+C.lambda*(x(:,xi)-mu_x).^2;
        g=vari./var_train;
        out(xi)=sum(b'*log(x(:,xi)./g)+c'*x(:,xi)./g+A);
    end    
end
end

