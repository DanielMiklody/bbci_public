function C = train_Schlauch3m(xTr, yTr, varargin)
% TRAIN_RLDASHRINK - Regularized LDA with automatic shrinkage selection
%
%Synopsis:
%   C = train_RLDAshrink(XTR, YTR)
%   C = train_RLDAshrink(XTR, YTR, OPTS)
%
%Arguments:
%   XTR: DOUBLE [NxM] - Data matrix, with N feature dimensions, and M training points/examples. 
%   YTR: INT [CxM] - Class membership labels of points in X_TR. C by M matrix of training
%                     labels, with C representing the number of classes and M the number of training examples/points.
%                     Y_TR(i,j)==1 if the point j belongs to class i.
%   OPT: PROPLIST - Structure or property/value list of optional
%                   properties. Options are also passed to clsutil_shrinkage.
%     'ExcludeInfs' - BOOL (default 0): If true, training data points with value 'inf' are excluded from XTR
%     'Prior' - DOUBLE (default ones(nClasses, 1)/nClasses): Empirical class priors
%     'StorePrior' - BOOL (default 0): If true, the prior will be stored with the classifier in C.prior
%     'Scaling' - BOOL (default 0): scale projection vector such that the distance between
%        the projected means becomes 2. Scaling only implemented for 2 classes so far. Using Scaling=1 will disable the use of a prior.
%     'StoreMeans' - BOOL (default 0): If true, the classwise means of the feature vectors
%        are stored in the classifier structure C. This can be used, e.g., for bbci_adaptation_pmean
%     'UsePcov' - BOOL (default 0): If true, the pooled covariance matrix is used instead of the average classwise covariance matrix.
%     'StoreCov',  - BOOL (default 0): If true, the covariance matrix will be stored with the classifier in C.cov
%     'StoreInvcov' - BOOL (default 0): If true, the inverse of the covariance matrix is stored in
%        the classifier structure C. This can be used, e.g., for bbci_adaptation_pcovmean
%     'StoreExtinvcov' - BOOL (default 0): If true, the extended inverse of the covariance will be stored with the classifier in C.extinvcov
%
%Returns:
%   C: STRUCT - Trained classifier structure, with the hyperplane given by
%               fields C.w and C.b.  C includes the fields:
%    'w' : weight matrix
%    'b' : FLOAT bias
%    'prior' : (optional) classwise priors
%    'means' :  (optional) classwise means
%    'cov' :  (optional) covariance matrix
%    'invcov' :  (optional) inverse of the covariance matrix
%    'extinvcov' : (optional) extended inverse of the covariance matrix
%
%Description:
%   TRAIN_RLDA trains a regularized LDA classifier on data X with class
%   labels given in LABELS. The shrinkage parameter is selected by the
%   function clsutil_shrinkage.
%
%   References: J.H. Friedman, Regularized Discriminant Analysis, Journal
%   of the Americal Statistical Association, vol.84(405), 1989. The
%   method implemented here is Friedman's method with LAMBDA==1. The
%   original RDA method is implemented in TRAIN_RDAREJECT.
%
%Examples:
%   train_RLDA(X, labels)
%   train_RLDA(X, labels, 'Target', 'D')
%   
%See also:
%   APPLY_SEPARATINGHYPERPLANE, CLSUTIL_SHRINKAGE, 
%   TRAIN_LDA, TRAIN_RDAREJECT

% Benjamin Blankertz
% 12-09-2012: revised to fit with new naming standards and automatic
% opt-type checking (Michael Tangermann)
% 08-2019 Daniel Miklody

props= {'Gamma'      'auto'  'CHAR(auto)|!DOUBLE'
       };


if size(yTr,1)==1, yTr= [yTr<0; yTr>0]; end
nClasses= size(yTr,1);
nChannels=size(xTr,2);

props_shrinkage= clsutil_shrinkage;

if nargin==0,
  %C= opt_catProps(props, props_shrinkage); 
  C= opt_catProps(props_shrinkage); 
  return
end

opt= opt_proplistToStruct(varargin{3:end});
[opt, isdefault]= opt_setDefaults(opt, props);
opt_checkProplist(opt, props, props_shrinkage);

if nargin<3||~isnumeric(varargin{1})
    error('Schlauch3 needs yo Eigenvalues')
end

vars=zeros(1,nChannels,nClasses);
for ii=1:nClasses
    vars(:,:,ii)=var(xTr(:,:,yTr(ii,:)==1),[],3);
end

C.D=varargin{1};
kest=[2*(1-C.D)./vars(:,:,1)',(2*C.D./vars(:,:,2)')];
kest(C.D>0.5,1)=kest(C.D>0.5,2);
kest=kest(:,1);
C.k=kest;

C.w=C.k./2.*(1./(1-C.D)-1./C.D);
C.b=(C.k'/2-1)*log(1./C.D-1);



