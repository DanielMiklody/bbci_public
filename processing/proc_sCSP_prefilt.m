function [dat, varargout]= proc_sCSP_prefilt(dat, varargin)
%PROC_CSSDP - Common Spatio-Frequency Decomposition Pattern (CSP) Analysis
%
%Synopsis:
% [DAT, CSSDP_W, CSSDP_A, SCORE]= proc_cssdp_prefilt(DAT, <OPT>);
%
%Arguments:
% DAT - data structure of epoched data
% OPT - struct or property/value list of optional properties:
%  'CovFcn' - Handle of a function to estimate a covariance matrix.
%             Can also be a CELL {@FCN, PARAM1, PARAM2, ...} that
%             provides further arguments to that function.
%             Default: @cov
%             Other options, e.g., @procutil_covShrinkage
%  'ScoreFcn' - Handle of a function that calculates a score of
%               the extracted components.
%               Can also be a CELL {@FCN, PARAM1, PARAM2, ...}.
%               Default: @procutil_score_eigenvalues
%               Other options, e.g., @procutil_score_ratioOfMedians
%  'SelectFcn' - Handle of a function that selects a subset of the
%                extracted components.
%                Can also be a CELL {@FCN, PARAM1, PARAM2, ...}.
%                Default: {@procutil_selectMinMax, 3}
%                Other options, e.g., {@procutil_selectMaxAbs, 6}
%  'Verbose' - Print warnings and other output if larger than 0. Default 1
% 
%Returns:
% DAT   - updated data structure
% CSP_W - CSP projection matrix (spatial filters, in the columns)
% CSP_A - estimated mixing matrix (activation patterns, in the columns)
% SCORE - score of each CSP component
%
%Description:
% calculate common spatial patterns (CSP).
% please note that this preprocessing uses label information. so for
% estimating a generalization error based on this processing, you must
% not apply csp to your whole data set and then do cross-validation
% (csp would have used labels of samples in future test sets).
% you should use the .proc (train/apply) feature in xvalidation, see
% demos/demo_validation_csp
%
%See also demos/demo_validate_csp

props= {'CovFcn'      {@cov}                            '!FUNC|CELL'
        'ScoreFcn'    {@score_eigenvalues}              '!FUNC|CELL'
        'SelectFcn'   {@cspselect_equalPerClass, 3}     '!FUNC|CELL'
        'Verbose'     1                                 'INT'
        'filterOrder'   3                               'INT'
        'ival'  []                               'DOUBLE[- -2]'
        'alpha'      1                              'DOUBLE'
        'chunksize'      10                              'DOUBLE'
       };

if nargin==0,
  dat= props; return
end

misc_checkType(dat, 'STRUCT(x y)'); 
misc_checkType(dat.y, 'DOUBLE[2 -]', 'dat.y');    % two classes only
opt= opt_proplistToStruct(varargin{:});
opt= opt_setDefaults(opt, props);
opt_checkProplist(opt, props);
dat= misc_history(dat);



% Calculate classwise covariance matrices
[covFcn, covPar]= misc_getFuncParam(opt.CovFcn);
nChans= size(dat.x, 2);
nEpo=size(dat.x,3);
C_c= zeros(nChans, nChans, 2);
for k= 1:2,
  X= permute(dat.x(:,:,dat.y(k,:)==1), [1 3 2]);
  X= reshape(X, [], nChans);
  C_c(:,:,k)= covFcn(X, covPar{:});
end


X= permute(dat.x, [1 3 2]);
X= reshape(X, [], nChans);
Ctr=covFcn(X, covPar{:});

X=dat.x(:,:,1:size(dat.x,3)-mod(size(dat.x,3),opt.chunksize));
X=reshape(X,size(X,1),size(X,2),opt.chunksize,size(X,3)/opt.chunksize);
X=permute(X,[1 3 2 4]);
X=reshape(X,[],size(X,3),size(X,4));
C_k=zeros(nChans,nChans,size(X,3));
for k= 1:size(X,3),
  C_temp= covFcn(X(:,:,k), covPar{:});
  if dat.y(1,k)
      C_temp=C_temp-C_c(:,:,1);
  else
      C_temp=C_temp-C_c(:,:,2);
  end
  [~,D_k]=eig(C_temp);
  C_k(:,:,k)=abs(C_temp*sign(D_k));
end
C_k=mean(C_k,3);
% ORIGINAL CODE FOR COMPUTING CSSDP IN CHANNEL SPACE
% % Do actual CSSDP computation as generalized eigenvalue decomposition
[W, D]= eig( C_c(:,:,1)-C_c(:,:,2), C_c(:,:,1)+C_c(:,:,2)+opt.alpha*(C_k) );

% Calculate score for each CSP channel
[scoreFcn, scorePar]= misc_getFuncParam(opt.ScoreFcn);
score= scoreFcn(dat, W, D, scorePar{:});

% Select desired CSSDP filters
[selectFcn, selectPar]= misc_getFuncParam(opt.SelectFcn);
if numel(selectPar{1})>1
    if chanind(dat,'*flt1*')
        selectPar{1}=selectPar{1}(1);
    elseif chanind(dat,'*flt2*')
        selectPar{1}=selectPar{1}(2);
    elseif chanind(dat,'*flt3*')
        selectPar{1}=selectPar{1}(3);
    end
end
idx= selectFcn(score, W, selectPar{:});
W= W(:,idx);
score= score(idx);

% Save old channel labels
if isfield(dat, 'clab'),
  dat.origClab= dat.clab;
end

% Apply CSP filters to time series
dat= proc_linearDerivation(dat, W, 'prependix','sCSP');

% Determine patterns according to [Haufe et al, Neuroimage, 87:96-110, 2014]
% http://dx.doi.org/10.1016/j.neuroimage.2013.10.067
A= Ctr * W / (W'*Ctr*W);

varargout= {W, A, score, Ctr};
