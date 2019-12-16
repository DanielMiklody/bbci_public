function [dat, varargout]= proc_sCSPAuto_orig(dat, varargin)
%PROC_CSPAUTO - Common Spatial Pattern Analysis with Auto Filter Selection
%
%Synopsis:
% [DAT, CSP_W, CSP_A, CSP_EIG]= proc_cspAuto(DAT, <OPT>);
% [DAT, CSP_W, CSP_A, CSP_EIG]= proc_cspAuto(DAT, NPATS);
%
%Arguments:
% DAT    - data structure of epoched data
% NPATS  - number of patterns per class to be calculated, 
%          default nChans/nClasses.
% OPT - struct or property/value list of optional properties:
%  .patterns - 'all' or matrix of given filters or number of filters. 
%      Selection depends on opt.selectPolicy.
%  .selectPolicy - determines how the patterns are selected. Valid choices are 
%      'all', 'equalperclass', 
%      'maxvalues' (attention: possibly not every class is represented!)
%      'maxvalueswithcut' (like maxvalues, but components with less than
%      half the maximal score are omitted),
%      'directorscut' (default, heuristic, similar to maxvalueswithcut, 
%      but every class will be represented), 
%      'matchfilters' (in that case opt.patterns must be a matrix),
%      'matchpatterns' (not implemented yet)
%  .covPolicy - 'normal', 'average' (default) or manually defined cov matrices
%      of size [nchans, nchans, 2] . 'average' calculates the average of the 
%      single-trial covariance matrices.
%  .score - 'eigenvalues', 'medianvar' (default) or 'auc' to
%       determine the components by the respective score
%
%Returns:
% DAT    - updated data structure
% CSP_W  - CSP projection matrix (spatial filters, in the columns)
% CSP_A  - estimated mixing matrix (activation patterns, in the columns)
% CSP_EIG- eigenvalue score of CSP projections 
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

% Author(s): Benjamin Blankertz
props= { 'CovFcn'      {@cov}                            '!FUNC|CELL'
         'patterns'     3           'INT|CHAR'
         'score'        'medianvar' '!CHAR(eigenvalues medianvar auc)'
         'covPolicy'    'average'   'CHAR|DOUBLE[- - 2]'
         'scaling'      'none'      'CHAR'
         'normalize'    0           'BOOL'
         'selectPolicy' 'directorscut'  'CHAR'
         'weight'       []        'DOUBLE'
         'weightExp'    1           'BOOL'
         'alpha'      1                              'DOUBLE'
         'chunksize'      10                              'DOUBLE'
         };

if nargin==0,
  dat = props; return
end

dat = misc_history(dat);
misc_checkType(dat, 'STRUCT(x clab y)'); 
if length(varargin)==1 && isnumeric(varargin{1}),
  opt= struct('patterns', varargin{1});
else
  opt= opt_proplistToStruct(varargin{:});
end
[opt, isdefault]= opt_setDefaults(opt, props);
opt_checkProplist(opt, props);
if isempty(opt.weight)
    opt.weight=ones(1,size(dat.y,2)) ;
end
[T, nChans, nEpochs]= size(dat.x);

if size(dat.y,1)~=2,
  error('this function works only for 2-class data.');
end


%% Calculate classwise covariance matrices
[covFcn, covPar]= misc_getFuncParam(opt.CovFcn);
nChans= size(dat.x, 2);
nEpo=size(dat.x,3);
C_c= zeros(nChans, nChans, 2);
for k= 1:2,
  X= permute(dat.x(:,:,dat.y(k,:)==1), [1 3 2]);
  X= reshape(X, [], nChans);
  C_c(:,:,k)= covFcn(X, covPar{:});
  C_c(:,:,k)= C_c(:,:,k)/trace(C_c(:,:,k));
end

%% do actual sCSP calculation as generalized eigenvalues
C_k=zeros(nChans,nChans,2);
for k=1:2
    X=dat.x(:,:,dat.y(k,:)==1);
    X=X(:,:,1:size(X,3)-mod(size(X,3),opt.chunksize));
    X=reshape(X,size(X,1),size(X,2),opt.chunksize,size(X,3)/opt.chunksize);
    X=permute(X,[1 3 2 4]);
    X=reshape(X,[],size(X,3),size(X,4));
    C_l=zeros(nChans,nChans,size(X,3));
    for l= 1:size(X,3),
        C_temp=covFcn(X(:,:,l), covPar{:})-C_c(:,:,k);
        [V,D_k]=eig(C_temp);
        %C_l(:,:,l)=abs(C_temp*sign(D_k));
        C_l(:,:,l)=V*abs(D_k)*V';
    end
    C_k(:,:,k)=mean(C_l,3);
end
C_k=sum(C_k,3);
C_k=C_k/trace(C_k);
% ORIGINAL CODE FOR COMPUTING CSSDP IN CHANNEL SPACE
% % Do actual CSSDP computation as generalized eigenvalue decomposition
[W, D]= eig( C_c(:,:,1)-C_c(:,:,2), C_c(:,:,1)+C_c(:,:,2)+opt.alpha*(C_k) );

%% calculate score for each CSP channel
switch(lower(opt.score)),
 case 'eigenvalues',
  score= diag(D);
 case 'medianvar',
  fv= proc_linearDerivation(dat, W);
  fv= proc_variance(fv);
  score= zeros(nChans, 1);
  c1= find(fv.y(1,:));
  c2= find(fv.y(2,:));
  for kk= 1:nChans,
    v1= median(fv.x(1,kk,c1),3);
    v2= median(fv.x(1,kk,c2),3);
    score(kk)= v2/(v1+v2);
  end
 case 'auc',
  fv= proc_linearDerivation(dat, W);
  fv= proc_variance(fv);
  fv= proc_aucValues(fv);
  score= -fv.x;
 otherwise,
  error('unknown option for score. It should be ''eigenvalues'', ''medianvar'' or ''auc'' ');
end
%% force score to be a column vector
score= score(:);

%% select patterns
if ischar(opt.patterns) & strcmpi(opt.patterns, 'all'),
  fi= 1:nChans;
elseif ischar(opt.patterns) & strcmpi(opt.patterns, 'auto'),
  if ~strcmpi(opt.selectPolicy, 'maxvalues'),
    score= max(score, 1-score);
  end
  perc= stat_percentiles(score, [20 80]);
  thresh= perc(2) + diff(perc);
  fi= find(score>thresh);
else
  switch(lower(opt.selectPolicy)),
   case 'all',
    fi= 1:nChans;
   case 'equalperclass',
    [dd,di]= sort(score);
    fi= [di(1:opt.patterns); di(end:-1:nChans-opt.patterns+1)];
   case 'minclass2',
    [dd,di]= sort(score);
    fi= di(1:opt.patterns);
   case 'minclass1',
    [dd,di]= sort(score);
    fi= di(end:-1:nChans-opt.patterns+1);
   case 'bestvalues',
    [dd,di]= sort(min(score, 1-score));
    fi= di(1:opt.patterns);
   case 'maxvalues',
    [dd,di]= sort(-score);
    fi= di(1:opt.patterns);
   case 'maxvalueswithcut',
    score= score/max(score);
    [dd,di]= sort(-score);
    iMax= 1:opt.patterns;
    iCut= find(-dd>=0.5);
    idx= intersect(iMax, iCut,'legacy');
    fi= di(idx);
   case 'directorscut',
    if ismember(opt.score, {'eigenvalues','medianvar'},'legacy'),
      absscore= 2*(max(score, 1-score)-0.5);
      [dd,di]= sort(score);
      Nh= floor(nChans/2);
      iC1= find(ismember(di, 1:Nh,'legacy'));
      iC2= flipud(find(ismember(di, [nChans-Nh+1:nChans],'legacy')));
      iCut= find(absscore(di)>=0.66*max(absscore));
      idx1= [iC1(1); intersect(iC1(2:opt.patterns), iCut,'legacy')];
      idx2= [iC2(1); intersect(iC2(2:opt.patterns), iCut,'legacy')];
      fi= di([idx1; flipud(idx2)]);
    else
      score= score/max(score);
      [dd,di]= sort(-score);
      Nh= floor(nChans/2);
      iC1= find(ismember(di, 1:Nh,'legacy'));
      iC2= find(ismember(di, [nChans-Nh+1:nChans],'legacy'));
      iCut= find(-dd>=0.5);
      idx1= [iC1(1); intersect(iC1(2:opt.patterns), iCut,'legacy')];      
      idx2= [iC2(1); intersect(iC2(2:opt.patterns), iCut,'legacy')];
      fi= di([idx1; idx2]);
    end
   case 'matchpatterns',  %% to be implemented
    error('to be implemented');
   case 'matchfilters',  %% greedy, not well implemented
    fi= zeros(1,size(opt.patterns,2));
    for ii= 1:size(opt.patterns,2),
      v1= opt.patterns(:,ii);
      v1= v1/sqrt(v1'*v1);
      sp= -inf*ones(1,nChans);
      for jj= 1:nChans,
        if ismember(jj, fi,'legacy'), continue; end
        v2= W(:,jj);
        sp(jj)= abs(v1'*v2/sqrt(v2'*v2));
      end
      [mm,mi]= max(sp);
      fi(ii)= mi;
    end
   otherwise,
    error('unknown selectPolicy');
  end
end

Wp= W(:,fi);
la= score(fi);

%% optional scaling of CSP filters to make solution unique
switch(lower(opt.scaling)),
 case 'maxto1',
  for kk= 1:size(Wp,2),
    [mm,mi]= max(abs(Wp(:,kk)));
    Wp(:,kk)= Wp(:,kk) / Wp(mi,kk);
  end
 case 'none',
 otherwise
  error('unknown scaling');
end


%% save old channel labels
if isfield(dat, 'clab')  && ~isfield(dat, 'origClab'),
  dat.origClab= dat.clab;
end

%% apply CSP filters to time series
dat= proc_linearDerivation(dat, Wp, 'prependix','sCSP');

%% arrange optional output arguments
if nargout>1,
    varargout{1}= Wp;
    if nargout>2,
        A= pinv(W)'; % return patterns in the columns of A
        varargout{2}= A(:,fi);
        if nargout>3,
            varargout{3}= la;
            if nargout>4,
                varargout{4}= fi;
            end
        end
    end
end
