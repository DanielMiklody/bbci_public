function R= procutil_covClasswise(dat, varargin)

props= {'covPolicy', 'average', 'CHAR(average normal)'
        'normalize', 0,         'INT'
        'weight',ones(1,size(dat.y,2)), 'DOUBLE[1 -]'
        'WeightExp',1 ,                 '!DOUBLE[1]'
        'shrinkage',1 ,                  'BOOL'
       };
opt= opt_proplistToStruct(varargin{:});
opt= opt_setDefaults(opt, props);

nClasses= size(dat.y, 1);
nChans= size(dat.x, 2);

if opt.shrinkage
    covfunc=@(x)clsutil_shrinkage(x');
else
    covfunc=@cov;
end
%% calculate classwise covariance matrices
R= zeros(nChans, nChans, nClasses);
switch(lower(opt.covPolicy)),
 case 'average',
  for c= 1:nClasses,
    C= zeros(nChans, nChans);
    idx= find(dat.y(c,:));
    for m= idx,
      C= C + (opt.weight(m)^opt.WeightExp)*covfunc(dat.x(:,:,m));
    end
    R(:,:,c)= C/length(idx);
  end
 case 'normal',
  T= size(dat.x, 1);
  for c= 1:nClasses,
    idx= find(dat.y(c,:));
    x= permute(dat.x(:,:,idx), [1 3 2]);
    x= reshape(x, T*length(idx), nChans);
    R(:,:,c) = covfunc(x);
  end
 otherwise,
  error('covPolicy not known');
end

t= zeros(nClasses, 1);
if opt.normalize,
  for cc= 1:nClasses,
    t(cc)= trace(R(:,:,cc));
  end
  if opt.normalize==1,
    for cc= 1:nClasses,
      R(:,:,cc)= R(:,:,cc)/t(cc);
    end
  elseif opt.normalize==2,
    for cc= 1:nClasses,
      R(:,:,1)= R(:,:,1)/mean(t);
    end
  else
    error('unknown option for ''normalize''.');
  end
end
