function [lowband,highband]= select_noisebandbroad(dat, sigband, varargin)
%band= select_noisebandbroad(cnt, mrk, ival, <opt>)
%band= select_noisebandbroad(epo, <opt>)

motor_areas= {{'FC5,3','CFC5,3','C5,3','CCP5,3','CP5,3'},
              {'FC1-2','CFC1-2','C1-2','CCP1-2','CP1-2'}, 
              {'FC4,6','CFC4,6','C4,6','CCP4,6','CP4,6'}};
          
done_laplace= regexpi(dat.clab,'^\w+ lap\w*');
done_laplace= any(cell2mat(done_laplace));


props= {'band'          [5 45]
        'scoreProc'    @proc_rSquare
        'areas'         motor_areas
        'doLaplace'    ~done_laplace
        'smoothSpectrum'   1        
        'threshold' 0.5};

if nargin==0,
  band = props; return
end


misc_checkType(dat, 'STRUCT(x clab fs)'); 
misc_checkType(sigband, 'DOUBLE[2]'); 

if ~isfield(dat,'y')&&~isstruct(varargin{1})
    error('Specify markers and time interval for unepoched data!')
elseif ~isfield(dat,'y')&&isstruct(varargin{1})
    mrk=varargin{1};
    ival=varargin{2}; 
    misc_checkType(mrk, 'STRUCT(time)'); 
    misc_checkType(ival, 'DOUBLE[2]');   
    opt= opt_proplistToStruct(varargin{3:end});
else
    opt= opt_proplistToStruct(varargin{:});
end



[opt, isdefault]= opt_setDefaults(opt, props);
opt_checkProplist(opt, props);

if isempty(opt.areas),
  opt.areas= {dat.clab};
end

[score_fcn, score_param]= misc_getFuncParam(opt.scoreProc);

if opt.doLaplace,
  dat= proc_laplacian(dat);
end

dat= proc_selectChannels(dat, cat(2, opt.areas{:}));
if ~isfield(dat,'y')
    dat= proc_segmentation(dat, mrk, ival);
end
spec= proc_spectrum(dat, opt.band);
%%spec= proc_spectrum(spec, opt.band, 'db_scaled',0);  %% isn't this better?
score= score_fcn(spec, score_param{:});

if opt.smoothSpectrum,
  score= proc_movingAverage(score, 3000/dat.fs, 'method','centered', 'window',[.5 1 .5]');
end

%% calculate score for each channel and choose good channels
%% (one from each area)
chanscore= sqrt(sum(score.x.^2, 1));
for aa= 1:length(opt.areas),
  ci= chanind(score, opt.areas{aa});
  [mm,mi]= max(chanscore(ci));
  chansel(aa)= ci(mi);
end

%% build  score accordingly
freqscore= score.x(:,chansel);
freqscore= mean(freqscore, 2);

sigscore=sum(freqscore(find(score.t==sigband(1)):find(score.t==sigband(2))));

lnbub=find(score.t==sigband(1))-1;
unblb=find(score.t==sigband(2))+1;

lnblb=lnbub-1;
unbub=unblb+1;
while (lnblb>1)&&(unbub<numel(freqscore))&&...
    (sum(freqscore([lnblb:lnbub unblb:unbub]))<opt.threshold*sigscore)
    lnblb=lnblb-1;
    unbub=unbub+2;
end

lowband=score.t([lnblb lnbub]);
highband=score.t([unblb unbub]);
% avspec=proc_average(spec);
% figure
% subplot(311)
% plot(score.t,freqscore)
% subplot(312)
% plot(spec.t,squeeze(mean(avspec.x(:,chansel,:),2)))
% subplot(313)
% plot(score.t,((score.t>=lowband(1))&(score.t<=lowband(2)))|((score.t>=highband(1))&(score.t<=highband(2))))
