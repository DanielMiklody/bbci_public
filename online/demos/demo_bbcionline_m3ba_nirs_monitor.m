% BBCI DEMO - Testing Matlab feedback using random signals.
%
%  The demo shows, how a Matlab-based feedback can be tested in simulated
%  online mode. As a source of signals, a random signal generator is used
%  (acquire fcn 'bbci_acquire_randomSignals').
%  For this demo, the data processing does not matter. Here, we define
%  a very simple processing chain with a random clasifiers.

%clab= {'NIRS1low','NIRS1high','NIRS2low', 'NIRS2high','NIRS3low', 'NIRS3high','NIRS4low', 'NIRS4high'};
clab= {'ACC1','ACC2','ACC3'};
C= struct('b', 0);
C.w= randn(length(clab), 1);

bbci= struct;
bbci.source.acquire_fcn= @bbci_acquire_m3ba_nirs;
bbci.source.acquire_param= {'clab', clab ,'port',7};
% bbci.source.port=7;
% bbci.source.acquire_fcn= @bbci_acquire_randomSignals;
% bbci.source.acquire_param= {'fs',250,'clab',clab, 'realtime', 2};
% bbci.feature.proc={ {@online_BeerLambert,[],'Baseline',ones(8,1)*8/100,'wavelengths',[750 850]}};

monitor_signalViewer(bbci,'Band',[],'Maximize',true,'CloseOnExit',true,'Range',1)