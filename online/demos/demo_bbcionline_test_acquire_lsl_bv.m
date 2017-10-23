% setup the bbci variable to define the online processing chain
bbci= struct;
bbci.source.acquire_fcn= @bbci_acquire_lsl;

% define the dummy electrode setting
% clab = str_cprintf('%d',1:35);
clab = {'Fp1' 'Fz' 'AF3' 'AF7' 'FT7' 'FC5' 'FC1' 'C3' 'C5' 'TP7' 'CP5' 'P1'...
    'Pz' 'P3' 'P7' 'O1' 'Oz' 'O2' 'P4' 'P8' 'TP8' 'CP6' 'CP2' 'Cz' 'C4' ...
    'C6' 'FT8' 'FC6' 'FC2' 'AF4' 'AF8' 'Fp2' 'Acc1' 'Acc2' 'Acc3'};

% clab = {'Fp1' 'AFz' 'AF3' 'F5' 'F9' '
%     'AF5' 'AF3' 'AF1' 'AFz' 'AF2' 'AF4' 'AF6' ...
%         'F5' 'F3', 'F1' 'Fz' 'F2' 'F4' 'F6' ...
%         'FC7' 'FC5'};
%
% filt1=[];
% [filt1.b ,filt1.a]=butter(5,[1 40 ]/250*2,'bandpass');
% Wps= [40 49]/250*2;
% [n, Ws]= cheb2ord(Wps(1), Wps(2), 3, 50);
% filt1=[];
% [filt1.b, filt1.a]= cheby2(n, 50, Ws);
% filt2=[];
% [filt2.b, filt2.a]= cheby2(5, 50, 1/250*2);
% filtHd= procutil_catFilters(filt1, filt2);

bbci.signal.clab = clab;
% provide clab and markerstreamname to lsl acquire function
% bbci.source.acquire_param = {'clab', clab, 'markerstreamname', 'LiveAmpSN-054206-0110-Sampled-Markers','fs',250};
% bbci.source.acquire_param = {'clab', clab, 'markerstreamname', 'MyMarkerStream','filtHd',filtHd};
bbci.source.acquire_param = {'clab', clab, 'markerstreamname', 'MyMarkerStream','FiltHd',filtHd_HP};
% bbci.source.acquire_param = { 'markerstreamname', 'MyMarkerStream'};

bbci.feature.proc= {@proc_variance, @proc_logarithm};
bbci.feature.ival= [-500 0];
bbci.log.output= 'screen&file';
bbci.log.clock= 1;
bbci.log.markers= 1;

bbci.quit_condition.running_time= 930;
bbci.quit_condition.marker= {'S 255'};

bbci_recordSignals(bbci, 'harmonicsCalibration');


