function [data, bbci]= bbci_apply_uni(bbci)
%BBCI_APPLY_UNI - Apply BBCI classifier to continuously acquired data
%
%This is a simplyfied version of bbci_apply. It allows only unimodal
%signals, one type of feature, one classifier etc., i.e., all fields of
%bbci (like bbci.source, bbci.feature) must NOT be arrays, but have length 1.
%Under this restriction, the functionality of this function is the same as
%that of the full bbci_apply.
%Note, that for this function bbci.source.min_blocklength must be > 0
%(setting min_blocklength to 0 makes sense only in the case of multiple
%sources anyway).
%
%Synopsis:
%  [DATA, BBCI]= bbci_apply_uni(BBCI)
%
%See: bbci_apply


% 02-2011 Benjamin Blankertz
% scaling=2*10^5;
scaling=60;
lib = lsl_loadlib();
bbci= bbci_apply_setDefaults(bbci);
[data, bbci]= bbci_apply_initData(bbci);
info = lsl_streaminfo(lib,'BBCI_Ch1','Markers',1,0,'cf_float32','BBCIonlinetoolbox');
outletCh1 = lsl_outlet(info);
% info = lsl_streaminfo(lib,'BBCI_Ch2','Markers',1,0,'cf_float32','BBCIonlinetoolbox');
% outletCh2 = lsl_outlet(info);
% datacollector=[];
run= true;
while run,
    [data.source, data.marker]= ...
        bbci_apply_acquireData(data.source, bbci.source, data.marker);
    outletCh1.push_sample(data.source.x(:,1)/scaling+0.5);
%     outletCh2.push_sample(data.source.x(:,2)/scaling+0.5);
    %   datacollector(end+1:end+size(data.source.x,1),1:size(data.source.x,2))=data.source.x;
    if ~data.source.state.running,
        break;
    end
    data.marker.current_time= data.source.time;
    data.signal= bbci_apply_evalSignal(data.source, data.signal, bbci.signal);
    data.control.time= data.source.time;
    events= bbci_apply_evalCondition(data.marker, data.control, bbci.control);
    data.control.lastcheck= data.marker.current_time;
    for ev= 1:length(events),
        data.event= events(ev);
        data.feature= ...
            bbci_apply_evalFeature(data.signal, bbci.feature, data.event);
        %     datacollector(end+1:end+size(data.feature.x,2),1:size(data.feature.x,1))=data.feature.x;
        data.classifier= ...
            bbci_apply_evalClassifier(data.feature.x, bbci.classifier);
        data.control= ...
            bbci_apply_evalControl(data.classifier.x, data.control, ...
            bbci.control, data.event, data.marker);
        data.feedback= ...
            bbci_apply_sendControl(data.control.packet, bbci.feedback, ...
            data.feedback);
        bbci_apply_logEvent(data, bbci, 1);
    end
    [bbci, data]= bbci_apply_adaptation(bbci, data);
    run= bbci_apply_evalQuitCondition(data.marker, bbci, data.log.fid);
end
% data.datacollector=datacollector;
bbci_apply_close(bbci, data);
