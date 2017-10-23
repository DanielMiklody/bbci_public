function [data, bbci]= bbci_apply(bbci)
%BBCI_APPLY - Apply BBCI classifier to continuously acquired data
%
%Synopsis:
%  [DATA, BBCI]= bbci_apply(BBCI)
%
%To get a description on the structures 'bbci' and 'data', type
%help bbci_apply_structures

% 02-2011 Benjamin Blankertz

scaling=10^5;%2*10^6;
lib = lsl_loadlib();
bbci= bbci_apply_setDefaults(bbci);
[data, bbci]= bbci_apply_initData(bbci);
info = lsl_streaminfo(lib,'BBCI_Ch1','Markers',1,0,'cf_float32','BBCIonlinetoolbox');
outletCh1 = lsl_outlet(info);
run= true;
while run,
  for k= 1:length(bbci.source),
    [data.source(k), data.marker]= ...
        bbci_apply_acquireData(data.source(k), bbci.source(k), data.marker);
    if (max(abs(data.source(k).x(:,25)))/scaling)>1
        scaling=max(abs(data.source(k).x(:,25)));
    end
    outletCh1.push_sample(data.source(k).x(:,25)/scaling);
  end
  if ~all(cellfun(@(x)getfield(x,'running'), {data.source(:).state})),
    break;
  end
  % markers are expected only from source #1
  data.marker.current_time= data.source(1).time;
  for k= 1:length(bbci.signal),
    in= bbci.signal(k).source;
    data.signal(k)= bbci_apply_evalSignal(data.source(in), ...
                                          data.signal(k), ...
                                          bbci.signal(k));
  end
  for ic= 1:length(bbci.control),
    src_list= bbci.control(ic).source_list;
    data.control(ic).time= max([data.source(src_list).time]);
    % if no new data is acquired for this control since last check -> continue
    if data.control(ic).time <= data.control(ic).lastcheck,
      continue;
    end
    events= bbci_apply_evalCondition(data.marker, data.control(ic), ...
                                     bbci.control(ic));
    data.control(ic).lastcheck= data.control(ic).time;
    for ev= 1:length(events),
      data.event= events(ev);
      cfy_list= bbci.control(ic).classifier;
      feat_list= [bbci.classifier(cfy_list).feature];
      for k= feat_list,
        if data.event.time > data.feature(k).time,
          signal= data.signal( bbci.feature(k).signal );
          data.feature(k)= ...
              bbci_apply_evalFeature(signal, bbci.feature(k), data.event);
        end
      end
      for cfy= cfy_list,
        fv= cat(1, data.feature(bbci.classifier(cfy).feature).x);
        data.classifier(cfy)= ...
            bbci_apply_evalClassifier(fv, bbci.classifier(cfy));
      end
      cfy_out= cat(1, data.classifier(cfy_list).x);
      cfy_out=2*(1./(1+exp(cfy_out))-0.5);
      data.control(ic)= ...
          bbci_apply_evalControl(cfy_out, data.control(ic), ...
                                 bbci.control(ic), data.event, data.marker);
      for k= 1:length(bbci.feedback),
        if ismember(ic, bbci.feedback(k).control),
          data.feedback(k)= ...
              bbci_apply_sendControl(data.control(ic).packet, ...
                                     bbci.feedback(k), ...
                                     data.feedback(k));
        end
      end
      bbci_apply_logEvent(data, bbci, ic);
    end
  end
  [bbci, data]= bbci_apply_adaptation(bbci, data);
  run= bbci_apply_evalQuitCondition(data.marker, bbci, data.log.fid);
end
bbci_apply_close(bbci, data);
