% Initiate Connection
sp1=m3ba_init(8);
sp2=m3ba_init(4);
pause(0.1);

% config for multi module preset
m3ba_config_preset_multi(sp1, sp2);

% start measurement, slave first
m3ba_start(sp2);
m3ba_start(sp1);

% Get data for 5 seconds
disp('starting measurement');
start=0;
for i=1:500
    if (i==1)
        [pack1, buf1] = m3ba_getdata(sp1, start, false);
        [pack2, buf2] = m3ba_getdata(sp2, start, false);
    else
        [pack1, buf1] = m3ba_getdata(sp1, buf1, false);
        [pack2, buf2] = m3ba_getdata(sp2, buf2, false);
    end
   % wait 10 ms
   pause(0.010);
end

m3ba_stop(sp1);
pause(0.01)
m3ba_stop(sp2);
disp('measurement done');
m3ba_close(sp1);
m3ba_close(sp2);

%Prepare Data Bins
EEG1.ch=[0;0;0;0];
EEG1.datdim='DOUBLE [#Channels Time]';
EEG1.packet=0;
NIR1.ch=[0;0;0;0;0;0;0;0;0;0];
NIR1.datdim='DOUBLE [#Channels(WL 1 and 2) Time]';
NIR1.packet=0;
ACC1.ch=[0;0;0];
ACC1.datdim='DOUBLE [Channels(XYZ) Time]';
ACC1.packet=0;
EEG2.ch=[0;0;0;0];
EEG2.datdim='DOUBLE [#Channels Time]';
EEG2.packet=0;
NIR2.ch=[0;0;0;0;0;0;0;0];
NIR2.datdim='DOUBLE [#Channels(WL 1 and 2) Time]';
NIR2.packet=0;
ACC2.ch=[0;0;0];
ACC2.datdim='DOUBLE [Channels(XYZ) Time]';
ACC2.packet=0;

% save received data and feed in to empty the data buffer (for offline read out)
received1=buf1;
received2=buf2;
idx1=1;
idx2=1;
while (length(buf1) > 100) || (length(buf2) > 100)
        %extract data (driver) for sp1
        [pbuf1, buf1]=m3ba_getdata(sp1, buf1, true);
        
        packet1{idx1}=pbuf1;
        idx1=idx1+1;
        switch pbuf1.type
            case 'EEG'
                EEG1.ch=[EEG1.ch pbuf1.data];
                EEG1.packet = [EEG1.packet pbuf1.nr];
            case 'ACC'
                ACC1.ch=[ACC1.ch pbuf1.data];
                ACC1.packet = [ACC1.packet pbuf1.nr];
            case 'NIR'
                NIR1.ch=[pbuf1.data NIR1.ch];
                NIR1.packet = [NIR1.packet pbuf1.nr];
            case 'MSG'
                disp(pbuf1.message)
            case 'FEB'    
                disp(pbuf1.feedback);
        end
        
        %extract data (driver) for sp2
        [pbuf2, buf2]=m3ba_getdata(sp2, buf2, true);
        
        packet2{idx2}=pbuf2;
        idx2=idx2+1;
        switch pbuf2.type
            case 'EEG'
                EEG2.ch=[EEG2.ch pbuf2.data];
                EEG2.packet = [EEG2.packet pbuf2.nr];
            case 'ACC'
                ACC2.ch=[ACC2.ch pbuf2.data];
                ACC2.packet = [ACC2.packet pbuf2.nr];
            case 'NIR'
                NIR2.ch=[pbuf2.data NIR2.ch];
                NIR2.packet = [NIR2.packet pbuf2.nr];
            case 'MSG'
                disp(pbuf2.message)
            case 'FEB'    
                disp(pbuf2.feedback);
        end
end

% figure
% plot(EEG1.ch')
% figure
% plot(EEG1.packet)
figure
plot(ACC1.ch')
figure
plot(NIR1.ch')
figure
plot(ACC2.ch')
figure
plot(NIR2.ch')




