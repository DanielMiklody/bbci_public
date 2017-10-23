% % Initiate Connection
% sp=m3ba_init(8);
% pause(0.1);
% 
% % Switch on Feedback 
% m3ba_send_ctrl(sp, 'FE');
% %and switch off NIRS and ACCEl
% % m3ba_send_ctrl(sp, 'MND');
% % m3ba_send_ctrl(sp, 'MAD');
% % configure EEG Channel 1-4 for test signal, 5 to GND and 6 to MVDD
% m3ba_send_ctrl(sp, 'CECAE');
% m3ba_send_ctrl(sp, 'CEC1I5');
% m3ba_send_ctrl(sp, 'CEC1G4');
% m3ba_send_ctrl(sp, 'CEC2I5');
% m3ba_send_ctrl(sp, 'CEC2G8');
% m3ba_send_ctrl(sp, 'CEC3I5');
% m3ba_send_ctrl(sp, 'CEC3G12');
% m3ba_send_ctrl(sp, 'CEC4I5');
% m3ba_send_ctrl(sp, 'CEC4G24');
% m3ba_send_ctrl(sp, 'CEC5I1');
% m3ba_send_ctrl(sp, 'CEC6I3');
% 
% % start measurement
% m3ba_start(sp);
% pause(0.1);
% 
% % Get data for 5 seconds
% disp('starting measurement');
% start=0;
% for i=1:500
%     if (i==1)
%         [pack, buf] = m3ba_getdata(sp, start, false);
%     else
%         [pack, buf] = m3ba_getdata(sp, buf, false);
%     end
%    % wait 10 ms
%    pause(0.010);
% end
% 
% m3ba_stop(sp);
% disp('measurement done');
% m3ba_close(sp);

%Prepare Data Bins
EEG.ch=[0;0;0;0];
EEG.datdim='DOUBLE [#Channels Time]';
EEG.packet=0;
NIR.ch=[0;0;0;0;0;0;0;0];
NIR.datdim='DOUBLE [#Channels(WL 1 and 2) Time]';
NIR.packet=0;
ACC.ch=[0;0;0];
ACC.datdim='DOUBLE [Channels(XYZ) Time]';
ACC.packet=0;

% save received data and feed in to empty the data buffer (for offline read out)
received=buf;
idx=1;
while length(buf) > 100
        %extract data (driver)
        [pbuf, buf]=m3ba_getdata(sp, buf, true);
        
        packet{idx}=pbuf;
        idx=idx+1;
        switch pbuf.type
            case 'EEG'
                EEG.ch=[EEG.ch pbuf.data];
                EEG.packet = [EEG.packet pbuf.nr];
            case 'ACC'
                ACC.ch=[ACC.ch pbuf.data];
                ACC.packet = [ACC.packet pbuf.nr];
            case 'NIR'
                NIR.ch=[pbuf.data NIR.ch];
                NIR.packet = [NIR.packet pbuf.nr];
            case 'MSG'
                disp(pbuf.message)
            case 'FEB'    
                disp(pbuf.feedback);
        end
end

figure
plot(EEG.ch')
figure
plot(ACC.ch')
figure
plot(NIR.ch')
figure
plot(EEG.packet)




