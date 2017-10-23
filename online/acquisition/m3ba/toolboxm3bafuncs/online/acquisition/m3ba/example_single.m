% Initiate Connection
[sp1, stat1]=m3ba_init(8);
% Configure
m3ba_config_preset_single(sp1)

%%
%Prepare Data Bins
EEG1.ch=[0;0;0;0];
EEG1.datdim='DOUBLE [#Channels Time]';
EEG1.packet=0;
NIR1.ch=[0;0;0;0;0;0;0;0];
NIR1.datdim='DOUBLE [#Channels(WL 1 and 2) Time]';
NIR1.packet=0;
ACC1.ch=[0;0;0];
ACC1.datdim='DOUBLE [Channels(XYZ) Time]';
ACC1.packet=0;

% start measurement
m3ba_start(sp1);

% Get data for 5 seconds
disp('starting measurement');
start=0;
it=0;
idx1=1;
tic
while it < 5000
    if (it==0)
        [pack1, buf1] = m3ba_getdata(sp1, stat1, start, true);
    else
        [pack1, buf1] = m3ba_getdata(sp1, stat1, buf1, true);
    end
    
    if ~(strcmp(pack1.type, 'ERR') || strcmp(pack1.type, 'NUL'))
        packet1{idx1}=pack1;
            idx1=idx1+1;
            switch pack1.type
                case 'EEG'
                    EEG1.ch=[EEG1.ch pack1.data];
                    EEG1.packet = [EEG1.packet pack1.nr];
                case 'ACC'
                    ACC1.ch=[ACC1.ch pack1.data];
                    ACC1.packet = [ACC1.packet pack1.nr];
                case 'NIR'
                    NIR1.ch=[pack1.data NIR1.ch];
                    NIR1.packet = [NIR1.packet pack1.nr];
                case 'MSG'
                    disp(pack1.message)
                case 'FEB'    
                    disp(pack1.feedback);
            end
    end
        
   % wait 2 ms
   %pause(0.002);
   it=it+1;
   time(it)=toc;
end

m3ba_stop(sp1);
disp('measurement done');

% figure
% plot(EEG1.ch')
% figure
% plot(EEG1.packet)
figure
plot(ACC1.ch')
figure
plot(NIR1.ch')


%%
% close connection
m3ba_close(sp1);
