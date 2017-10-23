function [ output_args ] = m3ba_config_preset_multi( sp1, sp2 )
%M3BA_CONFIG_PRESET configures two m3ba modules with serial port objects
%sp1 and sp2 to run as Master/Slave for concurrent fNIRS/ACCEL acquisition
%note: sp1 is assumed to be master
    
    %Enable Feedback messages
    m3ba_send_ctrl(sp1, 'FE');
    m3ba_send_ctrl(sp2, 'FE');
    
    %Enable NIRS and ACCEL, disable EEG
    m3ba_send_ctrl(sp1, 'MNE');
    m3ba_send_ctrl(sp2, 'MNE');
    m3ba_send_ctrl(sp1, 'MAE');
    m3ba_send_ctrl(sp2, 'MAE');
    m3ba_send_ctrl(sp1, 'MED');
    m3ba_send_ctrl(sp2, 'MED');
    
    % Set NIRS channel gain to 12 and LED intensity to 9
    m3ba_send_ctrl(sp1, 'CNCG12');
    m3ba_send_ctrl(sp2, 'CNCG12');
    m3ba_send_ctrl(sp1, 'CNL9');
    m3ba_send_ctrl(sp2, 'CNL9');
    
    % Master Configuration for sp1
    m3ba_send_ctrl(sp1, 'MMI1');
    m3ba_send_ctrl(sp1, 'MML1');
    m3ba_send_ctrl(sp1, 'MMN1');
    m3ba_send_ctrl(sp1, 'MMM');
    %Slave Configuration for sp2
    m3ba_send_ctrl(sp1, 'MMI2');
    m3ba_send_ctrl(sp1, 'MML1');
    m3ba_send_ctrl(sp1, 'MMN1');
    m3ba_send_ctrl(sp1, 'MMS');
    
    
end

