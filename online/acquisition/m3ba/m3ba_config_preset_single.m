function [ output_args ] = m3ba_config_preset_single( sp1 )
%M3BA_CONFIG_PRESET configures one M3BA module with the settings done below
    
%Enable Feedback messages, NIRS and ACCEL, disable EEG
m3ba_send_ctrl(sp1, 'FE');
m3ba_send_ctrl(sp1, 'MNE');
m3ba_send_ctrl(sp1, 'MAE');
m3ba_send_ctrl(sp1, 'MED');
%m3ba_send_ctrl(sp1, 'MEE');


% Set NIRS channel gain to 12 and LED intensity to 9
m3ba_send_ctrl(sp1, 'CNCG12');
m3ba_send_ctrl(sp1, 'CNL9');
    
    
end

