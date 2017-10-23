function [ output_args ] = m3ba_start( sp )
%M3BA_START starts the M3BA measurement cycle by sending an "S" command on
%the serial port

    m3ba_send_ctrl(sp, 'S');
    
end

