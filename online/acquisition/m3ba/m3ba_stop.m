function [ output_args ] = m3ba_stop( sp )
%M3BA_STOP stops the M3BA measurement cycle by sending an "P" command on
%the serial port

    m3ba_send_ctrl(sp, 'P');
end

