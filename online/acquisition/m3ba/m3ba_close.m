function [ success ] = m3ba_close( sp )
%M3BA_CLOSE stops the device and closes the serial connection

    m3ba_send_ctrl(sp, 'P');
    m3ba_send_ctrl(sp, 'D');

    fclose(sp);
    
    disp(['COM Port ' sp.Status]);
    if strcmp(sp.Status, 'open')
        success = false;
    else
        success = true;
    end
end

