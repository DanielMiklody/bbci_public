function [ output_args ] = m3ba_send_ctrl( sp, control )
%M3BA_SEND_CTRL Sends ASCII control string/characters to the device

    % write command + carriage return (0D / 13);
    fwrite(sp, [control 13]);
    % wait 100ms for instrument reply, get data and display
    pause(0.1);
    if (sp.BytesAvailable)
        reply = fscanf(sp, '%s', sp.BytesAvailable);
        disp(char(reply));
    end
end

