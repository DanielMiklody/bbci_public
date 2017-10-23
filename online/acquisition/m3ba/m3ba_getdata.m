function [ packet, buf_out ] = m3ba_getdata( sp, stat, buf_in, expackets )
%M3BA_GETDATA fetches data from the serial port and adds it to
% buf_in (binary data, byte array). If one complete data package(s) within the 
% binary input buffer data is found, it will be converted from binary to 
% integer data and a complete "packet" as well as the remaining "buf_out" 
% data is returned (use as input for the next iteration);

%   sp: serial port object, from m3ba_init()
%   stat: status struct for driver, from m3ba_init()
%   buf_in: byte array of type uint8
%   expackets: set to true if packets should be extracted while buffer is filled
%       (e.g. in online as opposed to offline data conversion )

  


    % Check whether new data has arrived and serial port is still open, append
    if (sp.BytesAvailable && strcmp(sp.Status, 'open'))
        % cast buffer, just to be safe
        buf = uint8(buf_in);
        received = fread(sp, sp.BytesAvailable, 'uint8');
        % add new data
        buf = [buf_in; received];

    else
        buf= buf_in;
    end
    
    % find end of package flag \t\r\n - 09 0D 0A - 9 13 10
        ind=strfind(buf', [9 13 10]);
        % Package terminator found?
        if (length(ind)>1 && expackets)

            % split off package (first in fifo), put rest in output buffer for next iteration
            packbuf = buf(1:ind(1)-1);
            buf_out = buf(ind(1)+3:end);

            % Do Consistency Check and convert/wrap up packet
            switch packbuf(1)
                case stat.A  % Accel Package
                    if length(packbuf) == 8
                        rdata=swapbytes(typecast(uint8(packbuf(2:7)), 'int16'));
                        % convert to +-2g units
                        packet.data=double(rdata)*0.0039;
                        packet.nr = packbuf(8);
                        packet.type='ACC';
                    else
                        packet.type='ERR';
                    end
                case stat.E  % EEG Package
                    if length(packbuf) == 26
                        % sort channels, sign extension and convert
                        rdata=reshape(packbuf(2:25),[3,8]);
                        z(1:8)=uint8(0);
                        zdata = [rdata; z];
                        for i=1:8
                            % arrange bits and bytes (from 24bit int to 32 bit
                            % int
                            tc(i)=swapbytes(typecast(zdata(:,i),'int32'));
                            % bitshift to the right by one byte and conversion
                            % to volt
                            dat(i)=double(bitshift(tc(i),-8))*4.5/(2^23);
                        end
                        packet.data= [dat(1:4)', dat(5:8)'];
                        % fill packet info
                        packet.type='EEG';
                        packet.nr = packbuf(26);

                    elseif length(packbuf) == 38
                        % sort channels, sign extension and convert
                        rdata=reshape(packbuf(2:37),[3,12]);
                        z(1:12)=uint8(0);
                        zdata = [rdata; z];
                        for i=1:12
                            % arrange bits and bytes (from 24bit int to 32 bit
                            % int
                            tc(i)=swapbytes(typecast(zdata(:,i),'int32'));
                            % bitshift to the right by one byte and conversion
                            % to volt
                            dat(i)=double(bitshift(tc(i),-8))*4.5/(2^23);
                        end
                        packet.data= [dat(1:6)', dat(7:12)'];
                        packet.type='EEG';
                        packet.nr = packbuf(38);
                    else
                        packet.type='ERR';
                    end
                case stat.N  % NIRS Package with 4 channels á 2 wavelengths
                    if length(packbuf) == 26 
                        % sort channels, sign extension and convert
                        rdata=reshape(packbuf(2:25),[3,8]);
                        z(1:8)=uint8(0);
                        zdata = [rdata; z];
                        for i=1:8
                            % arrange bits and bytes (from 24bit int to 32 bit
                            % int
                            tc(i)=swapbytes(typecast(zdata(:,i),'int32'));
                            % bitshift to the right by one byte and conversion
                            % to volt
                            dat(i)=double(bitshift(tc(i),-8))*4.5/(2^23);
                        end
                        packet.data= dat';
                        packet.type='NIR';
                        packet.nr = packbuf(26);
                        %NIRS Package with 4+1 channels á 2 wavelengths
                    elseif length(packbuf) == 32
                        % sort channels, sign extension and convert
                        rdata=reshape(packbuf(2:31),[3,10]);
                        z(1:10)=uint8(0);
                        zdata = [rdata; z];
                        for i=1:10
                            % arrange bits and bytes (from 24bit int to 32 bit
                            % int
                            tc(i)=swapbytes(typecast(zdata(:,i),'int32'));
                            % bitshift to the right by one byte and conversion
                            % to volt
                            dat(i)=double(bitshift(tc(i),-8))*4.5/(2^23);
                        end
                        packet.data= dat';
                        packet.type='NIR';
                        packet.nr = packbuf(32);
                    else
                        packet.type='ERR';
                    end
                case stat.I  % Impedance Package- TBD
                    if length(packbuf) == 49 

                        packet.type='IMP';
                    else
                        packet.type='ERR';
                    end    
                case stat.M  % # Message Package
                    packet.type='MSG';
                    packet.message=char(packbuf(2:end));
                case stat.F  % * Feedback Package
                    packet.type='FEB';
                    packet.feedback=char(packbuf(2:end));
                otherwise
                    packet.type='ERR';
            end

        % Else put all data in output buffer for next iteration
        else
            buf_out=buf;
            packet.type='NUL';
        end

end

