function [ sp, stat ] = m3ba_init( spn )
%M3BA_INIT configures and opens the virtual COM port and connects to the m3ba device
%
%   spn: number of serial com port with which to connect
%
%
% RETURNS:
%   sp: Serial Port Object
%   stat:   status struct with info and definitions for driver

% Driver data and Packet Definitions
stat.A = uint8(65);  % Accel Package
stat.E = uint8(69);  % EEG Package
stat.N = uint8(78);  % NIRS Package
stat.I = uint8(73);  % Impedance Package
stat.M = uint8(35);  % System Message (#)
stat.F = uint8(42);  % Feedback Message (*)



% try to find connected cognionics port
  m3ba_ports= instrfind('Tag', 'M3BA');
  if isempty(m3ba_ports),
    
      % Create and Configure Serial Port Object
    sp = serial(['COM' num2str(spn)], 'BaudRate', 230400, 'FlowControl', 'none', ...
        'InputBufferSize', 1024, 'DataBits', 8, 'StopBit', 1, 'Parity', 'none', ...
        'Timeout', 3, 'Name', ['M3BAport' num2str(spn)], 'Tag', ...
        ['M3BA' num2str(spn)], 'Terminator', 'CR'); %
    
  elseif length(m3ba_ports)==1,
    sp= m3ba_ports;
    if strcmp(sp.Status, 'open'),
      fclose(sp);
    end
  else
    error('multiple ports connected to M3BA found');
  end
    
  fopen(sp);
  
    
end

