function header = load_header(header_filename)
% read the .csv header info from Kistler sensors
% The script takes a single argument - the name of the header file. This
% can be specified as a complete filename or as the base filename (e.g.
% without the suffix _headerinfo.csv
% The function returns a data structure with some of the information in the
% header file:
% header.serial_number      # LabAmp serial numbers in the order recorded
% header.sampling_rates     # sampling rate in samples/second
% header.pressure_sensor_serial_numbers     # serial number of the pressure sensor
% header.temperature_sensor_serial_numbers  # serial number of the temperature sensor

if ~endsWith(header_filename,'_headerinfo.csv')
    header_filename = [header_filename,'_headerinfo.csv'];
end

header = struct();
fh = fopen(header_filename,'r');
section = 'none'; % section of the header file
while ~feof(fh)
    line = fgetl(fh);
    if length(line) >= 18 && strcmp(line(1:18),'Device Information')
        section = 'Device Information';
        % get the labamp serial numbers
    elseif length(line) >= 7 && strcmp(line(1:7),'Channel')
        section = 'Channel';
        this_channel = sscanf(line,'Channel %d',1);
    else
        switch section
            case 'Device Information'
                if startsWith(line,'Serial Number')
                    fields = split(line,',');
                    header.serial_number = cellfun(@(x) str2num(x),fields(2:end));
                    % header.serial_number = sscanf(line,'Serial Number,%d,%d',2);
                elseif startsWith(line,'Sampling Rate')
                    fields = split(line,',');
                    header.sampling_rates = cellfun(@(x) str2num(x),fields(2:end));
                end
            case 'Channel'
                if startsWith(line,'Name')
                    fields = split(line,',');
                    if startsWith(fields{2},'Temperature')
                        sensor_type = 'T';                        
                        sensor_numbers = cellfun(@(x) sscanf(x,'Temperature_%d',1),fields(2:end));
                    elseif startsWith(fields{2},'Pressure')
                        sensor_type = 'P';
                        sensor_numbers = cellfun(@(x) sscanf(x,'Pressure_%d',1),fields(2:end));                        
                    end
                elseif startsWith(line,'Serial Number')
                    fields = split(line,',');
                    serial_numbers = cellfun(@(x) sscanf(x,'%d',1),fields(2:end));                        
                    if sensor_type == 'P'
                        disp(serial_numbers)
                        header.pressure_sensor_serial_numbers(sensor_numbers) = serial_numbers;
                    elseif sensor_type == 'T'
                        header.temperature_sensor_serial_numbers(sensor_numbers) = serial_numbers;

                    end
                end % end serial number
            otherwise
                %.. do nothing
        end % end switch on sections
        
    end
end
