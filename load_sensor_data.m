function [header,P,T] = load_sensor_data(filename,calibration_table)
% load data from kistler labamps
% INPUT
% filename should be the filename, with or without the trailing .bin
% the calibration table should look something like this. Each row
% corresponds to one sensor.
% calibration_table = [
% SSN    % P offset (bar)   T_offset     T_slope
% ];

header = load_header(filename);
if ~endsWith(filename,'.bin');
    filename = [filename '.bin'];
end

% read binary
fd=fopen(filename,"rb");

% header
FileVersion=fread(fd,1,'int32');
Fs=fread(fd,1,'float');
DevCount=fread(fd,1,'uint32');

for i=1:DevCount
    DevID(i)=fread(fd,1,'int32');
    SNL(i)=fread(fd,1,'uint32');
    SN(i,:)=fread(fd,SNL(i),'*char');
    NameL(i)=fread(fd,1,'uint32');
    Name(i,:)=fread(fd,NameL(i),'*char');
    NumEnChan(i)=fread(fd,1,'uint32');

    for j=1:NumEnChan(i)
        ChanNum(i,j)=fread(fd,1,'int32');
    end
end

% count the number of aggregated data frames
data_start = ftell(fd);
fseek(fd,0,"eof");
file_size = ftell(fd);
fseek(fd,data_start,"bof");
ftell(fd)
c=1;
Trel=0;
nframe = 0;
Nt_total = 0; % total number of time steps
status = true;
while status
    for i=1:DevCount
        TSec=fread(fd,1,'uint64');
        TNsec=fread(fd,1,'uint32');
        NS=fread(fd,1,'uint32');
        Nt=NS/NumEnChan(i);

        % d=fread(fd,[NumEnChan(i),Nt],'float32')'; % [number of enabled channels, number of times]
        % Don't read the data - just seek forward by the data size and get
        % the size of the next dataframe
        if i==1
            Nt_total = Nt_total + Nt;
        end
        if ftell(fd) + NS*4 < file_size
            fseek(fd,NS*4,'cof');
        elseif ftell(fd) + NS*4 == file_size
            assert(i==DevCount)
            status=false;
            break;
        else
            error('file read error - unexpected size');
        end
    end

    c=c+Nt;
    nframe = nframe+1;

end	% endwhile

Nt_total
nframe

data = zeros([Nt_total,max(NumEnChan),DevCount],'single');

% now we have a known number of aggregated data frames
fseek(fd,data_start,"bof");

status=0; 	% EOF marker
c=1;
Trel=0;
while status==0
    for i=1:DevCount
        TSec=fread(fd,1,'uint64');
        TNsec=fread(fd,1,'uint32');
        NS=fread(fd,1,'uint32');
        Nt=NS/NumEnChan(i);

        d=fread(fd,[NumEnChan(i),Nt],'float32')'; % [number of enabled channels, number of times]
        data(c:c+Nt-1,1:NumEnChan(i),i)=d;
    end

    c=c+Nt;
    status=feof(fd);
end	% endwhile

fclose(fd);
T = zeros([sum(NumEnChan)/2,Nt_total],'single');
P = zeros([sum(NumEnChan)/2,Nt_total],'single');
% convert voltage to temperature, assuming chan 2 and 4 on each device are temp voltage
for i=1:DevCount
    T((i-1)*2+1,:)=data(:,2,i);
    T((i-1)*2+2,:)=data(:,4,i);
    P((i-1)*2+1,:)=data(:,1,i);
    P((i-1)*2+2,:)=data(:,3,i);
end

% apply calibrations
for i=1:DevCount*2
    ind = find(calibration_table(:,1)==header.pressure_sensor_serial_numbers(i));
    if ~isempty(ind)
        T(i,:)=(T(i,:)-calibration_table(ind,3))/calibration_table(ind,4)+25;
        P(i,:)=P(i,:)+calibration_table(ind,2);
    else
        error('sensor not found')
    end
end

end