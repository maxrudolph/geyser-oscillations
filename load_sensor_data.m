function load_sensor_data(filename)
% read binary
fd=fopen(filename,"rb");

% header
FileVersion=fread(fd,1,'int32');
Fs=fread(fd,1,'float');
DevCount=fread(fd,1,'uint32');

for i=1:DevCount
    DevID(i)=fread(fd,1,'int32');
    SNL(i)=fread(fd,1,'uint32');
    SN(i,:)=fread(fd,SNL(i),'char');% read as bytes
    NameL(i)=fread(fd,1,'uint32');
    Name(i,:)=fread(fd,NameL(i),'*char');
    NumEnChan(i)=fread(fd,1,'uint32');

    for j=1:NumEnChan(i)
        ChanNum(i,j)=fread(fd,1,'int32');
    end

end
SN = char(SN);%convert to character string
% Get data acquisition system box order:
order_str = Name(:,end);
% order_name = np.argsort(Name[:, -1])
order_num = str2num(order_str);
[~,sort_order] = sort(order_num);
num_sens = 2*DevCount;
order = zeros(num_sens,1); % order = np.empty(self.num_sens, 'int')
for i = 1:length(order_num)
    order(2*(i-1)+1) = 2*(sort_order(i)-1)+1;
    order(2*(i-1)+2) = 2*(sort_order(i)-1)+2;
end
% for i in range(len(order_name)):
% order[2*i] = 2*order_name[i]
% order[2*i+1] = 2*order_name[i]+1
% print("order_name=",order_name)
disp(['order=' , num2str(order)])

% SN = SN[order_name, :]
% Name = Name[order_name, :]




% unknown number of aggregated data frames
c=1;
Trel=0;

while ~feof(fd)
    for i=1:DevCount
        TSec=fread(fd,1,'uint64');
        TNsec=fread(fd,1,'uint32');
        NS=fread(fd,1,'uint32');
        Nt=NS/NumEnChan(i);
        d=fread(fd,[NumEnChan(i),Nt],'float32');
        data(c:c+Nt-1,1:NumEnChan(i),i)=d';
    end
    c=c+Nt;    
end	% endwhile

fclose(fd);

% convert voltage to temperature, assuming chan 2 and 4 on each device are temp voltage
for i=1:DevCount
    T((i-1)*2+1,:)=(data(:,2,i)-1.478)/0.01746+25;
    T((i-1)*2+2,:)=(data(:,4,i)-1.478)/0.01746+25;
    P((i-1)*2+1,:)=data(:,1,i);
    P((i-1)*2+2,:)=data(:,3,i);
end

t=0:1/Fs:(length(T(1,:))-1)/Fs; % relative time vector

% apply pressure calibrations

% Array with sensor SN with index corresponding to position
sens_used = [5122778, 5122769, 5940428, 5122770, 5122777, 5940430];% serial numbers of sensors 1-6


