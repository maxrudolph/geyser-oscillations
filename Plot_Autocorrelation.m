% script to read binary output files from Kistler LabAmps and plot time series and spectral data
%
% RAS, 4/22
%

clear all
close all
format compact

% set params
% NW=2;	% time-bandwidth product for multi-taper spectral estimation

% viewing order
vo=[6 4 3 2 1 5];

% load data

% file name
% iS=['/Volumes/LoneStar/2023/
% iS=['./05-18/Eruption1_1inch_topconst_20230518_15_28_59.bin'];
% iS = ['./05-18/Eruption4_1inch_topconst_midconst-20230518-18-19-41.bin'];
% iS = ['./05-17/Eruption3_1inch_midconstrict1a_-20230517-16-20-38.bin']
iS = ['./05-18/Eruption1_1inch_topconst-20230518-15-31-15.bin'];

% read binary
fd=fopen(iS,"rb");

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

% unknown number of aggregated data frames
status=0; 	% EOF marker
c=1;
Trel=0;

while status==0

  for i=1:DevCount
    TSec=fread(fd,1,'uint64');
    TNsec=fread(fd,1,'uint32');
    NS=fread(fd,1,'uint32');
    Nt=NS/NumEnChan(i);
    d=fread(fd,[NumEnChan(i),Nt],'float32');
    data(c:c+Nt-1,1:NumEnChan(i),i)=d';
  end

  c=c+Nt;
  status=feof(fd);

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
%% Set min/max time
% define window to view
t1=100;	% start time in seconds, start of record = 0
t2=400;	% end time in seconds, end of record = 0

if t1 == 0
  t1=1;
else
  t1=t1*Fs;
end
if t2 == 0
  t2=length(t);
else
  t2=t2*Fs;
end

%% plot time series

% pressure
figure(1)
clf

for i=1:6
  tS=['a1',num2str(i),'=subplot(6,1,i);'];
  eval(tS)
  yyaxis left
  plot(t(t1:t2),P(vo(i),t1:t2),'b')
  tS=['Sensor ',num2str(vo(i))];
  ylabel('P, bar')
  yyaxis right
  plot(t(t1:t2),T(vo(i),t1:t2),'r')
  ylabel('T, °C')
  title(tS)
  if t1==1
    set(gca,'XLim',[0 t2/Fs]);
  else
    set(gca,'XLim',[t1/Fs t2/Fs]);
  end
end

%linkaxes([a11 a12],'x')
linkaxes([a11 a12 a13 a14 a15 a16],'x')

% temperature
figure(2)
clf
orient tall

for i=1:6
  tS=['a2',num2str(i),'=subplot(6,1,i);'];
  eval(tS)
  plot(t,T(i,:),'r')
  tS=['Sensor ',num2str(i)];
  ylabel('T, C')
  title(tS)
end

%% Play with computing an autocorrelation function
signal = detrend(P(3,:));
% filter signal
myfilt=designfilt('lowpassfir','PassbandFrequency',10,'StopbandFrequency',20,'PassbandRipple',1,'StopbandAttenuation',60,'SampleRate',10000);
signal_filt = filtfilt(myfilt,signal);

% signal_filt = lowpass(signal,10,Fs);
% window length, in samples
window_length = 5.0*Fs;
shift = window_length;
N = length(signal_filt);
clear tstart ind
istart = 1:shift:(N-window_length)
max_lag = 4*Fs; % set the maximum lag (in samples)
iout=1;
for ind = istart
    signal1 = detrend(signal_filt(ind:ind+window_length));
    [C,lags] = xcorr(signal1,max_lag,'normalized');
    
    if iout ==1
        Cout = zeros(length(lags(lags>0)),length(istart));
    end
        Cout(:,iout) = C(lags>0);
        iout = iout+1;
end
figure;
ax1=subplot(2,1,1)
plot(t,signal)
hold on
plot(t,signal_filt)
ax1 = gca();
ax2=subplot(2,1,2)
pcolor(t(istart)+window_length/Fs/2,lags(lags>0)/Fs,Cout);
shading flat
colorbar();
ax2 = gca();
ax1.Position([1 3]) = ax2.Position([1 3]);
linkaxes([ax1 ax2],'x');

