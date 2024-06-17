clear
% close all

% viewing order - sort sensors top to bottom in lab apparatus.
vo=[6 4 3 2 1 5];

% load data
% note - this eruption has a nice pre-eruption period from 150-200 s.
% iS = ['./05-18/Eruption1_1inch_topconst-20230518-15-31-15.bin'];
% start and end of window (in seconds)
% t1 = 150; t2 = 200;
% t1 = 350; t2 = 500; % during eruption ... This shows peak of 0.29 Hz
% t1 = 210; t2 = 235; %pre-eruptive
% t1=237; t2=500; % eruption proper

iS = ['05-19/Eruption1_1inch_topconst-20230519-15-34-50.bin'];
t1=0; t2=250; %filling stage
% t1=100; t2=200; % also filling stage.
% t1=310; t2=350; % system nearly full, pre-eruption

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


if ~exist('t1') || t1 == 0
  t1=1;
else
  t1=t1*Fs;
end
if ~exist('t2') || t2 == 0
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

%% Examine signal of sensor 2 in isolation
% % design a finite impulse response low pass filter
signal = detrend(P(2,t1:t2));
% myfilt=designfilt('lowpassfir','PassbandFrequency',2.0,'StopbandFrequency',4.0,'PassbandRipple',1,'StopbandAttenuation',60,'SampleRate',10000);
% save('lowpass_2-4Hz.mat',"myfilt")
load('lowpass_2-4Hz.mat');
signal_filt = filtfilt(myfilt,signal);
figure, 
plot(t(t1:t2),signal,'k');
hold on
plot(t(t1:t2),signal_filt,'r')

%% Spectrum of the filtered signal
[Pxx,Exx,pX,frequency]=mt_phs(signal,Fs,1);
figure, plot(frequency, Pxx)
set(gca,'XLim',[0 4])
legend([num2str(t1) '--' num2str(t2)])
