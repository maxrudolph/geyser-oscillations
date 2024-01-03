% script to read binary output files from Kistler LabAmps and examine pre-eruption wobbles as a function
% of water level
%
% RAS, 12/23
%
clear;
close all;

% define window to view
t1=0; % start time in seconds, start of record = 0
t2=190; % end time in seconds, end of record = 0
% viewing order
vo=[6 4 3 2 1 5];

inst=4; % sensor # for spectragram

% is data already loaded?
exist P ;

if ans ~= 1

% load data
% file name
  % iS=['/Volumes/LoneStar/2023/Eruption2b_1inch_topconst_20230518_16_21_51.bin'];
  iS=['./05-18/Eruption2b_1inch_topconst-20230518-16-21-51.bin'];
  pamb=1.032;	% ambient pressure, bar


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

% assign SN to sensors
    if SN(i,:)=='5954878'
      SSN((i-1)*2+1)=5122778;
      SSN((i-1)*2+2)=5122769;
    elseif SN(i,:)=='5954876'
      SSN((i-1)*2+1)=5940428;
      SSN((i-1)*2+2)=5122770;
    elseif SN(i,:)=='5954877'
      SSN((i-1)*2+1)=5122777;
      SSN((i-1)*2+2)=5940430;
    end

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
    T((i-1)*2+1,:)=data(:,2,i);
    T((i-1)*2+2,:)=data(:,4,i);
    P((i-1)*2+1,:)=data(:,1,i);
    P((i-1)*2+2,:)=data(:,3,i);
  end

% apply calibrations
  for i=1:DevCount*2
    if SSN(i)==5122778
      T(i,:)=(T(i,:)-1.475)/0.01748+25;
      P(i,:)=P(i,:)-0.22586478;
    elseif SSN(i)==5122770
      T(i,:)=(T(i,:)-1.470)/0.01742+25;
      P(i,:)=P(i,:)-0.2425978;
    elseif SSN(i)==5122769
      T(i,:)=(T(i,:)-1.465)/0.01738+25;
      P(i,:)=P(i,:)-0.24022625;
    elseif SSN(i)==5940428
      T(i,:)=(T(i,:)-1.477)/0.01744+25;
      P(i,:)=P(i,:)-0.2364321;
    elseif SSN(i)==5122777
      T(i,:)=(T(i,:)-1.480)/0.01739+25;
      P(i,:)=P(i,:)-0.22593766;
    elseif SSN(i)==5940430
      T(i,:)=(T(i,:)-1.484)/0.01734+25;
      P(i,:)=P(i,:)-0.25612138;
    end

  end

else
  clear Pd Pxx f txx td watlev 

end	% end if data exists

t=0:1/Fs:(length(T(1,:))-1)/Fs; % relative time vector

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

% plot raw data first
figure(1)
clf
orient tall
for i=1:6
  tS=['a1',num2str(i),'=subplot(6,1,i);'];
  eval(tS)
  plot(t(t1:t2),P(vo(i),t1:t2),'b')
  tS=['Sensor ',num2str(vo(i))];
  ylabel('P, bar')
  if i==1
    tS=[iS,' ',tS];
  end
  title(tS)
  if t1==1
    set(gca,'XLim',[0 t2/Fs]);
  else
    set(gca,'XLim',[t1/Fs t2/Fs]);
  end
  grid
end
xlabel('Time, s')
linkaxes([a11 a12 a13 a14 a15 a16],'x')
drawnow

% extract segment and decimate
dec=2000;
fs=Fs/dec;
td=decimate(t(t1:t2),dec);

for i=1:6
  x=detrend(P(i,t1:t2));
  Pd(i,:)=decimate(x,dec);
end

figure(2)
clf
orient tall
for i=1:6
  tS=['b1',num2str(i),'=subplot(6,1,i);'];
  eval(tS)
  plot(td,Pd(vo(i),:),'r')
  grid
  tS=['Sensor ',num2str(vo(i))];
  title(tS)
  ylabel('P, bar')
  if t1==1
    set(gca,'XLim',[0 t2/Fs]);
  else
    set(gca,'XLim',[t1/Fs t2/Fs]);
  end
end

xlabel('Time, s')
linkaxes([b11 b12 b13 b14 b15 b16],'x')

% get spectra for discrete time windows
winlens=30;	% window length, sec
winlen=winlens*fs;
[j,Nwin]=size(Pd);
Nwin=floor(Nwin/winlen);
nfft=512;
Pxx=zeros(Nwin,nfft/2+1);
txx=0:winlens:(Nwin-1)*winlens;

for n=1:Nwin
  x=Pd(inst,(1+(n-1)*winlen:n*winlen));
  [Pxx(n,:),f]=pmtm(x-mean(x),2,nfft,fs);
  %Pxx(n,:)=Pxx(n,:)/max(Pxx(n,:));
  y=P(3,t1+(n-1)*winlens*Fs:t1+n*winlens*Fs);
  watlev(n)=80+(mean(y')-pamb)/9.81e-04;	% water level, cm
end 

f=f';
Pxx=Pxx';

%% plot spectragram
figure(3)
clf
ax1=subplot(4,1,1:3);
surf(txx,log10(f),10*log10(Pxx),'edgecolor','none');
view(0,90);
set(gca,'YLim',[-0.7 0.0]);
cmin=10*log10(min(min(Pxx)));
cmax=10*log10(max(max(Pxx)));
crange=cmax-cmin;
cmax2=cmin+0.80*crange;
cmin2=cmax-0.6*crange;
set(gca,'CLim',[cmin2 cmax2]);
%set(gca,'YTickMode','manual');
%set(gca,'YTick',[10^-0.7 10^-0.6 10^-0.5 10^-0.4 10^-0.3 10^-0.2 10^-0.1 10^0])
set(gca,'YTickLabel',[{'0.20','0.25','0.32','0.40','0.50','0.63','0.79','1.0'}]);
ylabel('Frequency, Hz')
tS=['Power spectra (dB), Sensor ',num2str(inst)];
title(tS)

subplot(4,1,4)
plot(txx,watlev)
xlabel('Time, s')
ylabel('Wat Lev, cm')
grid
set(gca,'YLim',[100 260]);

%% set up options structure for lab geyser
geometry = 1;
conduit_diameter = 1; % 1 -> 1 inch, 2 -> 2 inches
switch conduit_diameter
    case 1
        par.sc = pi*(1/2*2.54/100)^2; %cross-section of 1-inch tube in cm^2
        par.H = 31.8e-2; % "new" short 1 inch pickup in m        
        %par.H = 58e-2; % old, longer 1 inch pickup in m
    case 2
        par.sc = pi*(2/2*2.54/100)^2; %cross-section of 2-inch tube in cm^2
        par.H = (53-1.8)*1e-2; %Bubble trap height in m
end

switch geometry
    case 1
        fprintf('Using new lab dimensions. \n')
        par.sb = 1852e-4;
        par.sl = 1; %cross-section of lateral connector in m^2
        par.L = 0; % length of lateral connector
end
% sensor 6 - top of conduit
% sensor 4 - middle of conduit
% sensor 3 - bottom of conduit
% sensor 2 - top of tank
% sensor 1 - bottom of tank
t1p = 70; t2p = 100; % begin and end times for analysis
mask = t>=t1p & t <= t2p;

sensor_heights = [0.0 0.0 0.0 0.0 0.0 0.0];% vector of sensor heights, measured up from bottom of tank.

% tank height
tank_height = 0.76;% internal height of tank, in m.

% plot Rob's water level and the pressure from P1
watlev4 = (P(4,mask)-pamb)*1e5/9.81/1e3; % in m above bottom of tank
watlev3 = (P(3,mask)-pamb)*1e5/9.81/1e3; % in m above bottom of tank
watlev1 = (P(1,mask)-pamb)*1e5/9.81/1e3; % in m above bottom of tank
%p2 = p1 - rho_w*g*(water level above bottom of tank) - rho_v*g*(depth of
%steam)
tank_level = -(P(2,:) - P(1,:))*1e5/1000/9.81;


figure(); 
subplot(2,1,1)
plot(t(mask),watlev1)
hold on
plot(t(mask),watlev3)
plot(t(mask),watlev4)
legend('1','3','4')
subplot(2,1,2);
plot(t(mask),tank_level(mask))
delta_4_1 = mean(watlev4-watlev1)
delta_3_1 = mean(watlev3-watlev1)
delta_4_3 = mean(watlev4-watlev3)


