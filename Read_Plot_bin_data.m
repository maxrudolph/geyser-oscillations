% script to read binary output files from Kistler LabAmps and plot time series and spectral data
%
% RAS, 4/22
%

clear all
format compact

% set params
NW=2;	% time-bandwidth product for multi-taper spectral estimation

% define window to view
t1=0;	% start time in seconds, start of record = 0
t2=0;	% end time in seconds, end of record = 0

% viewing order
vo=[4 6 3 2 1 5];

% load data

% file name
% iS=['/Volumes/LoneStar/2023/Eruption1_1inch_topconst_20230518_15_28_59.bin'];
iS = ['./05-18/Eruption4_1inch_topconst_midconst-20230518-18-19-41.bin'];

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

%linkaxes([a21 a22],'x')
linkaxes([a21 a22 a23 a24 a25 a26],'x')

% estimate spectra
for i=1:6
  [Pxx(i,:),f]=pmtm(P(i,t1:t2)-mean(P(i,t1:t2)),NW,length(P(i,t1:t2)),Fs);
%  [Txx(i,:),f]=pmtm(T(i,:)-mean(T(i,:)),NW,length(T(i,:)),Fs);
end

% plot spectra
figure(3)
clf

for i=1:6
  subplot(2,3,i)
  semilogx(f,10*log10(Pxx(vo(i),:)),'b')
  grid
  xlabel('Frequency, Hz')
  ylabel('Power, dB')
  tS=['Sensor ',num2str(vo(i))];
  title(tS)
end

%figure(4)
%clf

%for i=1:6
%  subplot(1,6,i)
%  semilogx(f,10*log10(Txx(i,:)),'r')
%  grid
%  xlabel('Frequency, Hz')
%  ylabel('Power, dB')
%  tS=['Sensor ',num2str(i)];
%  title(tS)
%end

%figure(5)
%clf
%% 
for i=4:4
  figure(i+3)
  clf
  
  [S,F,TT,PP]=spectrogram(P(vo(i),t1:t2),5*Fs,0,[],Fs,'yaxis');
  subplot(7,1,1)
  plot(t(t1:t2),P(vo(i),t1:t2),'k')
  ylabel('P, bar')
  set(gca,'XLim',[t1/Fs t2/Fs]);
  tS=['Sensor ',num2str(vo(i))];
  title(tS)
  subplot(7,1,2:4)
  pcolor(TT,log10(F),10*log10(PP));%,'edgecolor','none')
  view(0,90)
  ii=find(F >= 10 & F < 1000);
  cmin=10*log10(min(min(PP(ii))));
  cmax=10*log10(max(max(PP(ii))));
  crange=cmax-cmin;
  cmax2=cmin+1.0*crange;
  cmin2=cmax-0.5*crange;
  set(gca,'Clim',[cmin2 cmax2]);
  set(gca,'YLim',[log10(10) log10(1000)]);
  set(gca,'XLim',[0 (t2-t1)/Fs]);
  shading flat
%   [S,F,TT,PP]=spectrogram(P(vo(i),t1:t2),hamming(25*Fs),0,[],Fs);
%   subplot(7,1,5:7)
%   surf(TT,log10(F),10*log10(PP),'edgecolor','none')
%   view(0,90)
%   ii=find(F < 10 & F > 0.1);
%   cmin=10*log10(min(min(PP(ii))));
%   cmax=10*log10(max(max(PP(ii))));
%   crange=cmax-cmin;
%   cmax2=cmin+1.0*crange;
%   cmin2=cmax-0.5*crange;
%   set(gca,'Clim',[cmin2 cmax2]);
%   set(gca,'YLim',[log10(0.1) log10(10)]);
%   set(gca,'XLim',[0 (t2-t1)/Fs]);
%   shading flat
end
