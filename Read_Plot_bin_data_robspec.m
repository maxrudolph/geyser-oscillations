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
% iS = ['./05-18/Eruption1_1inch_topconst-20230518-15-31-15.bin'];


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
%% try to use multitaper code to plot phase between signals
signal1 = P(6,t1:t2);
signal2 = P(4,t1:t2);
[e,v] = dpss(length(signal1),2);
[Pxx,Pyy,Pxy,Txy,Cxy,Exx,Eyy,ETxy,pX,pY,pXY,f,DOFxys,zsl,Er_p,Er_c,Er_Tf,r] = mt_cspek_phs(signal1,signal2,Fs,e,v,1);
%% 
figure();
ax1=subplot(3,1,1);
plot(f,pXY);
set(gca,'XScale','log')
ax2=subplot(3,1,2);
plot(f,Cxy)
hold on
plot(f,zsl,'r')
set(gca,'XScale','log')
ax3=subplot(3,1,3)
v = 1.05*2*pi*(f./pXY);
plot(f,v);
set(gca,'XScale','log');
set(gca,'Ylim',[-20 20]);

linkaxes([ax1,ax2,ax3],'x')


%% try Rob's multitaper code to make a spectrogram

% [Pxx,Exx,pX,f]=mt_cspek_phs(P(2,t1:t2),Fs,1);
% figure()
% plot(f,Pxx);
window_length = 50*Fs;
overlap = window_length*0.9;
pp = P(2,t1:t2);
N1 = length(pp);

ind=1;
start_indices = 1:(window_length-overlap):(N1-window_length-1);
start = start_indices(ind);
signal = pp(start:(start+window_length-1));
%window = hamming(window_length)';
window = ones(size(signal));

[Pxx,Exx,pX,frequency]=mt_phs(signal.*window,Fs,1);
all_Pxx = zeros(length(Pxx),length(start_indices));
all_Pxx(:,ind) = Pxx;

parfor ind = 2:length(start_indices)
    start = start_indices(ind);
    signal = pp(start:(start+window_length-1));
    [Pxx,Exx,pX,frequency]=mt_phs(window.*signal,Fs,1);
    if ind==1
       all_Pxx(:,ind) = Pxx;
    else
       all_Pxx(:,ind) = Pxx; 
    end
end


%% spectrogram figure
% figure('Position',[535 210 2158 701]);
figure();
fig=gcf();
fig.Position(3) = fig.Position(4)*2;
subplot(2,1,1);
plot(t(t1:t2),pp);
ax1 = gca();
subplot(2,1,2);
pcolor(t(t1)+(start_indices+window_length/2-1)/Fs,frequency,all_Pxx); shading flat; 
% hold on
% yline(0.28,'k--')
ax2 = gca();
set(ax2,'XLim',ax1.XLim)
set(gca,'ColorScale','log'); colorbar; drawnow;
set(ax2,'YScale','log');
set(gca,'YLim',[max(1/window_length/Fs,min(frequency)) 100]);
linkaxes([ax1,ax2],'x')
ax1.Position(3) = ax2.Position(3);
yline(1/7,'r--')


%% very slow power spectrum
%linkaxes([a21 a22],'x')
% linkaxes([a21 a22 a23 a24 a25 a26],'x')

% estimate spectra
% for i=1:6
%   [Pxx(i,:),f]=pmtm(P(i,t1:t2)-mean(P(i,t1:t2)),NW,length(P(i,t1:t2)),Fs);
%  [Txx(i,:),f]=pmtm(T(i,:)-mean(T(i,:)),NW,length(T(i,:)),Fs);
% end

% plot spectra
% figure(3)
% clf

% for i=1:6
%   subplot(2,3,i)
%   semilogx(f,10*log10(Pxx(vo(i),:)),'b')
%   grid
%   xlabel('Frequency, Hz')
%   ylabel('Power, dB')
%   tS=['Sensor ',num2str(vo(i))];
%   title(tS)
% end

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
% for i=4:4
%   figure(i+3)
%   clf
% 
%   [S,F,TT,PP]=spectrogram(P(vo(i),t1:t2),5*Fs,0,[],Fs,'yaxis');
%   subplot(7,1,1)
%   plot(t(t1:t2),P(vo(i),t1:t2),'k')
%   ylabel('P, bar')
%   set(gca,'XLim',[t1/Fs t2/Fs]);
%   tS=['Sensor ',num2str(vo(i))];
%   title(tS)
%   subplot(7,1,2:4)
%   pcolor(TT,log10(F),10*log10(PP));%,'edgecolor','none')
%   view(0,90)
%   ii=find(F >= 10 & F < 1000);
%   cmin=10*log10(min(min(PP(ii))));
%   cmax=10*log10(max(max(PP(ii))));
%   crange=cmax-cmin;
%   cmax2=cmin+1.0*crange;
%   cmin2=cmax-0.5*crange;
%   set(gca,'Clim',[cmin2 cmax2]);
%   set(gca,'YLim',[log10(10) log10(1000)]);
%   set(gca,'XLim',[0 (t2-t1)/Fs]);
%   shading flat
% %   [S,F,TT,PP]=spectrogram(P(vo(i),t1:t2),hamming(25*Fs),0,[],Fs);
% %   subplot(7,1,5:7)
% %   surf(TT,log10(F),10*log10(PP),'edgecolor','none')
% %   view(0,90)
% %   ii=find(F < 10 & F > 0.1);
% %   cmin=10*log10(min(min(PP(ii))));
% %   cmax=10*log10(max(max(PP(ii))));
% %   crange=cmax-cmin;
% %   cmax2=cmin+1.0*crange;
% %   cmin2=cmax-0.5*crange;
% %   set(gca,'Clim',[cmin2 cmax2]);
% %   set(gca,'YLim',[log10(0.1) log10(10)]);
% %   set(gca,'XLim',[0 (t2-t1)/Fs]);
% %   shading flat
% end
