% script to read sac files containing Old Faithful P/T data from Kedar's experiment in 1994
% and do some preliminary analyses
%
% RAS, 11/15
%
% revised to make figure for Rudolph et al. paper 2024

clear;
close all;
format compact
addpath ..
addpath ../..

% Pressure
p(1)=rdsac('cyc.1.sac');
p(2)=rdsac('cyc.2.sac');
p(3)=rdsac('cyc.3.sac');

% velocity - Station BD00
v(1)=rdsac('bd2.00.a.cyc.1.sac');
v(2)=rdsac('bd2.00.a.cyc.2.sac');
v(3)=rdsac('bd2.00.a.cyc.3.sac');
% velocity - Station BD60
v(4)=rdsac('bd2.60.a.cyc.1.sac');
v(5)=rdsac('bd2.60.a.cyc.2.sac');
v(6)=rdsac('bd2.60.a.cyc.3.sac');
% convert time to datetime object
for i=1:3
    p(i).tdate = datetime(p(i).t,'ConvertFrom','datenum');
    p(i).ts = (0:length(p(i).d)-1)*p(i).HEADER.DELTA;
end
for i=1:6
    v(i).tdate = datetime(v(i).t,'ConvertFrom','datenum');
    v(i).ts = (0:length(v(i).d)-1)*v(i).HEADER.DELTA;
end

%% plot raw data:
figure(101); clf;
subplot(3,1,1);
h1 = gca();
for i=1:2
    plot(p(i).tdate,p(i).d,'DisplayName',['P' num2str(i)]);
    hold on
end
legend()
subplot(3,1,2);
h2 = gca();
for i=1:3
    plot(v(i).tdate,v(i).d);
    hold on
end
subplot(3,1,3);
for i=4:6
    plot(v(i).tdate,v(i).d);
    hold on
end
h3 = gca();
linkaxes([h1,h2,h3],'x')
%% plot raw data (second time - this time use seconds)
figure(102); clf;
subplot(3,1,1);
h1 = gca();
for i=1:2
    plot(p(i).ts,p(i).d,'DisplayName',['P' num2str(i)]);
    hold on
end
legend()
subplot(3,1,2);
h2 = gca();
for i=1:3
    plot(v(i).ts,v(i).d);
    hold on
end
subplot(3,1,3);
for i=4:6
    plot(v(i).ts,v(i).d);
    hold on
end
h3 = gca();
linkaxes([h1,h2,h3],'x');
tstart = 500;
tend = 1500;
set(gca,'XLim',[tstart,tend]);


%% high pass filter the seismic velocity
pass_frequency = 0.1;
fs = 250;
for i=1:6
    vfilt(i).d = highpass(v(i).d,pass_frequency,fs);
    vfilt(i).ts = v(i).ts;
    dt = v(i).HEADER.DELTA;
    fs = 1/dt;
end
%% plot filtered seismic data together with the pressure signal
figure(102); clf;
t=tiledlayout(7,1,'padding','none','TileSpacing','tight');
% Pressure
nexttile(t)
xlim = [1400 1410];
% subplot(2,1,1);
plot(p(1).ts,p(1).d);
hold on
plot(p(2).ts,p(2).d);
h=[];
h(1) = gca();
% Seismics
for i=1:3
    nexttile(t);
    plot(vfilt(i).ts,vfilt(i).d);
    % set(gca,'XLim',xlim);
    h(i+1) = gca();
    
end
% subplot(2,1,2);
for i=4:6
    nexttile(t);    

    plot(vfilt(i).ts,v(i).d);

    % set(gca,'XLim',xlim);
    % plot(p(1).ts,p(1).d);
    h(i+1) = gca();

end
linkaxes(h,'x');
set(gca,'XLim',xlim);
set(gcf,'Color','none')
exportgraphics(t,'pressure_and_seismic.pdf')

%% next... plot a spectrum for part of the pressure time series
% start time for data
t0=p(1).t(1);
startime=datevec(t0)
% end time for data
xl=length(p(1).t);
t1=p(1).t(xl);
endtime=datevec(t1)
% sample rate
dt=(p(1).t(2)-p(1).t(1))*(24*3600);
fs=round(1/dt);
dt=1/fs;

% excerpt data for analysis
% t1a=datenum(1994, 10, 20, 18, 08, 00);
% t2a=datenum(1994, 10, 20, 18, 16, 00);
t1a=datenum(1994, 10, 20, 18, 05, 00);
t2a=datenum(1994, 10, 20, 18, 25, 00);
% t2a=datenum(1994, 10, 20, 18, 19, 20);
ti=find(p(1).t >= t1a & p(1).t <= t2a);
t=p(1).t(ti);
p1a=p(1).d(ti);
p2a=p(2).d(ti);
ts=0:dt:(length(t)-1)*dt;

figure(103); clf;
subplot(2,1,1);
h1=plot(t,p1a,'r');
hold on
h2=plot(t,p2a,'b');
datetick('x','HH:MM')
xlabel('Time, Oct. 20, 1994')
ylabel('Pressure signal')

% estimate spectra
p1u=detrend(p1a,1); % remove linear trend from data
p2u=detrend(p2a,1);
[ee,vv]=dpss(length(p1u),3);
unc=1;
[Pxx,Pyy,Pxy,Txy,Cxy,Exx,Eyy,ETxy,pX,pY,pXY,f,DOFxys,zsl,Ep_gf,EC_gf,r]= mt_cspek_phs(p1u,p2u,fs,ee,vv,unc);
subplot(2,1,2);
plot(f,Pxx);
hold on
plot(f,Pyy);
set(gca,'XLim',[0 20])
set(gca,'YLim',[0 0.3e-7])
xlabel('Frequency (Hz)');
ylabel('Power (db/Hz)?')

%% compute spectrograms for pressure and vertical component broadband
window_size=fs*1*30;
overlap = window_size*0.95;
% t1a and t2a are set in the previous cell.
for i=1:2
    mask = p(i).t >= t1a & p(i).t <= t2a;
    
    data = detrend(p(i).d(mask),1);
    [P2,F2,T2] = mt_spectrogram(data,window_size,overlap,fs);
    p(i).P2 = P2;
    p(i).F2 = F2;
    p(i).T2 = T2;
    p(i).mask = mask;
end
for i = [1 4]
    mask = v(i).t >= t1a & v(i).t <= t2a;
    
    data = detrend(vfilt(i).d(mask),1);
    [P2,F2,T2] = mt_spectrogram(data,window_size,window_size*0.95,fs);
    v(i).P2 = P2;
    v(i).F2 = F2;
    v(i).T2 = T2;
    v(i).tmask = mask;

end

% [S2,F2,T2,P2]=spectrogram(p2u,hamming(window_size),[window_size*0.90],[],fs,'centered');


% apply correction to convert Kedar pressure units to equivalent meters of water
calibrate_p = @(x) x*3.33e07/(1000*9.81)-5;
%%
fh=figure(201)
fh.Position(3:4) = [655 738];
clf;
% top panel - pressure signal
% next panel - pressure spectrogram
% next panel - vv
vchan=1; % which velocity seismogram to plot
% bottom panel - vv spectrogram

% set(gcf,'Position',[560   204   560 360]);
t=tiledlayout(6,1,'TileSpacing','tight');

fontsize=14;
h=[];
nexttile(t);
pchan = 1;
mask = p(pchan).mask;
tplot = p(pchan).ts( mask );
pplot = calibrate_p(p(pchan).d(mask));
plot(tplot,pplot);
h(1) = gca();
ylabel('Conduit liquid level (m)')
title(['Pressure Sensor ' num2str(pchan)]);
hold on;

nexttile(t,[2 1]);
pcolor(p(pchan).T2+tplot(1),p(pchan).F2,10*log10(p(pchan).P2)); shading flat;
% set(gcam)
% h4=pcolor(T2,F2,10*log10(P2));
% pcolor(T2,F2,P2);
shading flat
set(gca,'ColorScale','linear');
set(gca,'CLim',[-150 -60])
set(gca,'YScale','log');

% set(gca,'ylim',[0.1 50]);
set(gca,'YLim',[1e-2 1e2])
h(2) = gca();


ylabel('Frequency (Hz)')
% set(gca,'FontSize',fontsize);
% set(gca,'FontName','Helvetica');
ha2 = gca;
hcb=colorbar;
hcb.Label.String = 'PSD (dB/Hz)';

nexttile(t);
mask = v(vchan).tmask;
tplot = v(vchan).ts(mask);
plot(tplot,v(vchan).d(mask));
ylabel('Velocity');
h(3) = gca();

nexttile(t,[2 1]);
pcolor(v(vchan).T2+tplot(1),v(vchan).F2,10*log10(v(vchan).P2)); shading flat;
set(gca,'ColorScale','linear')
set(gca,'CLim',[-150 -40])
set(gca,'YScale','log')
set(gca,'YLim',[1e-2 1e2])
h(4) = gca();
ha2 = gca;
hcb=colorbar;
hcb.Label.String = 'PSD (dB/Hz)';
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gcf,'Color','white')

exportgraphics(t,['Spectrogram_' num2str(window_size) '.pdf'])%,'ContentType','vector')