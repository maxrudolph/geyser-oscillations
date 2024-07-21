% script to read sac files containing Old Faithful P/T data from Kedar's experiment in 1994
% and do some preliminary analyses
%
% RAS, 11/15
%
% revised to make figure for Sohn and Rudolph paper, 5/16

clear;
close all;
format compact

p1=rdsac('cyc.1.sac');
p2=rdsac('cyc.2.sac');
p3=rdsac('cyc.3.sac');

% start time for data
t0=p1.t(1);
startime=datevec(t0)
% end time for data
xl=length(p1.t);
t1=p1.t(xl);
endtime=datevec(t1)
% sample rate
dt=(p1.t(2)-p1.t(1))*(24*3600);
fs=round(1/dt);
dt=1/fs;

% excerpt data for analysis
t1a=datenum(1994, 10, 20, 18, 00, 40);
t2a=datenum(1994, 10, 20, 18, 20, 20);
% t2a=datenum(1994, 10, 20, 18, 19, 20);
ti=find(p1.t >= t1a & p1.t <= t2a);
t=p1.t(ti);
p1a=p1.d(ti);
p2a=p2.d(ti);
ts=0:dt:(length(t)-1)*dt;

figure(1)
clf
h1=plot(t,p1a,'r');
hold on
h2=plot(t,p2a,'b');
datetick('x','HH:MM')
xlabel('Time, Oct. 20, 1994')
ylabel('Pressure')

% estimate spectra
p1u=detrend(p1a,0);
p2u=detrend(p2a,0);
[e,v]=dpss(length(p1u),3);
unc=1;
%[Pxx,Pyy,Pxy,Txy,Cxy,Exx,Eyy,ETxy,pX,pY,pXY,f,DOFxys,zsl,Ep_gf,EC_gf,r]= mt_cspek_phs(p1u,p2u,fs,e,v,unc);

%% compute spectrogram
window_size=fs*1*240;
% [S2,F2,T2,P2]=spectrogram(p2u,hamming(window_size),[window_size*0.90],[],fs,'centered');
[P2,F2,T2] = mt_spectrogram(p2u,window_size,window_size*0.95,fs);

% apply correction to convert Kedar pressure units to equivalent meters of water
emw2=p2a*3.33e07/(1000*9.81)-5;
%%
figure(2)
clf;
set(gcf,'Position',[560   204   560 360]);
subplot(2,1,1);
fontsize=14;

tz1 = 137;
tz2 = 147;

% clf
% orient portrait
% subplot(211)
h3=plot(ts,emw2,'k');
set(gca,'xlim',[50 max(T2)]);
%set(gca,'ylim',[-1e-03 2e-03]);
set(gca,'ylim',[5 20]);
set(gca,'FontSize',fontsize);
set(gca,'FontName','Helvetica');
ylabel('Conduit liquid level (m)')
title('Kedar et al. "Middle" Sensor');
hold on;
% plot(tz1*[1 1],get(gca,'ylim'),'r:');
% plot(tz2*[1 1],get(gca,'ylim'),'r:');
yl = get(gca,'ylim');
x = [tz1 tz1 tz2 tz2];
y = [yl(2) yl(1) yl(1) yl(2)];
% hp=patch(x,y,0.8*[1 1 1]);
% hp.EdgeColor=hp.FaceColor;
set(gca,'children',flipud(get(gca,'children')))
ha1 = gca;
text(0.01,0.9,'A','FontSize',16,'Units','normalized');
subplot(2,1,2)
% h4=pcolor(T2,F2,10*log10(P2));
pcolor(T2,F2,P2);
shading flat
set(gca,'ColorScale','log');
set(gca,'YScale','log');
% shading('interp')
% cmin=10*log10(min(min(P2)));
% cmax=10*log10(max(max(P2)));
% crange=cmax-cmin
% cmax2=cmin+0.8*crange;
% cmin2=cmax-0.35*crange;
% set(gca,'Clim',[cmin2 cmax2]);
set(gca,'xlim',[50 max(T2)]);
set(gca,'ydir','normal');
set(gca,'ylim',10.^[-1. 1.]);
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'FontSize',fontsize);
set(gca,'FontName','Helvetica');
ha2 = gca;
hcb=colorbar;
hcb.Label.String = 'PSD (dB/Hz)';
tzoom=find(ts > 137 & ts <= 147);
pzoom=emw2(tzoom);
drawnow;
ha1.Position(3) = ha2.Position(3)
text(0.01,0.9,'B','FontSize',16,'Units','normalized');
% subplot(3,1,1);
% % clf
% h5=plot(ts(tzoom),pzoom,'k')
% xlabel('Time (s)')
% ylabel('Conduit liquid level (m)')
% set(gca,'FontSize',fs);
% set(gca,'FontName','Helvetica');
% ha3=gca;
% ha3.Position(3) = ha2.Position(3)
% set(gcf,'Color','w');
% text(0.01,0.9,'A','FontSize',16,'Units','normalized');
exportgraphics(gcf,'figure4_proposal_version.pdf');
% % export_fig('figure4.png');

