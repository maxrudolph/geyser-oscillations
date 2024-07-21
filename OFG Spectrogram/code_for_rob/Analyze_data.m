% script to read sac files containing Old Faithful P/T data from Kedar's experiment in 1994
% and do some preliminary analyses
%
% RAS, 11/15

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
t1a=datenum(1994, 10, 20, 17, 59, 0);
t2a=datenum(1994, 10, 20, 18, 21, 0);
ti=find(p1.t >= t1a & p1.t <= t2a);
t=p1.t(ti);
p1a=p1.d(ti);
p2a=p2.d(ti);
p3a=p3.d(ti);
ts=0:dt:(length(t)-1)*dt;

figure(1)
clf
h1=plot(t,p1a,'r');
hold on
h2=plot(t,p2a,'b');
h3=plot(t,p3a,'g');
datetick('x','HH:MM')
xlabel('Time, Oct. 20, 1994')
ylabel('Pressure')

% estimate spectra
p1u=detrend(p1a,0);
p2u=detrend(p2a,0);
p3u=detrend(p3a,0);
[e,v]=dpss(length(p1u),3);
unc=1;
% [Pxx,Pyy,Pxy,Txy,Cxy,Exx,Eyy,ETxy,pX,pY,pXY,f,DOFxys,zsl,Ep_gf,EC_gf,r]= mt_cspek_phs(p1u,p2u,fs,e,v,unc);

% compute spectrogram
[S1,F1,T1,P1]=spectrogram(p1u,hamming(fs*3*60),[fs*3*60*0.95],[],fs);
[S2,F2,T2,P2]=spectrogram(p2u,hamming(fs*3*60),[fs*3*60*0.95],[],fs);
[S3,F3,T3,P3]=spectrogram(p3u,hamming(fs*3*60),[fs*3*60*0.95],[],fs);

%% Make first figure
Splot = S2;
Fplot = F2;
Tplot = T2;
Pplot = P2;

figure(99);
h=gcf;
h.Position(3:4) = [576   800];
subplot(3,1,3);
fmask = (Fplot <= 10) & (Fplot >= 0.1);
% contourf(T1,F1(fmask),log10(P1(fmask,:)),200,'Color','none');
pcolor(Tplot,Fplot(fmask),10*log10(Pplot(fmask,:)));
hcb=colorbar;
hcb.Label.String = 'PSD (db/Hz)';
shading flat
set(gca,'YScale','log');
set(gca,'YLim',[0.1 10]);
ha1 = gca;
ha1.FontSize=14;
ha1.FontName='Helvetica';
set(gca,'Layer','top');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

subplot(3,1,1);
plot((t-t(1))*3600*24,p2u);
ha2=gca;
ha2.XLim = ha1.XLim;
ha2.Position(3) = ha1.Position(3);
ha2.FontSize=14;
ha2.FontName='Helvetica';

subplot(3,1,2);
tplot = (t-t(1))*3600*24;
mask = tplot>1000 & tplot<1010;
plot(tplot(mask),p2u(mask));
ha3=gca;
ha3.Position(3) = ha1.Position(3);
ha3.FontSize=14;
ha3.FontName='Helvetica';

%%

% % plot
% figure(2)
% clf
% h1=semilogx(f,10*log10(Pxx),'r');
% hold on
% h2=semilogx(f,10*log10(Pyy),'b');
% xlabel('Frequency, Hz')
% ylabel('Power')
% 
% figure(3)
% clf
% subplot(211)
% semilogx(f,Cxy,'k')
% xlabel('Frequency, Hz')
% ylabel('Coherency, gamma-squared')
% set(gca,'Xlim',[1e-04 1e02]);
% subplot(212)
% semilogx(f,pXY,'k')
% xlabel('Frequency, Hz')
% ylabel('Phase, rad')
% set(gca,'Xlim',[1e-04 1e02]);
% 
% %%
% figure(4)
% clf
% orient tall
% subplot(611)
% h1=plot(ts,detrend(p1a,0),'r');
% set(gca,'xlim',[0 max(T1)]);
% set(gca,'ylim',[-1e-03 2e-03]);
% ylabel('Pressure')
% title('Top Sensor')
% subplot(6,1,2:3)
% % h2=pcolor(T1,log10(F1),10*log10(P1));
% contourf(T1,log10(F1(2:end)),10*log10(P1(2:end,:)),256,'Color','none');
% % shading('interp')
% cmin=10*log10(min(min(P1)));
% cmax=10*log10(max(max(P1)));
% crange=cmax-cmin
% cmax2=cmin+0.85*crange;
% cmin2=cmax-0.45*crange;
% set(gca,'Clim',[cmin2 cmax2]);
% set(gca,'xlim',[0 max(T1)]);
% set(gca,'ydir','normal');
% set(gca,'ylim',[log10(min(F1)) 1.]);
% ylabel('Log10 frequency, Hz')
% subplot(614)
% h3=plot(ts,detrend(p2a,0),'b');
% set(gca,'xlim',[0 max(T2)]);
% set(gca,'ylim',[-1e-03 2e-03]);
% ylabel('Pressure')
% title('Bottom Sensor')
% subplot(6,1,5:6)
% % h4=pcolor(T2,log10(F2),10*log10(P2));
% contourf(T2,log10(F2(2:end)),10*log10(P2(2:end,:)),256,'Color','none');
% shading('interp')
% cmin=10*log10(min(min(P2)));
% cmax=10*log10(max(max(P2)));
% crange=cmax-cmin
% cmax2=cmin+0.85*crange;
% cmin2=cmax-0.45*crange;
% set(gca,'Clim',[cmin2 cmax2]);
% set(gca,'xlim',[0 max(T2)]);
% set(gca,'ydir','normal');
% set(gca,'ylim',[log10(min(F2)) 1.]);
% xlabel('Time, seconds')
% ylabel('Log 10 frequency, Hz')
% 
% [b,a]=butter(4,.05/125);
% p1f=filtfilt(b,a,detrend(p1a,0));
% p2f=filtfilt(b,a,detrend(p2a,0));
% %%
% figure(5)
% clf
% plot(ts,p1f,'r')
% hold on
% plot(ts,p2f,'b')
