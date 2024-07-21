% script to read binary output files from Kistler LabAmps and examine pre-eruption wobbles as a function
% of water level
%
% RAS, 12/23
% MLR, 1/24
%
clear;
close all;

% define window to view
t1=0; % start time in seconds, start of record = 0
t2=0; % end time in seconds, end of record = 0
% viewing order
% vo=[6 4 3 2 1 5];
vo = [6 4 3 2 1 5];

% whether to do calibration?
do_calibration = false;
if do_calibration
    % example of what the calibration table is expected to look like:
    calibration_table = [
       % SSN    % P offset (bar)   T_offset     T_slope
       5122778  0.22586478
       5122770  0.2425978
       5122769  0.24022625
       5940428  0.2364321
       5122777  0.22593766
       5940430  0.25612138
    ];

    calibration_table(:,2) = 0;
else
    load('calibration_05222024.mat');        
end

inst=1; % sensor # for spectrogram

% is data already loaded?
exist P ;

if ans ~= 1
    
    pamb = 1.0110777892; % calibration day (5/22) 29.85711 in Hg
    % iS = ['05-22-2024-calibration/sensor_test_empty_tank-20240522-17-11-00.bin']
    iS = ['05-22-2024-calibration/sensor_test_full_tank-20240522-17-52-00.bin']

    % pamb=1.01364772;	% ambient pressure, bar (5/30)
    % iS = ['05-30-2024/HeatingUp_fulltank-20240530-14-38-13.bin']    
    % iS = ['05-30-2024/SteamInupt-20240530-16-28-33.bin']

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
            SSN((i-1)*2+1)=5940430;% P1?
            SSN((i-1)*2+2)=5940428;% P2?
        elseif SN(i,:)=='5954876'
            SSN((i-1)*2+1)=5122770;% P3?
            SSN((i-1)*2+2)=5122769	;%5122777;% P4?
            % Data Logger '5954877' not recording
            % elseif SN(i,:)=='5954877'
            %   SSN((i-1)*2+1)=5122777;
            %   SSN((i-1)*2+2)=5940430;
        else
            error('no such box')           
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
        ind = find(calibration_table(:,1)==SSN(i));
        if ~isempty(ind)
            % T(i,:)=(T(i,:)-1.475)/0.01748+25;
            P(i,:)=P(i,:)-calibration_table(ind,2);
        else
            error('sensor not found')
        end

        if SSN(i)==5122778
            T(i,:)=(T(i,:)-1.475)/0.01748+25;
            % P(i,:)=P(i,:)-0.22586478;
        elseif SSN(i)==5122770
            T(i,:)=(T(i,:)-1.470)/0.01742+25;
            % P(i,:)=P(i,:)-0.2425978;
        elseif SSN(i)==5122769
            T(i,:)=(T(i,:)-1.465)/0.01738+25;
            % P(i,:)=P(i,:)-0.24022625;
        elseif SSN(i)==5940428
            T(i,:)=(T(i,:)-1.477)/0.01744+25;
            % P(i,:)=P(i,:)-0.2364321;
        elseif SSN(i)==5122777
            T(i,:)=(T(i,:)-1.480)/0.01739+25;
            % P(i,:)=P(i,:)-0.22593766;
        elseif SSN(i)==5940430
            T(i,:)=(T(i,:)-1.484)/0.01734+25;
            % P(i,:)=P(i,:)-0.25612138;
        else
            error('no such sensor');
        end
    end

else
    clear Pd Pxx f txx td watlev

end	% end if data exists
%%
if do_calibration
    figure();
    np = size(P,1);
    for i=1:np
        subplot(np,1,i);
        plot(P(i,:));
        hold on
        plot([1 length(P(i,:))],[1 1]*mean(P(i,:)));
    end
    calibration_table = [SSN' mean(P,2)-pamb];
end

%% 

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

%% plot raw data first
figure(1)
clf
orient tall
for i=1:4
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
% linkaxes([a11 a12 a13 a14 a15 a16],'x')
linkaxes([a11 a12 a13 a14],'x')
drawnow

%% extract segment and decimate
dec=1;% this is the decimation factor
fs=Fs/dec;
td=decimate(t(t1:t2),dec);

for i=1:4
    x=detrend(P(i,t1:t2));
    tmp = decimate(x,dec);
    if i==1
        Pd = zeros(4,length(tmp));
    end
    Pd(i,:)=tmp;
end

figure(2)
clf
orient tall
for i=1:4
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
linkaxes([b11 b12 b13 b14],'x')

%% get spectra for discrete time windows
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
par = struct(); % initialize structure
par.g = 9.80665;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
% par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
% par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 = pamb*1e5; %atmospheric pressure at equilibrium in Pa

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
        % plastic tank, 2nd version
        fprintf('Using new lab dimensions. \n')
        par.sb = 1852e-4;
        par.sl = 1; %cross-section of lateral connector in m^2
        par.L = 0; % length of lateral connector
    case 2
        % metal tank
        par.sb = pi*(0.46/2)^2; % 46 cm inner diameter
        par.sl = 1.0;
        par.L = 0;
end
% sensor 6 - top of conduit
% sensor 4 - middle of conduit
% sensor 3 - bottom of conduit
% sensor 2 - top of tank
% sensor 1 - bottom of tank

% sensor 4
% sensor 3
% sensor 2
% sensor 1

t1p = 0; t2p = Inf; % begin and end times for analysis
mask = t>=t1p & t <= t2p;

% sensor_heights = [0.0 0.0 0.0 0.0 0.0 0.0];% vector of sensor heights, measured up from bottom of tank.

% tank height
tank_height = 0.685;% internal height of tank, in m. (metal tank)

% for 4/5/2023, 
% sensor 4 appears to be garbage

% plot Rob's water level and the pressure from P1
T_tank = T(2,:); % select the temperature record for density determination
T_lookup = linspace(10,105,1000);
for i=1:1000
    rho_lookup(i) = XSteam('rhol_T',T_lookup(i));
end
rho_timeseries = interp1(T_lookup,rho_lookup,T_tank);

watlev4 = (P(4,:)-pamb)*1e5/par.g./rho_timeseries; % in m above bottom of tank
watlev3 = (P(3,:)-pamb)*1e5/par.g./rho_timeseries; % in m above bottom of tank
watlev2 = (P(2,:)-pamb)*1e5/par.g./rho_timeseries; % in m above bottom of tank
watlev1 = (P(1,:)-pamb)*1e5/par.g./rho_timeseries; % in m above bottom of tank
%p2 = p1 - rho_w*g*(water level above bottom of tank) - rho_v*g*(depth of
%steam)
tank_level = -(P(3,:) - P(1,:))*1e5/par.g ./rho_timeseries-0.0;
conduit_level = watlev1; % water level inferred from bottom sensor, ok +/- 1 cm?

figure();
subplot(3,1,1)
plot(t(mask),watlev1(mask))
hold on
plot(t(mask),watlev2(mask))
plot(t(mask),watlev3(mask))
plot(t(mask),watlev4(mask))
ylabel('level above sensor (m)')
% plot(t(mask),watlev4(mask))
legend('1','2','3','4');%
subplot(3,1,2);
plot(t(mask),tank_level(mask))
hold on
plot(t(mask),smoothdata(tank_level(mask),"gaussian",100),'k')
ylabel('Tank level (m)')
delta_2_1 = mean(watlev2(mask)-watlev1(mask))
delta_3_1 = mean(watlev3(mask)-watlev1(mask))
delta_3_2 = mean(watlev3(mask)-watlev2(mask))

plot(t(mask),rho_timeseries(mask)/1000,'r');

subplot(3,1,3);
plot(t(mask),T(1,mask)); hold on
plot(t(mask),T(2,mask));
plot(t(mask),T(3,mask));
plot(t(mask),T(4,mask));
ylabel('T (C)')
xlabel('time')

    %% loop over windows.
addpath Steam/
addpath XSteam_Matlab_v2.6/
window_s = 15; % length of window, in seconds
step_s = 1;
window_start = t1p:step_s:t2p;
Nwin = length(window_start);
dec = 20;

% arrays to hold predicted frequencies
frequencies = zeros(size(window_start));
uncertainties = zeros(size(window_start));

Pxx = zeros();
% get spectra for discrete time windows
winlens=40;	% window length, sec
winlen=winlens*fs;

% nfft=512;
% Pxx=zeros(Nwin,nfft/2+1);


iwin = 0;
for win_start = window_start
    iwin = iwin + 1;
    mask = t>= win_start & t <= win_start+window_s;

    % set xbar and ybar
    par.xbar = mean(tank_level(mask)) - (tank_height-par.H);
    par.ybar = mean(conduit_level(mask)) - (tank_height-par.H);
    par.delxy = par.ybar-par.xbar;
    [f,u]=steam_frequency(par);
    frequencies(iwin) = f;
    uncertainties(iwin) = u;

    % compute a spectrogram
    win_Pd = decimate(detrend(P(inst,mask)),dec); % decimate sensor (inst) by decimation factor.
    win_td = decimate(t(mask),dec);
    Fsd = Fs/dec;
    % [Pxx(iwin,:),f]=pmtm(win_Pd,2,nfft,Fsd);
    window = ones(size(win_Pd));
    if iwin==1
        [~,~,~,f]=mt_phs(win_Pd.*window,Fsd,1);
        Pxx = zeros(Nwin,length(f));
    end
    [Pxx(iwin,:),~,~,f]=mt_phs(win_Pd.*window,Fsd,1);

    % y=P(3,t1+(n-1)*winlens*Fs:t1+n*winlens*Fs);
    % watlev(n)=80+(mean(y')-pamb)/9.81e-04;	% water level, cm


    % figure()
    % plot(t(mask),tank_level(mask))
end

figure()
mask = t>=t1p & t<=t2p;

ax2=subplot(4,1,1);
plot(t(mask),P(inst,mask));

ax1=subplot(4,1,2:4);
pcolor(window_start,(f),log10(Pxx'));%,'edgecolor','none');
set(gca,'YScale','log');
set(gca,'ColorScale','log');
shading flat
% view(0,90);
set(gca,'YLim',[-0.7 0.0]);
% cmin=10*log10(min(min(Pxx)));
% cmax=10*log10(max(max(Pxx)));
% crange=cmax-cmin;
% cmax2=cmin+0.80*crange;
% cmin2=cmax-0.6*crange;
% set(gca,'CLim',[cmin2 cmax2]);
%set(gca,'YTickMode','manual');
%set(gca,'YTick',[10^-0.7 10^-0.6 10^-0.5 10^-0.4 10^-0.3 10^-0.2 10^-0.1 10^0])
% set(gca,'YTickLabel',[{'0.20','0.25','0.32','0.40','0.50','0.63','0.79','1.0'}]);
ylabel('Frequency, Hz')
xlabel('Time (s)')
tS=['Power spectra (dB), Sensor ',num2str(inst)];
title(tS)
hold on;
plot(window_start,frequencies,'k','LineWidth',1.5);
plot(window_start,frequencies+uncertainties,'k--')
plot(window_start,frequencies-uncertainties,'k--')
linkaxes([ax1,ax2],'x');
% loop over Rob's windows and determine water level in each one and make a
% frequency calculation.
% window_length = 30; % seconds
% for n=1:Nwin
%     ind = (1+(n-1)*winlen:n*winlen); %indices for this window
%     xbar = mean( tank_level(ind) )-;
%
% end



