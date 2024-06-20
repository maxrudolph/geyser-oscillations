clear
% close all

% load and plot lab data with calibrations...
p_amb = 1.03707;
% filename = ['06-18-2024-calibration/NewSensorTests_WaterLevel69p0_Room1033p77-20240618-15-54-04']
filename = '06-19-2024/AddingWaterFromTop_1p5kWheat_RP1036p70-20240619-14-49-00';
calibration_file = 'calibration-test.mat';
sensors_plot = [1 2 3 4 5];
% p_amb = 0;
% filename = ['06-18-2024-calibration/NewSensorTests-20240618-11-49-47']
% calibration_file = 'calibration_06182024.mat';

load(calibration_file)

[header,P,T] = load_sensor_data(filename,calibration_table);
fs = header.sampling_rates(1);
nsensor = length(header.pressure_sensor_serial_numbers);
t = (1:size(P,2))/fs;

%% 
dec = 3;
n1 = length( decimate(P(1,:),dec));
Pd = zeros(size(P,1),n1);
Td = Pd;
for i=1:size(P,1)
    Td(i,:) = decimate(T(i,:),dec);
    Pd(i,:) = decimate(P(i,:),dec);
end
fsd = fs/dec;
td = (1:size(Pd,2))/fsd;

%% 
figure();
subplot(3,1,1);
for i=sensors_plot
    label = sprintf('P%d-%d',i,header.pressure_sensor_serial_numbers(i));
    h(i) = plot(P(i,:)-p_amb,'DisplayName',label);
    hold on
end
legend();
title('Calibrated Pressure')
ylabel('Pressure (bar)')

subplot(3,1,2);
for i=sensors_plot
    [f,xi] = ksdensity(P(i,:)-p_amb);
    plot(xi,f,'Color',get(h(i),'Color'));
    hold on
end
title('Pressure Histogram')
xlabel('Pressure (bar)')

subplot(3,1,3);
for i=sensors_plot
    label = sprintf('T%d-%d',i,header.pressure_sensor_serial_numbers(i));
    plot(T(i,:),'DisplayName',label);
    hold on
end
legend();
title('Temperature (C)')

%% plot using water equivalent
T_sensor = 1;
T_tank = T(T_sensor,:);
T_lookup = linspace(10,105,1000);
for i=1:1000
    rho_lookup(i) = XSteam('rhol_T',T_lookup(i));
end
rho_timeseries = interp1(T_lookup,rho_lookup,T_tank);

water_level = (P-p_amb)*1e5/9.81./rho_timeseries;
% mean_others = mean(mean(water_level(2:6),2));
% water_level = water_level-mean_others;

figure()
subplot(3,1,1);
for i=sensors_plot
    label = sprintf('P%d-%d',i,header.pressure_sensor_serial_numbers(i));
    h(i) = plot(water_level(i,:),'DisplayName',label);
    hold on
end
legend();
title('Calibrated Pressure')
ylabel('level (m)')

subplot(3,1,2);
for i=sensors_plot
    [f,xi] = ksdensity(water_level(i,:));
    plot(xi,f,'Color',get(h(i),'Color'));
    hold on
end
title('Level Histogram')
xlabel('Level (m)')

subplot(3,1,3);
for i=sensors_plot
    label = sprintf('T%d-%d',i,header.pressure_sensor_serial_numbers(i));
    plot(T(i,:),'DisplayName',label);
    hold on
end
legend();
title('Temperature (C)')

%% Spectrogram
addpath Steam/
addpath XSteam_Matlab_v2.6/
window_s = 30; % length of window, in seconds
step_s = 1;
inst=3;
dec = 1;
td = (1:size(Pd,2))*(1/fsd);
t1p = 0;
t2p = td(end)-window_s;
window_start = t1p:step_s:t2p;

Nwin = length(window_start);

% arrays to hold predicted frequencies
frequencies = zeros(size(window_start));
uncertainties = zeros(size(window_start));

Pxx = zeros();
% get spectra for discrete time windows

iwin = 0;
for win_start = window_start
    iwin = iwin + 1;
    mask = td>= win_start & td <= win_start+window_s;

    % set xbar and ybar
    % par.xbar = mean(tank_level(mask)) - (tank_height-par.H);
    % par.ybar = mean(conduit_level(mask)) - (tank_height-par.H);
    % par.delxy = par.ybar-par.xbar;
    % [f,u]=steam_frequency(par);
    % frequencies(iwin) = f;
    % uncertainties(iwin) = u;

    % compute a spectrogram
    win_Pd = decimate(detrend(Pd(inst,mask)),dec); % decimate sensor (inst) by decimation factor.
    win_td = decimate(td(mask),dec);
    Fsd = fsd/dec;
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
plot(td(mask),Pd(inst,mask));

ax1=subplot(4,1,2:4);
pcolor(window_start,(f),log10(Pxx'));%,'edgecolor','none');
set(gca,'YScale','log');
set(gca,'ColorScale','log');
shading flat
% view(0,90);
set(gca,'YLim',[-0.7 2]);
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