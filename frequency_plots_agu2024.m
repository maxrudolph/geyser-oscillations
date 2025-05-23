clear
close all
addpath XSteam_Matlab_v2.6/
addpath Steam/

% load and plot lab data with calibrations...
% p_amb = 1.03707;
% filename = ['06-18-2024-calibration/NewSensorTests_WaterLevel69p0_Room1033p77-20240618-15-54-04']
% filename = '06-19-2024/AddingWaterFromTop_1p5kWheat_RP1036p70-20240619-14-49-00';
% calibration_file = 'calibration-test.mat';
% sensors_plot = [1 2 3 4 5];
% filename = '06-19-2024/ConeEruption_AllHeaters_1Lpmin-20240619-17-00-56'
% calibration_file = '06-19-2024/calibration-06-19-2024.mat'

% p_amb = 1.03683;
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/06-20-2024/ConeEruptionTopConstriction_1p2LpMin_AllHeaters_-20240620-15-17-17';
% calibration_file = '/Volumes/GeyserData/NSFGeyserProject/SensorData/06-20-2024/calibration-EmptyTank_RP1036p83-20240620-09-46-40.mat';
% dec = 1;

p_amb = 1.03683;
filename = '06-20-2024/ConeEruptionTopConstriction_1p2LpMin_AllHeaters_-20240620-15-17-17';
calibration_file = '06-20-2024/calibration-EmptyTank_RP1036p83-20240620-09-46-40.mat';
dec=1;

% calibration_file = '10-25-2024/calibration-EmptyTank_1019p01-20241025-10-14-35.mat'
% filename = '10-25-2024/MidConstriction_Cone_Stage1-20241025-11-44-34'
% filename = '10-25-2024/MidConstriction_Cone_Stage5-20241025-14-19-02'

% p_amb = 0.99423;
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-08-2024/cycles_no_constrictions-20241108-12-10-06'
% calibration_file = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-08-2024/calibration-EmptyTank_0994p23-11-08-2024.mat'


% p_amb = 0.98989;
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-18-2024/Eruption_Pool_TopConstriction_try2-20241118-14-24-43'
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-18-2024/AbovePoolCone_BelowPoolConstriction_try1-20241118-15-43-22'
% calibration_file = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-18-2024/calibration-EmptyTank-Room989p89-20241118-12-28-55.mat'

sensors_plot = [1 2 3 4 5 ];
% p_amb = 0;
% filename = ['06-18-2024-calibration/NewSensorTests-20240618-11-49-47']
% calibration_file = 'calibration_06182024.mat';

load(calibration_file)

[header,P,T] = load_sensor_data(filename,calibration_table);
fs = header.sampling_rates(1);
nsensor = length(header.pressure_sensor_serial_numbers);
t = (1:size(P,2))/fs;

%%
% dec = 1;
n1 = length( decimate(double(P(1,:)),dec));
Pd = zeros(size(P,1),n1);
Td = Pd;
for i=1:size(P,1)
    Td(i,:) = decimate(double(T(i,:)),dec);
    Pd(i,:) = decimate(double(P(i,:)),dec);
end
fsd = fs/dec;
td = (1:size(Pd,2))/fsd;

clear P T;

%%
figure();
subplot(3,1,1);
for i=sensors_plot
    label = sprintf('P%d-%d',i,header.pressure_sensor_serial_numbers(i));
    h(i) = plot(Pd(i,:)-p_amb,'DisplayName',label);
    hold on
end
legend();
title('Calibrated Pressure')
ylabel('Pressure (bar)')

subplot(3,1,2);
for i=sensors_plot
    % [f,xi] = ksdensity(Pd(i,:)-p_amb);
    % plot(xi,f,'Color',get(h(i),'Color'));
    hold on
end
title('Pressure Histogram')
xlabel('Pressure (bar)')

subplot(3,1,3);
for i=sensors_plot
    label = sprintf('T%d-%d',i,header.pressure_sensor_serial_numbers(i));
    plot(Td(i,:),'DisplayName',label);
    hold on
end
legend();
title('Temperature (C)')

%% use P4 to compute a 1S moving window - to detect eruptions
% movRMS = dsp.MovingRMS(fs);
% winlen = 1000;
% P4d = Pd(4,:);
%
% % win_start = 1:shift:(length(P4d)-winlen);
% p4_rms = zeros(size(P4d));
%
% for j=1:length(P4d)
%     if j < winlen/2+1
%         chunk = P4d(1:(winlen+1));
%     elseif j > length(P4d) - winlen/2
%         chunk = P4d(end-winlen:end);
%     else
%         chunk = P4d(j-winlen/2:j+winlen/2);
%     end
%     assert(length(chunk) == winlen+1)
%     % chunk = P4d(start:start+winlen-1);
%     p4_rms(j) = sqrt( sum( (chunk-mean(chunk)).^2 )/(winlen+1) );
%
% end
% eruption_threshold = .006;
%
% erupt_mask = movmedian(p4_rms,5000) > eruption_threshold;
%
%
erupt_mask = movmedian(Td(4,:),5000) > 65;
figure, plot(td,erupt_mask);


%% plot using water equivalent
T_sensor = 4;
T_tank = Td(T_sensor,:);
T_lookup = linspace(10,110,1000);
for i=1:1000
    rho_lookup(i) = XSteam('rhol_T',T_lookup(i));
end
rho_timeseries = interp1(T_lookup,rho_lookup,T_tank);

water_level = (Pd-p_amb)*1e5/9.81./rho_timeseries;
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
    plot(Td(i,:),'DisplayName',label);
    hold on
end
legend();
title('Temperature (C)')

%% Set up parameters for lab geyser
geometry = 2;
conduit_diameter = 1; % 1 -> 1 inch, 2 -> 2 inches
par = struct(); % initialize structure
par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
% par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
% par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 = p_amb*1e5; %atmospheric pressure at equilibrium in Pa

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
    case 2 % 2024 - metal tank configuration
        par.sb = pi*0.23^2;
        par.sl = 5.0671e-04;
        par.L = 0.1665; % double check
end
conduit_level = water_level(1,:);
tank_level = -(water_level(2,:) - water_level(1,:)); % measured upward from bottom of tank


%% Spectrogram
addpath Steam/
addpath XSteam_Matlab_v2.6/
inst=[2] %sensors_plot 3=bottom of conduit
window_s = 10; % length of window, in seconds
step_s = 5;
% inst=4;
dec = 1; % this is an additional decimation beyond the first decimation.
td = (1:size(Pd,2))*(1/fsd);
t1p = 0;
t2p = td(end)-window_s;
window_start = t1p:step_s:t2p;

Nwin = length(window_start);

% arrays to hold predicted frequencies
frequencies = zeros(size(window_start));
uncertainties = zeros(size(window_start));
f_acoustic = zeros(4,Nwin);
cs = zeros(size(window_start));
Pxx = zeros();
% get spectra for discrete time windows

iwin = 0;
for win_start = window_start
    iwin = iwin + 1;
    mask = td>= win_start & td <= win_start+window_s;

    % set xbar and ybar
    % tank_level = -(P(2,:) - P(1,:))*1e5/1000/9.81;
    tank_height = 0.76;% internal height of tank, in m.


    par.xbar = mean(tank_level(mask)) - (tank_height-par.H);
    par.ybar = mean(conduit_level(mask)) - (tank_height-par.H);
    par.delxy = par.ybar-par.xbar;
    if ~isnan(par.xbar)
        [f,u]=steam_frequency(par);
        % if imag(f) ~= 0
        % keyboard;
        % end
        frequencies(iwin) = f;
        uncertainties(iwin) = u;
    else
        frequencies(iwin) = NaN;
        uncertainties(iwin) = NaN;
    end
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
end
%% model predictions - acoustic modes:
iwin = 0;
tank_levels = zeros(length(window_start),1);
conduit_levels = zeros(length(window_start),1);
window_erupt = zeros(size(window_start));
Lsave = zeros(length(window_start),1);
for win_start = window_start
    iwin = iwin + 1;
    mask = td>= win_start & td <= win_start+window_s;

    tank_height = 0.76;% internal height of tank, in m.

    %acoustic modes
    % compute sound speed
    cs(iwin) = XSteam('w_pT', mean(Pd(inst,mask)),mean(Td(inst,mask)));
    this_tank_level = mean(tank_level(1,mask));
    tank_levels(iwin) = this_tank_level;
    conduit_levels(iwin) = mean(conduit_level(mask));
    window_erupt(iwin) = median(erupt_mask(mask));

    % f_acoustic(:,iwin) = [1 3 5 7].*cs(iwin)./(4*(mean(conduit_level(mask))+(tank_height-par.H)));
    L = (mean(conduit_level(mask) + par.L) + 0.33*0 + this_tank_level +0 ) + 0.0;
    % f_acoustic(:,iwin) = (2*[1 2 3 4]-1) * cs(iwin) / (4*L); % forced at one end, closed at the other
    L = 1.7960*(mean(conduit_level(mask)) - 0.5063);
    Lsave(iwin) = L;
    f_acoustic(:,iwin) = ([1 2 3 4]) * cs(iwin) / (4*L); % forced at one end, closed at the other
    % forced at one end, open at the other
    % a = sqrt(par.sb/pi);% pi r^2 = sc -> r = sqrt(sc/pi)
    % f_acoustic(:,iwin) = [1 2 3 4] * cs(iwin) / (2 * (L + 0.6*a));

    %Helmholtz resonator - doesn't do a good job.
    % a = sqrt(par.sb/pi);% pi r^2 = sc -> r = sqrt(sc/pi) - radius of resonant cavity
    % Lprime = mean(conduit_level(mask)) + 1.4*a;
    % V = mean(tank_level(mask))*par.sb;
    % omega_0 = cs(iwin) * sqrt(par.sc/Lprime/V);
    % f_acoustic(1,iwin) = omega_0/(2*pi);
end
%% make the figure
figure()
set(gcf,'name',['Sensor ' num2str(inst)]);
mask = td>=t1p & td<=t2p;

ax2=subplot(4,1,1);
plot(td(mask),Pd(inst,mask),'k');
set(gca,'XTickLabel',[]);
ylabel('Pressure (bar)','FontSize',16)
colormap parula

ax1=subplot(4,1,2:4);
pcolor(window_start,(f),(Pxx'));%,'edgecolor','none');
set(gca,'YScale','log');
set(gca,'ColorScale','log');
shading flat
% view(0,90);
% set(gca,'YLim',[-0.7 2]);
% cmin=10*log10(min(min(Pxx)));
% cmax=10*log10(max(max(Pxx)));
% crange=cmax-cmin;
% cmax2=cmin+0.80*crange;
% cmin2=cmax-0.6*crange;
% set(gca,'CLim',[cmin2 cmax2]);
%set(gca,'YTickMode','manual');
%set(gca,'YTick',[10^-0.7 10^-0.6 10^-0.5 10^-0.4 10^-0.3 10^-0.2 10^-0.1 10^0])
% set(gca,'YTickLabel',[{'0.20','0.25','0.32','0.40','0.50','0.63','0.79','1.0'}]);
ylabel('Frequency (Hz)','FontSize',16)
xlabel('Time (s)','FontSize',16)
% tS=['Power spectra (dB), Sensor ',num2str(inst)];
% title(tS)

hold on;

frequencies(window_erupt > 0.5) = NaN;
plot(window_start,frequencies,'k','LineWidth',1.5);
plot(window_start,frequencies+uncertainties,'k--')
plot(window_start,frequencies-uncertainties,'k--')
f_acoustic(:,window_erupt > 0.5) = NaN;
plot(window_start,f_acoustic,'k');
hcb=colorbar();
hcb.Label.String = 'Power (dB/Hz)';
hcb.Label.FontSize = 14;

clim([1e-12 1e-4]);
linkaxes([ax1,ax2],'x');
set(ax2,'XLim',[3200 9400])
% drawnow();
ax2.Position(3) = ax1.Position(3);
set(gca,'Layer','top')
tmp = ax2.Position(2) + ax2.Position(4);
ax2.Position(2) = ax1.Position(2) + ax1.Position(4) + 0.02;
ax2.Position(4) = tmp-ax2.Position(2);
% plot the conduit 'effective length'
figure
% plot(td,)
plot(window_start,Lsave);
title('effective length')
hold on
plot(window_start,conduit_levels,'r')


%% plot conduit and reservoir pressures and temperatures
plot_order = [5 4 3 2 1];
labels = {'Conduit Top','Conduit Middle','Conduit Bottom','Bubble Trap','Tank Bottom'};
figure()
for i=1:5
    subplot(5,1,i);
    yyaxis left
    plot(td,Pd(plot_order(i),:),'-','DisplayName','P')
    % hold on
    % plot(td,Pd(2,:),'--','DisplayName','P-Tank Top')
    ylabel('Pressure (bar)');
    yyaxis right
    plot(td,Td(plot_order(i),:),'DisplayName','Bottom')
    % hold on
    % plot(td,Td(2,:),'DisplayName','T-Tank Top')
    % legend()
    set(gca,'XLim',[3200 9500])
    ylabel('Temperature (K)')
    text(0.05,0.75,labels{i},'Units','normalized','FontSize',16)
    set(gca,'XTickLabel',[]);
end
% exportgraphics(gcf,'all-signals.pdf')
% xlabel('Time (s)')
% end