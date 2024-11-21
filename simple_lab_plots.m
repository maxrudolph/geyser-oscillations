clear
close all
addpath XSteam_Matlab_v2.6/

% load and plot lab data with calibrations...
% p_amb = 1.03683;
% filename = ['06-18-2024-calibration/NewSensorTests_WaterLevel69p0_Room1033p77-20240618-15-54-04']
% filename = '06-20-2024/EmptyTank_RP1036p83-20240620-09-46-40';
% calibration_file = '06-20-2024/calibration-EmptyTank_RP1036p72-20240620-09-20-03.mat';

% filename = '06-19-2024/ConeEruption_AllHeaters_1Lpmin-20240619-17-00-56'
% calibration_file = '06-19-2024/calibration-06-19-2024.mat'
sensors_plot = [1 2 3 4 5 6];
% p_amb = 1.03263;
% filename = '06-21-2024/EmptyTank_S5inConstOnTop_RP1032p63mv-20240621-09-45-58'
% calibration_file = '06-20-2024/calibration-EmptyTank_RP1036p72-20240620-09-20-03.mat';

% p_amb = 0;
% filename = ['06-18-2024-calibration/NewSensorTests-20240618-11-49-47']
% calibration_file = 'calibration_06182024.mat';



% p_amb = 0.98989;
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-18-2024/Eruption_Pool_TopConstriction_try2-20241118-14-24-43'
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-18-2024/AbovePoolCone_BelowPoolConstriction_try1-20241118-15-43-22'
% calibration_file = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-18-2024/calibration-EmptyTank-Room989p89-20241118-12-28-55.mat'


% p_amb = 0.99382;
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-19-2024/Pool_MidConstriction_-20241119-12-08-17'
% calibration_file = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-19-2024/calibration-EmptyTank_Room993p82-20241119-10-23-10.mat'

filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-20-2024/cone-topconstriction-try3-20241120-15-33-06'
calibration_file = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-20-2024/calibration-EmptyTank_Room988p87-20241120.mat'
p_amb = 0.98887;


load(calibration_file)

[header,P,T] = load_sensor_data(filename,calibration_table);
nsensor = length(header.pressure_sensor_serial_numbers);
%% 
dec = 10;
Pd = decimate(double(P(1,:)),dec);
n1 = length(Pd);
Pd = zeros(size(P,1),n1);
Td = Pd;
for i=1:size(Pd,1)
    Td(i,:) = decimate(double(T(i,:)),dec);
    Pd(i,:) = decimate(double(P(i,:)),dec);
end
dt = 1/header.sampling_rates(1)*dec;
td = dt*(0:(size(Pd,2)-1));

clear P T;

%% 
figure();
subplot(3,1,1);
for i=sensors_plot
    label = sprintf('P%d-%d',i,header.pressure_sensor_serial_numbers(i));
    h(i) = plot(td,Pd(i,:)-p_amb,'DisplayName',label);
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
    plot(td,Td(i,:),'DisplayName',label);
    hold on
end
legend();
title('Temperature (C)')

%% plot using water equivalent
T_sensor = 2;
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
% for i=sensors_plot
    % [f,xi] = ksdensity(water_level(i,:));
    % plot(xi,f,'Color',get(h(i),'Color'));
    % hold on
% end
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

%% calculate some statistics
mean(water_level(1,:))*100
sqrt(var(water_level(1,:)*100))