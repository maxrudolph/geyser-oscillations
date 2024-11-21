clear
close all

% This script does the empty tank calibration. It reads an empty tank
% racord and writes out the calibration table to a .mat file.

% Calibration conditions: Pressure in bar:
% pamb = 1.0110777892; % calibration day (5/22) 29.85711 in Hg
% Calibration data file:
% filename = ['05-22-2024-calibration/sensor_test_empty_tank-20240522-17-11-00'];
% output_filename = ['calibration_TEST.mat'];
% sensor_plot = 1:6
% pamb = 1.03722; % during 6/18 AM calibration
% filename = ['06-20-2024/EmptyTank_RP1036p83-20240620-09-46-40'];
% output_filename = ['06-20-2024/calibration-EmptyTank_RP1036p83-20240620-09-46-40.mat']

% pamb = 1.03263;
% filename = '06-21-2024/EmptyTank_S5inConstOnTop_RP1032p63mv-20240621-09-45-58'
% output_filename = ['06-21-2024/calibration-EmptyTank_S5inConstOnTop_RP1032p63mv-20240621-09-45-58.mat']
% sensor_plot =1:6;

% pamb = 1.01901;
% filename = '10-25-2024/calibration_1019p01-20241025-10-14-35'
% output_filename = ['10-25-2024/calibration-EmptyTank_1019p01-20241025-10-14-35.mat']

% pamb = 0.99423;
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-08/calibration_Room_994p23-20241108-10-34-34';
% output_filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-08/calibration-EmptyTank_0994p23-11-08-2024'

% pamb=0.98989;
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-18-2024/Calibration_Room989p89-20241118-12-28-54'
% output_filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-18-2024/calibration-EmptyTank-Room989p89-20241118-12-28-55'

% pamb=0.99382;
% filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-19-2024/Calibration_Room993p82-20241119-10-23-10'
% output_filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-19-2024/calibration-EmptyTank_Room993p82-20241119-10-23-10.mat'

pamb = 0.98887;
filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-20-2024/calibration-EmptyTank-Room988p87-20241120-10-19-34';
output_filename = '/Volumes/GeyserData/NSFGeyserProject/SensorData/11-20-2024/calibration-EmptyTank_Room988p87-20241120';

sensor_plot=1:6;

% sensor_plot = [1 2 3 5 6]
% calibration table, with temperature values:
calibration_table = [
    % old sensors:
    % SSN    % P offset (bar)   T_offset     T_slope
    5122778  0      1.475   0.01748
    5122770  0      1.470   0.01742
    5122769  0      1.465   0.01738
    5940428  0      1.477   0.01744
    5122777  0      1.480   0.01739
    5940430  0      1.484   0.01734
    % new sensors:
    5940434	0 1.476	0.01690
    5961388	0 1.473	0.01688
    5940436	0 1.475	0.01701
    5961392	0 1.484	0.01698
    5940432	0 1.474	0.01693
    5940431	0 1.482	0.01688
    ];

% zero out the P values
calibration_table(:,2) = 0.0;

[header,P,T] = load_sensor_data(filename,calibration_table);

nsensor = length(header.pressure_sensor_serial_numbers);

% make the calibration plots
figure();
subplot(3,1,1);
for i=sensor_plot
    label = sprintf('P%d-%d',i,header.pressure_sensor_serial_numbers(i));
    h(i) = plot(P(i,:),'DisplayName',label);
    hold on
end
legend();
title('Uncalibrated pressure')

subplot(3,1,2);
for i=sensor_plot
    [f,xi] = ksdensity(P(i,:));
    plot(xi,f,'Color',get(h(i),'Color'));
    hold on
end
title('Pressure Histogram')
xlabel('Pressure (bar)')


subplot(3,1,3);
for i=sensor_plot
    label = sprintf('T%d-%d',i,header.pressure_sensor_serial_numbers(i));
    plot(T(i,:),'DisplayName',label);
    hold on
end
legend();

Pmean = mean(P,2);

P_offset = Pmean-pamb;

for i=1:nsensor
    ind(i) = find( calibration_table(:,1) == header.pressure_sensor_serial_numbers(i))
end

calibration_table(ind,2) = P_offset;
% only write out the rows correponding to calibrated sensors
calibration_table = calibration_table(ind,:);

figure()
% subplot(3,1,3);
for i=sensor_plot
    label = sprintf('T%d-%d',i,header.pressure_sensor_serial_numbers(i));
    plot(P(i,:)-P_offset(i),'DisplayName',label);
    hold on
end
legend();

save(output_filename,"calibration_table");

%% 
% re-load the sensor data using the calibration table:
load(output_filename);
[header,P,T] = load_sensor_data(filename,calibration_table);
figure()
subplot(3,1,1);
for i=1:nsensor
    label = sprintf('P%d-%d',i,header.pressure_sensor_serial_numbers(i));
    h(i) = plot(P(i,:),'DisplayName',label);
    hold on
end
legend();
title('Calibrated Pressure')
ylabel('Pressure (bar)')

subplot(3,1,2);
for i=1:nsensor
    [f,xi] = ksdensity(P(i,:));
    plot(xi,f,'Color',get(h(i),'Color'));
    hold on
end
title('Pressure Histogram')
xlabel('Pressure (bar)')

subplot(3,1,3);
for i=1:nsensor
    label = sprintf('T%d-%d',i,header.pressure_sensor_serial_numbers(i));
    plot(T(i,:),'DisplayName',label);
    hold on
end
legend();
title('Temperature (C)')

