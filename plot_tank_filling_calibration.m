clear
close all

file_dir = '06-18-2024-calibration';

file_list = {
    'NewSensorTests_WaterLevel0_Room1033p88-20240618-14-29-32',1.03388,0.0;
    'NewSensorTests_WaterLevel10p7_Room1033p89-20240618-14-41-04',1.03389,10.7;
    'NewSensorTests_WaterLevel20p5_Room1033p94-20240618-14-50-34',1.03394,20.5;
    'NewSensorTests_WaterLevel31p0_Room1033p89-20240618-15-01-57',1.03389,31.0;
    'NewSensorTests_WaterLevel41p0_Room1034p01-20240618-15-14-07',1.03501,41.0;
    'NewSensorTests_WaterLevel51p0_Room1033p79-20240618-15-36-56',1.03379,51.0;
    'NewSensorTests_WaterLevel61p2_Room1033p75-20240618-15-45-43',1.03375,61.2;
    'NewSensorTests_WaterLevel69p0_Room1033p77-20240618-15-54-04',1.03377,69.0;   
    }

calibration_file = 'calibration_06182024-Sensor1-installed.mat';
load(calibration_file);

nfile = size(file_list,1);
g = 9.81;
for ifile = 1:nfile
    filename = [file_dir '/' file_list{ifile,1}];
    pamb = file_list{ifile,2};

    [header,P,T] = load_sensor_data(filename,calibration_table);
    nsensor = length(header.pressure_sensor_serial_numbers);

    T_sensor = 1;
    T_tank = T(T_sensor,:);
    T_lookup = linspace(10,105,1000);
    for i=1:1000
        rho_lookup(i) = XSteam('rhol_T',T_lookup(i));
    end
    rho_timeseries = interp1(T_lookup,rho_lookup,T_tank);

    P1_mean(ifile) = mean(P(1,:)-pamb);
    mean_level(ifile) = mean( (P(1,:)-pamb)*1e5./rho_timeseries/g )*100;
end

fill_level = [file_list{:,3}];
figure,
% plot(fill_level,P1_mean,'x');
plot(mean_level,fill_level,'x');
hold on
plot(mean_level,mean_level,'--')

%% fit a line to the data
[coef,ci] = polyfit(mean_level,fill_level,1);
hold on
plot(mean_level,coef(1)*mean_level+coef(2),'b:')