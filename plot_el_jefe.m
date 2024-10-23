clear
close all

filedir = 'El Jefe P and T time series';
% Files MUST be ordered E,N,Z
files = {'Pts_EJ_2012_10_20_E.mat',
    'Pts_EJ_2012_10_20_N.mat',
    'Pts_EJ_2012_10_20_Z.mat'};
% Calibrations provided by Carolina Munoz-Saez by email 10/1/2024
% The calibration curve is a linear curve Mx+N
%
% ME=0.003493614;
% NE=-6650.939658;
% MZ=0.00353;
% NZ=-3196.423419;
% MN=0.005938286;
% NN=855.0403333;
%
% N is 1.5m
% Z  is 1.2m
% E  is 0.9m

slope = [0.003493614,0.005938286,0.00353];
intercept = [-6650.939658,855.0403333,-3196.423419];
label = {'E-0.9m','N-1.5m','Z-1.2m'};
figure(101);
for i=1:3
    alldata = load([filedir '/' files{i}]);
    fieldname = files{i}(1:end-4);
    data = alldata.(fieldname);
    dt = data.Time(2) - data.Time(1);
    time{i} = data.Time;
    data_uncal = data.Data;
    data_cal{i} = intercept(i) + data_uncal * slope(i);

    plot(time{i},data_cal{i},'DisplayName',label{i}); hold on;
end
legend()

% compute 1s rms
i=1
j=1;
window_length = 200; % samples
for j=1:length(data_cal{i})-window_length
    data_window = data_cal{i}(j:j+window_length);
    data_window = data_window - mean(data_window);
    data_rms(j) = sqrt( sum( data_window.^2 )/window_length );
end
drms_filt = medfilt1(data_rms,1000);
figure, plot(data_rms);
hold on
plot(drms_filt)
title('rms')

eruption_threshold = 125;
figure
subplot(2,1,1);
plot(data_rms);
ax1 = gca();
subplot(2,1,2);
plot(drms_filt > eruption_threshold);
ax2 = gca();
linkaxes([ax1 ax2],'x');

%%
drms = diff(data_rms);
figure, 
plot(data_rms);
hold on
plot(drms,'r');
[peaks,ind] = findpeaks(drms,'MinPeakProminence',2);


%%

[peaks,ind] = findpeaks(data_rms,'MinPeakProminence',600);
figure, plot(data_rms)
hold on
plot(ind,peaks,'x');

% stack on peaks
stack = zeros( min(diff(ind)),1 );
for j=1:length(ind)-1
    n = ind(j+1)-ind(j);
    stack(1:n) = data_cal{j}( ind(j):ind(j)+n-1 );
end
stack = stack/j;

figure, plot(stack);