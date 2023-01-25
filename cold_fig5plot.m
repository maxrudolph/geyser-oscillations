format long;
clear;
%close all;

addpath XSteam_Matlab_v2.6\;
addpath Ideal_Gas\;
addpath Steam\
conduit_diameter = 1; % 1 -> 1 inch, 2 -> 2 inches
if conduit_diameter == 1
    filename='cold-water-freq-filtered.csv';
elseif conduit_diameter ==2
    filename='cold-water-freq-2in-filtered.csv';
end
opts = detectImportOptions(filename);
opts.SelectedVariableNames = [1,2,6,8,14];
expdata = readmatrix(filename,opts);

fprintf('Using new lab dimensions. \n')
par.sb = 1852e-4;
par.sl = 1; %cross-section of lateral connector in m^2
par.L = 0; % length of lateral connector
par.H = 74.7e-2;
xbars = linspace(0, 74.25,length(expdata(1:27,1)))*1e-2;
correction = 2.5e-2;
if strcmp(filename,'cold-water-freq-2in-filtered.csv')
    par.sc = 20.26e-4; % cross-section of (2-inch) column in m^2
    exp_xbar = (expdata(:,1)-3)*1e-2;
    exp_ybar = (expdata(:,2)-23)*1e-2+par.H;
    %exp_xbar(9:length(expdata))=expdata(9:end,1)*1e-2 -correction;
    freqscale=2;
elseif strcmp(filename,'cold-water-freq-filtered.csv')
    par.sc = 5.07e-4; % cross-section of (1-inch) column in m^2
    exp_xbar = (expdata(1:27,1)-0)*1e-2;
    exp_ybar = (expdata(1:27,2)-0)*1e-2+par.H;
    %exp_xbar(28:length(expdata))=expdata(28:end,1)*1e-2-correction;
    freqscale=1;
end
exp_err = expdata(1:27,5)*2*pi;
exp_amplitude = expdata(1:27,3);
exp_frequencies = expdata(1:27,4)/(2*pi);
delxys = exp_ybar-exp_xbar;
par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in Pa
%ybars = delxys+par.xbar; %mean water height in conduit in m, must be greater than Sb/Sc*x0
numPar.x0 = 1e-3;
numPar.v0= 0; %initial velocity in m/s
numPar.tf = 5;
delta = 1e-6;%bar
numfrequencies = zeros(length(delxys),length(xbars));
labfrequencies = zeros(length(delxys), length(exp_xbar));
propagated_errors=zeros(length(delxys), length(exp_xbar));
for j =1:length(delxys)
    par.delxy = delxys(j);
    for k = 1:length(xbars)+length(exp_xbar)
        if k<= length(xbars)
            par.xbar = xbars(k);
            par.ybar = par.delxy + par.xbar;
        elseif k>length(xbars)
            par.xbar = exp_xbar(k-length(xbars));
            par.delxy = delxys((k-length(xbars)));
            par.ybar = exp_ybar(k-length(xbars));
        end
        
        fprintf("Using y-x= %f \n y= %f \n x= %f \n", par.delxy*1e2, par.ybar*1e2, par.xbar*1e2);
        par.Pb0= par.rho*par.g*par.delxy +par.Pa0;
        if k<= length(xbars)
            numfrequencies(j,k) = coldfreq(par);
        elseif k>length(xbars)
            labfrequencies(j,k-length(xbars))=coldfreq(par);
            %propagated_errors(j,k-length(xbars)) = 0.01; %garbage value
            propagated_errors(j,k-length(xbars)) = propagation(par,0);
        end
    end
end
% figure;
% plot(xbars*1e2,numfrequencies, 'MarkerSize',15); hold on
% for l=1:length(exp_frequencies)
%     errorbar(exp_xbar(l)*1e2,exp_frequencies(l),propagated_errors(l,l), '.', 'MarkerSize', 15);
% end
% hold off;
% title('Predicted Oscillation Frequencies, Cold Water')
% ylabel('Frequency (Hz)')
% xlabel('xbar (cm)')
% %legend('y-x=40.2 cm','y-x=61.3 cm','y-x=35.1 cm','y-x=36.7 cm','y-x=27.2 cm','y-x=49.7 cm','y-x=17.7 cm','Location','northeast')
% ylim([-0.1,1])
% xlim([0,par.H*1e2])

%% plot of predicted vs. observed frequencies
figure;
p2=errorbar(labfrequencies(1,:), exp_frequencies,exp_err,exp_err,propagated_errors(1,:),propagated_errors(1,:),'.'); hold on
p2(1).MarkerSize = 15;

cmap = colormap('jet'); % retrieve the jet colormap
% map delxy onto colors
color_variable = delxys*1e2;
colors = interp1(linspace(min(color_variable),max(color_variable),length(cmap)),cmap,color_variable);
caxis([min(color_variable) max(color_variable)])

scatter(labfrequencies(1,:),exp_frequencies,[],colors, 'filled');
a=colorbar();
ylabel(a, '$\bar y - \bar x$ (cm)','fontsize',14,'interpreter','latex', 'Rotation',90)

xave=linspace(0.4,.65*freqscale,51); plot(xave,xave,'--');
set(gca,'FontSize',14)
xlabel('Predicted Frequency (Hz)','fontsize',14)%,'interpreter','latex')
ylabel('Observed Frequency (Hz)','fontsize',14)%,'interpreter','latex')
title('Predicted vs Observed Frequency, Cold Water','fontsize',14)%'interpreter','latex')
hold off
