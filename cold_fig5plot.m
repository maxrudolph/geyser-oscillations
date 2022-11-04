format long;
clear;
close all;

addpath XSteam_Matlab_v2.6;
addpath Ideal_Gas;
addpath Steam;
opts = detectImportOptions('cold-water-freq-filtered.csv');
opts.SelectedVariableNames = [1,2,8,14];
expdata = readmatrix('cold-water-freq-filtered.csv',opts);

fprintf('Using new lab dimensions. \n')
par.sb = 1928e-4;
par.sc = 5.07e-4;   % cross-section of (1-inch) column in m^2
par.sl = 1;         %cross-section of lateral connector in m^2
par.L = 0;          % length of lateral connector (m)
par.H = 74.7e-2;    %height of bubble trap (m)
xbars = linspace(0, 74.25,101)*1e-2;
exp_xbar = expdata(:,1)*1e-2;
exp_ybar = expdata(:,2)*1e-2+par.H;
delxys = exp_ybar-exp_xbar;
par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in Pa
%ybars = delxys+par.xbar; %mean water height in conduit in m, must be greater than Sb/Sc*x0
numPar.x0 = 1e-4;
numPar.v0= 0; %initial velocity in m/s
numPar.tf = 5;
delta = 1e-6;%bar
numfrequencies = zeros(length(delxys),length(xbars));
labfrequencies = zeros(length(delxys), length(exp_xbar));
propagated_errors=zeros(length(delxys), length(exp_xbar));
for j =1:length(delxys)
    par.delxy = delxys(j);
    fprintf("Using %f \n", par.delxy*1e2);
    for k = 1:length(xbars)+length(exp_xbar)
        if k<= length(xbars)
            par.xbar = xbars(k);
        elseif k>length(xbars)
            par.xbar = exp_xbar(k-length(xbars));
            par.delxy = delxys((k-length(xbars)));
        end
        par.ybar = par.delxy + par.xbar;
        par.Pb0= par.rho*par.g*par.delxy +par.Pa0;
        if k<= length(xbars)
            numfrequencies(j,k) = coldfreq(par);
        elseif k>length(xbars)
            labfrequencies(j,k-length(xbars))=coldfreq(par);
            propagated_errors(j,k-length(xbars)) = 0.01; %garbage value
        end
    end
end
exp_err = expdata(:,4)*2*pi;

% make a plot showing predicted frequencies
figure;
exp_frequencies = expdata(:,3)/(2*pi);
plot(xbars*1e2,numfrequencies, 'MarkerSize',15); hold on
for l=1:length(exp_frequencies)
    errorbar(exp_xbar(l)*1e2,exp_frequencies(l),propagated_errors(l,l), '.', 'MarkerSize', 15);
end
hold off;
title('Predicted Oscillation Frequencies, Cold Water')
ylabel('Frequency (Hz)')
xlabel('xbar (cm)')
legend('y-x=40.2 cm','y-x=61.3 cm','y-x=35.1 cm','y-x=36.7 cm','y-x=27.2 cm','y-x=49.7 cm','y-x=17.7 cm','Location','northeast')
ylim([-0.1,1])
xlim([0,par.H*1e2])

%% plot of predicted vs. observed frequencies
figure;
p2=errorbar(labfrequencies(1,:), exp_frequencies,exp_err,exp_err,propagated_errors(1,:),propagated_errors(1,:),'.'); hold on
p2(1).MarkerSize = 15;

cmap = colormap('jet'); % retrieve the jet colormap
% map delxy onto colors
color_variable = delxys*100;
colors = interp1(linspace(min(color_variable),max(color_variable),length(cmap)),cmap,color_variable);
caxis([min(color_variable) max(color_variable)])

scatter(labfrequencies(1,:),exp_frequencies,[],colors, 'filled');
a=colorbar();
ylabel(a, '$\overline y$ - $\overline x$ (cm)','fontsize',14,'interpreter','latex', 'Rotation',90)
xave=linspace(0.4,0.7,51); plot(xave,xave,'--');
xlabel('Predicted Frequency (Hz)')
ylabel('Observed Frequency (Hz)')
title('Predicted vs Observed Frequency, Cold Water')
hold off
