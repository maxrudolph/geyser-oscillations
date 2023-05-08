format long;
clear;
%close all;

addpath XSteam_Matlab_v2.6\;
addpath Ideal_Gas\;
addpath Steam\
geometry = 1;
conduit_diameter = 1; % 1 -> 1 inch, 2 -> 2 inches
switch conduit_diameter
    case 1
        filename='hot-water-freq-filtered-1in_short.csv';
        par.sc = pi*(1/2*2.54/100)^2; %cross-section of 1-inch tube in cm^2
        par.H = 31.8e-2; % "new" short 1 inch pickup in m
        %par.H = 58e-2; % old, longer 1 inch pickup in m
    case 2
        filename='hot-water-freq-filtered-2in_v3.csv';
        par.sc = pi*(2/2*2.54/100)^2; %cross-section of 2-inch tube in cm^2
        par.H = (53-1.8)*1e-2; %Bubble trap height in m
end

switch geometry
    case 0
        fprintf('Using old lab dimensions. \n');
        par.sb = 613e-4; % cross-section of bubble trap in m^2
        par.sc = 5.07e-4; % cross-section of (1-inch) column in m^2
        par.sl = 7.92e-4; %cross-section of lateral connector in m^2
        par.L = 7.37e-2; % length of lateral connector in m
        par.H = 37.15e-2; % height of bubble trap in m
    case 1
        fprintf('Using new lab dimensions. \n')
        par.sb = 1852e-4;
        par.sl = 1; %cross-section of lateral connector in m^2
        par.L = 0; % length of lateral connector
end

opts = detectImportOptions(filename);
opts.SelectedVariableNames = [1,2,5,7,8];
expdata = readmatrix(filename,opts);

xbars = linspace(0, par.H,101)*1e-2;
labxbars = (expdata(:,1))*1e-2;
labybars = (expdata(:,2))*1e-2;
labsensors = expdata(:,5);
delxys = labybars-labxbars;


par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in Pa

delta = 1e-6;%bar
numfrequencies = zeros(1,length(xbars));
labfrequencies = zeros(1, length(labxbars));
propagated_errors=zeros(1, length(labxbars));
for j =1:length(delxys)
    par.delxy = delxys(j);
    for k = 1:length(xbars)+length(labxbars)
        if k<= length(xbars)
            par.xbar = xbars(k);
        elseif k>length(xbars)
            par.xbar = labxbars(k-length(xbars));
            par.delxy = delxys((k-length(xbars)));
        end
        par.ybar = par.delxy + par.xbar;
        %% Make a lookup table with pre-computed thermodynamic properties.
        par.Vol_0 = par.sb*(par.H-par.xbar); %initial (equilibrium) volume in m^3
        P_0 = par.rho*par.g*(par.delxy)+par.Pa0; %initial pressure in (Pascal)

        xv = 1 - 1e-2;
        sv=XSteam('sV_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
        sl=XSteam('sL_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
        s = xv*sv + (1-xv)*sl; % gives specific entropy of water+vapor mixture in J/kg/C
        V0 = XSteam('V_ps',P_0/1e5,s/1e3); % volume per unit mass
        par.m = par.Vol_0 / V0; %vapor mass in kg
        
        P=linspace((P_0-9e4)/1e5,(P_0+9e4)/1e5,51); %Range of pressure for lookup table (bar)
        vV=zeros(1,length(P));  % vapor volume (m^3)
        du_dp = zeros(1,length(P)); %dU/dP in (kJ/K/Pa)
        dv_dp = zeros(1,length(P)); %dV/dP in (m^3/Pa)
        for i=1:length(P)
            vV(i) = XSteam('v_ps',P(i),s/1e3)*par.m; %volume in m^3
            du_dp(i) = (XSteam('u_ps',P(i)+delta,s/1e3)-XSteam('u_ps',P(i)-delta,s/1e3))*1e3*par.m/(2*delta*1e5);% enthalpy in mks units, pressure in Pa
            dv_dp(i) = (XSteam('v_ps',P(i)+delta,s/1e3)-XSteam('v_ps',P(i)-delta,s/1e3))*par.m/(2*delta*1e5);% m^3/Pa
        end
        [~,i] = sort(vV);
        par.Fdudp = griddedInterpolant(vV(i),du_dp(i),'linear','none');
        par.FdPdv = griddedInterpolant(vV(i),1./dv_dp(i),'linear','none');
        if k<= length(xbars)
            numfrequencies(1,k) = frequency_calculator(par);
        elseif k>length(xbars)
            [labfrequencies(1,k-length(xbars)), propagated_errors(1,k-length(xbars)) ]=frequency_calculator(par);
            %propagated_errors(j,k-length(xbars)) = propagation(par,1);
        end
    end
    fprintf('Using %d out of %d \n',j, length(delxys));
end
err = expdata(:,4);
exp_frequencies = expdata(:,3)/(2*pi);

% figure;
% plot(xbars*1e2,numfrequencies, 'MarkerSize',15); hold on
% for l=1:6
%     errorbar(labxbars(l)*1e2,exp_frequencies(l),propagated_errors(l,l), '.', 'MarkerSize', 15);
% end
% hold off;
% title('Predicted Oscillation Frequencies')
% ylabel('Frequency (Hz)')
% xlabel('xbar (cm)')
% %legend('y-x=40.2 cm','y-x=61.3 cm','y-x=35.1 cm','y-x=36.7 cm','y-x=27.2 cm','y-x=49.7 cm','y-x=17.7 cm','Location','northeast')
% ylim([-0.1,2])
% xlim([0,par.H*1e2])

figure;
p2=errorbar(labfrequencies(1,:), exp_frequencies,err,err,propagated_errors(1,:),propagated_errors(1,:),'.'); hold on
p2(1).MarkerSize = 15;
cmap = colormap('jet'); % retrieve the jet colormap
% map delxy onto colors
color_variable = (delxys*1e2);
colors = interp1(linspace(min(color_variable),max(color_variable),length(cmap)),cmap,color_variable);
caxis([min(color_variable) max(color_variable)])

scatter(labfrequencies(1,:),exp_frequencies,[],colors, 'filled');
a=colorbar();
ylabel(a, 'ybar-xbar','fontsize',14,'interpreter','latex', 'Rotation',90)
xave=linspace(min(labfrequencies)-0.2,max(labfrequencies)+0.2,5); plot(xave,xave,'--');
xlabel('Predicted Frequency (Hz)','fontsize',14,'interpreter','latex')
ylabel('Observed Frequency (Hz)','fontsize',14,'interpreter','latex')
title('Predicted vs Observed Frequency','fontsize',14,'interpreter','latex')
hold off
