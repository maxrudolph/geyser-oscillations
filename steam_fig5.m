format long;
clear;
close all;

addpath XSteam_Matlab_v2.6\;
addpath Ideal_Gas\;
addpath Steam\
geometry = 1;
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
        par.sb = 2026.83e-4;
        par.sc = 5.07e-4; % cross-section of (1-inch) column in m^2
        par.sl = 1; %cross-section of lateral connector in m^2
        par.L = 0; % length of lateral connector
        par.H = 74.7e-2;
end
xbars = linspace(0, 74.25,101)*1e-2;
labxbars = [38,51,51,51,60,66]*1e-2;
labybars = [3.5,37.6,11.4,13,12.5,41]*1e-2 + par.H;
delxys = [40.2,61.3,35.1,36.7,27.2,49.7]*1e-2;
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
labfrequencies = zeros(length(delxys), length(labxbars));
propagated_errors=zeros(length(delxys), length(labxbars));
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
        par.m=par.Vol_0 / XSteam('vV_p', P_0/1e5); %vapor mass in kg

        xv = 1 - 1e-2;
        sv=XSteam('sV_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
        sl=XSteam('sL_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
        s = xv*sv + (1-xv)*sl; % gives specific entropy of water+vapor mixture in J/kg/C
        V0 = XSteam('V_ps',P_0/1e5,s/1e3); % volume per unit mass
        par.m = par.Vol_0 / V0;

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
            numfrequencies(j,k) = frequency_calculator(par);
        elseif k>length(xbars)
            labfrequencies(j,k-length(xbars))=frequency_calculator(par);
            propagated_errors(j,k-length(xbars)) = propagation(par,1);
        end
    end
end
err = [0.0003075117811,0.0006765997287,0.0003365434388,0.000366547911,0.0005258494879,0.0004347155603]*2*pi;
figure;
exp_frequencies = [3.82,3.19,3.69,3.64,3.75,3.27]/(2*pi);
plot(xbars*1e2,numfrequencies, 'MarkerSize',15); hold on
for l=1:6
    errorbar(labxbars(l)*1e2,exp_frequencies(l),propagated_errors(l,l), '.', 'MarkerSize', 15);
end
hold off;
title('Predicted Oscillation Frequencies')
ylabel('Frequency (Hz)')
xlabel('xbar (cm)')
%legend('y-x=40.2 cm','y-x=61.3 cm','y-x=35.1 cm','y-x=36.7 cm','y-x=27.2 cm','y-x=49.7 cm','y-x=17.7 cm','Location','northeast')
ylim([-0.1,2])
xlim([0,par.H*1e2])

figure;
p2=errorbar(labfrequencies(1,:), exp_frequencies,err,err,propagated_errors(1,:),propagated_errors(1,:),'.'); hold on
p2(1).MarkerSize = 15;
cmap = colormap('jet'); % retrieve the jet colormap
% map delxy onto colors
color_variable = delxys*1e2;
colors = interp1(linspace(min(color_variable),max(color_variable),length(cmap)),cmap,color_variable);
caxis([min(color_variable) max(color_variable)])

scatter(labfrequencies(1,:),exp_frequencies,[],colors, 'filled');
a=colorbar();
ylabel(a, '$\bar y - \bar x$ (cm)','fontsize',14,'interpreter','latex', 'Rotation',90)
xave=linspace(0.4,0.65,51); plot(xave,xave,'--');
xlabel('Predicted Frequency (Hz)','fontsize',14,'interpreter','latex')
ylabel('Observed Frequency (Hz)','fontsize',14,'interpreter','latex')
title('Predicted vs Observed Frequency','fontsize',14,'interpreter','latex')
hold off
