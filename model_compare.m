format long;
clear;
close all;

addpath XSteam_Matlab_v2.6\;
addpath Ideal_Gas\;
addpath Steam\
par.sb = 2026.83e-4;
par.sc = 5.07e-4; % cross-section of (1-inch) column in m^2
par.sl = 1; %cross-section of lateral connector in m^2
par.L = 0; % length of lateral connector
par.H = 74.7e-2;
xbars = linspace(0,73.5,101)*1e-2;
delxys = [44.6]*1e-2;
par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in Pa
numPar.x0 = 1e-5;
numPar.v0= 0; %initial velocity in m/s
numPar.tf = 5;
delta = 1e-6;%bar
IG_frequencies = zeros(length(delxys),length(xbars));
Steam_frequencies = zeros(length(delxys),length(xbars));
for j =1:length(delxys)
    par.delxy = delxys(j);
    for k = 1:length(xbars)
        par.xbar = xbars(k);
        par.ybar = xbars(k)+par.delxy;
        par.Pb0= par.rho*par.g*par.delxy +par.Pa0;
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

    Steam_frequencies(j,k) = frequency_calculator(par);
    fprintf('The Steam solution oscillates at a frequency of %f. \n',Steam_frequencies(j,k));

    IG_frequencies(j,k) = coldfreq(par);
    fprintf('The Ideal Gas solution oscillates at a frequency of %f. \n',IG_frequencies(j,k));
    end
end
%%post processing
figure;
plot(xbars*1e2, Steam_frequencies, xbars*1e2, IG_frequencies)
title('Frequency vs x-bar')
legend('Steam','Ideal Gas')
xlabel('x-bar (cm)')
ylabel('Frequency (Hz)')
IG_approximation_error = abs((Steam_frequencies-IG_frequencies)./IG_frequencies)*100;
figure;
plot(xbars*1e2, IG_approximation_error)
title('Steam vs Ideal Gas Frequencies')
xlabel('xbar (cm)');
ylabel('Percent Difference (%)');
KIG = 1./(par.gamma*P);
Ksteam = -dv_dp./vV;
figure;
plot(P*1e5,Ksteam*1e5,P*1e5,KIG)
title('Compressibility of Ideal Gas and Steam')
legend('Steam', 'Ideal Gas')
xlabel('Pressure (Pascal)')
ylabel('Compressibility (dP/dV)')