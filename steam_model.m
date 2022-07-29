format long;clc;clear;addpath ./XSteam_Matlab_v2.6/;
P=linspace(0.007,2,5000); %Range of pressure for lookup table (bar)
vV=zeros(1,length(P));  % vapor volume (m^3/kg)
hH = zeros(1,length(P));% enthalpy (J/kg)
dh_dp = zeros(1,length(P)); %dH/dP in (kJ/K/Pa)
dv_dp = zeros(1,length(P)); %dV/dP in (m^3/Pa)
par.sb = 0.0613; % cross-section of bubble trap in m^2
par.sc = 5.07e-4;% cross-section of (1-inch) column in m^2
par.sl = 7.92e-4; %cross-section of lateral connector in m^2
par.L = 7.37e-2; % length of lateral connector in m
par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.H = 37.15e-2; % height of bubble trap in m
par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in Pa
%fill = 0.9; %fraction of chamber above lateral connector that contains water
par.xbar = 30e-2; %mean water height in bubble trap in m, must be in (0,H)
par.ybar=15e-2; %mean water height in conduit in m, must be greater than Sb/Sc*x0
par.delxy = par.ybar-par.xbar;

numPar.x0= 0.0005*1e-2; %initial displacement in m
numPar.v0= 0; %initial velocity in m/s
numPar.h=0.001; %time step (s)
numPar.tf = 3;
numPar.time=0: numPar.h : numPar.tf;
numPar.j=length(numPar.time);

%% Make a lookup table with pre-computed thermodynamic properties.
par.Vol_0 = par.sb*(par.H-par.xbar); %initial (equilibrium) volume in m^3
P_0 = par.rho*par.g*(par.delxy)+par.Pa0; %initial pressure in (Pascal)
par.m=par.Vol_0 / XSteam('vV_p', (P_0*10^(-5))); %vapor mass in kg
s=XSteam('sV_p',(P_0*10^(-5)))*1000; %specific entropy J*kg^-1*degC^-1
for i=1:length(P)
    vV(i) = XSteam('v_ps',P(i),s/1e3); %specific volume in m^3*kg^-1
    hH(i) = XSteam('h_ps',P(i),s/1e3)*1000; %enthalpy in J*kg^-1
    delta = 1e-6;%bar
    dh_dp(i) = (XSteam('h_ps',P(i)+delta,s/1e3)-XSteam('h_ps',P(i)-delta,s/1e3))*1000/(2*delta)/1e5;% enthalpy in mks units, pressure in Pa
    dv_dp(i) = (XSteam('v_ps', P(i)+delta,s/1e3)-XSteam('v_ps', P(i)-delta,s/1e3))*1000/(2*delta)/1e5;
end

[~,i] = sort(vV);
par.Fdhdp = griddedInterpolant(vV(i),dh_dp(i));
par.FdPdv = griddedInterpolant(vV(i),1./dv_dp(i));

[x,v] = steam_solve(par,numPar);
%%Plotting
%{
b = linspace(5e-2,37.15e-2,100);
figure; plot(b,par.Fdhdp((1/par.m).*par.sb.*(par.H-par.xbar-b)));title('dhdp vs x')

t = numPar.time;
figure;
plot(t,x.*1e2);
title('Solution x(t) of IVP')
xlabel('Time (sec)'); grid on
ylabel('X (cm)');

figure;
plot(x,v);
title("Phase Diagram of IVP");
xlabel('x');ylabel('dx/dt');


figure();
y = par.ybar - par.sb/par.sc*x;
plot(t,y.*1e2);
title('Solution y(t)');
xlabel('Time (sec)');
ylabel('Y (cm)')
%%compute v(t)


figure; plot(b,par.FdPdv((1/par.m).*par.sb.*(par.H-par.xbar-b)));title('dpdv vs x')
vt = par.sb*(par.H - par.xbar - x)/par.m;
pt = interp1(vV,P,vt);
xt = zeros(size(pt));
for i = 1:length(pt)
    xt(i) = XSteam('x_ps',pt(i),s/1e3);
end
figure();
plot(t,xt);
%}

[~,locs1] = findpeaks(x);
T=(locs1(2)-locs1(1))*numPar.h;
nmfreq = 1/T;
fprintf('The numerical solution oscillates at a frequency of %f. \n',nmfreq);


anfreq = frequency_calculator(par);
fprintf('The analytical frequency is %f. \n', anfreq);
