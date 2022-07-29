format long
%clc
%clear
addpath ./XSteam_Matlab_v2.6/


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
par.xbar = 20e-2; %mean water height in bubble trap in m, must be in (0,H)
par.ybar=15e-2; %mean water height in conduit in m, must be greater than Sb/Sc*x0
par.delxy = par.ybar-par.xbar; %height difference between xbar and ybar in m

sbsc= par.sb ./ par.sc;
sbsl= par.sb ./ par.sl;

P=linspace(0.007,2,5000); %Range of pressure for lookup table (bar)
vV=zeros(1,length(P));  % vapor volume (m^3/kg)
hH = zeros(1,length(P));% enthalpy (J/kg)
dh_dp = zeros(1,length(P)); %dH/dP in (kJ/K/Pa)
dv_dp = zeros(1,length(P)); %dV/dP in (m^3/Pa)
Vol_0 = par.sb*(par.H-par.xbar); %initial (equilibrium) volume in m^3
P_0 = par.rho*par.g*(par.delxy)+par.Pa0; %initial pressure in (Pascal)
par.m=Vol_0 / XSteam('vV_p', (P_0*10^(-5))); %vapor mass in kg
s=XSteam('sV_p',(P_0*10^(-5)))*1000; %specific entropy J*kg^-1*degC^-1
delta = 1e-6;%bar
for i=1:length(P)
    vV(i) = XSteam('v_ps',P(i),s/1e3); %specific volume in m^3*kg^-1
    hH(i) = XSteam('h_ps',P(i),s/1e3)*1000; %enthalpy in J*kg^-1
    dh_dp(i) = (XSteam('h_ps',P(i)+delta,s/1e3)-XSteam('h_ps',P(i)-delta,s/1e3))*1000/(2*delta)/1e5;% enthalpy in mks units, pressure in Pa
    dv_dp(i) = (XSteam('v_ps', P(i)+delta,s/1e3)-XSteam('v_ps', P(i)-delta,s/1e3))*1000/(2*delta)/1e5;
end

[~,i] = sort(vV);
par.Fdhdp = griddedInterpolant(vV(i),dh_dp(i));
par.FdPdv = griddedInterpolant(vV(i),1./dv_dp(i));

k0 = -Vol_0.*par.FdPdv(Vol_0/par.m);
xi0 = par.Fdhdp(Vol_0/par.m);
delta = 1e-3;
dxidx0 = (par.Fdhdp((1/par.m)*par.sb*(par.H-par.xbar-delta))-par.Fdhdp((1/par.m)*par.sb*(par.H-par.xbar+delta) )) /(2*delta);

A = par.xbar+sbsc*par.ybar+sbsl*par.L;
B = par.g.*(1+sbsc);
C = (1/(par.rho*par.sb))*((-xi0*k0)/((par.H-par.xbar)^2) + (dxidx0*k0)/(par.H-par.xbar));

freq = 2*pi/(sqrt((B-C)/A))