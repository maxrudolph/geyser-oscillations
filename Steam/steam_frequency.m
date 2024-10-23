function [frequency,uncertainty] = steam_frequency(par)
% this function is meant to be a wrapper for frequency_calculator
% it includes the construction of a lookup table to compute steam
% compressibility from the steam tables.
% par is a data structure with the following fields:
% sb
% sc
% H
% xbar
% delxy
% g
% rho
% Pa0

delta = 1e-3; % step size for finite difference calculation.

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

[frequency,uncertainty] = frequency_calculator(par);
