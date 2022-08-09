format long;
clc;clear;
close all;

addpath ../XSteam_Matlab_v2.6/;
reps = 5; %number of delta values to solve for
numerical_frequency = zeros(1,reps);
anfreq = zeros(1,reps);
maxfreq = zeros(1,reps);
fill = 0.9; %fraction of chamber above lateral connector that contains water
geom=0;
switch geom
    case 0
        fprintf('Using laboratory dimensions. \n');
        par.sb = 613e-4; % cross-section of bubble trap in m^2
        par.sc = 5.07e-4; % cross-section of (1-inch) column in m^2
        par.sl = 7.92e-4; %cross-section of lateral connector in m^2
        par.L = 7.37e-2; % length of lateral connector in m
        par.H = 37.15e-2; % height of bubble trap in m
        par.xbar = fill*par.H; %mean water height in bubble trap in m, must be in (0,H)
        delxys = linspace(0,par.H-par.xbar,reps);
        numPar.x0= -1e-9; %initial displacement in m
    case 1
        fprintf('Using Old Faithful''s dimensions. \n ')
        par.sb = 80; % cross-section of bubble trap in m^2
        par.sc = 1; % cross-section of geyser column in m^2
        par.H = 7; % height of bubble trap in m
        par.L = 2; % length of lateral connector in m
        par.sl = 7; %cross-section of lateral connector in m^2
        par.xbar = fill*par.H; %mean water height in bubble trap in m, must be in (0,H)
        numPar.x0= -par.sc/par.sb; %initial displacement in m
        delxys = linspace(0,par.H-par.xbar,reps);
end

par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in Pa
ybars = delxys+par.xbar; %mean water height in conduit in m, must be greater than Sb/Sc*x0
%par.delxy = par.ybar-par.xbar;
numPar.v0= 0; %initial velocity in m/s
numPar.h=0.00025; %time step (s)
numPar.tf = 3;
numPar.time=0: numPar.h : numPar.tf;
numPar.j=length(numPar.time);
figure;
for k = 1:length(delxys)
    %par.xbar = xbars(k);
    par.ybar = ybars(k);
    par.delxy = delxys(k);
    %% Make a lookup table with pre-computed thermodynamic properties.
    par.Vol_0 = par.sb*(par.H-par.xbar); %initial (equilibrium) volume in m^3
    P_0 = par.rho*par.g*(par.delxy)+par.Pa0; %initial pressure in (Pascal)
    par.m=par.Vol_0 / XSteam('vV_p', P_0/1e5); %vapor mass in kg
    
    xv = 1 - 1e-1;
    sv=XSteam('sV_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
    sl=XSteam('sL_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
    s = xv*sv + (1-xv)*sl;
    
    P=linspace((P_0-1e4)/1e5,(P_0+1e4)/1e5,51); %Range of pressure for lookup table (bar)
    vV=zeros(1,length(P));  % vapor volume (m^3)
    hH = zeros(1,length(P));% enthalpy (J/kg)
    dh_dp = zeros(1,length(P)); %dH/dP in (kJ/K/Pa)
    dv_dp = zeros(1,length(P)); %dV/dP in (m^3/Pa)
    for i=1:length(P)
        vV(i) = XSteam('v_ps',P(i),s/1e3)*par.m; %volume in m^3
        hH(i) = XSteam('h_ps',P(i),s/1e3)*1e3; %enthalpy in J*kg^-1
        xS(i) = XSteam('x_ps',P(i),s/1e3);
        delta = 1e-6;%bar
        dh_dp(i) = (XSteam('h_ps',P(i)+delta,s/1e3)-XSteam('h_ps',P(i)-delta,s/1e3))*1e3*par.m/(2*delta*1e5);% enthalpy in mks units, pressure in Pa
        dv_dp(i) = (XSteam('v_ps', P(i)+delta,s/1e3)-XSteam('v_ps', P(i)-delta,s/1e3))*par.m/(2*delta*1e5);% m^3/kg/Pa
    end
    %{
    figure()
    subplot(2,1,1);
    plot(P,dh_dp);
    xlabel('P ((bar)');
    ylabel('dh/dP (J/Pa)');
    subplot(2,1,2);
    plot(P,dh_dp.*(1./dv_dp));
    hold on
    plot(P_0/1e5*[1 1],get(gca,'YLim'));
    xlabel('P (bar)');
    ylabel('dh/dv');
    %}
    [~,i] = sort(vV);
    par.Fdhdp = griddedInterpolant(vV(i),dh_dp(i),'linear','none');
    par.FdPdv = griddedInterpolant(vV(i),1./dv_dp(i),'linear','none');
    
    [t1,x1,v1] = solve_15s(0,par,numPar);
    [x2,v2] = steam_solve(par,numPar);
    %%Plotting
    %figure;
    plot(t1,x1.*1e2);%numPar.time,x2*1e2);
    title('Solution x(t) of IVP')
    xlabel('Time (sec)'); grid on
    ylabel('X (cm)');hold on;
    %drawnow();
    % figure;
    % plot(x,v);
    % title("Phase Diagram of IVP");
    % xlabel('x');ylabel('dx/dt');
    %
    %
    % figure();
    % y = par.ybar - par.sb/par.sc*x;
    % plot(t,y.*1e2);
    % title('Solution y(t)');
    % xlabel('Time (sec)');
    % ylabel('Y (cm)')
    % %%compute v(t)
    %
    %
    %
    % vt = par.sb*(par.H - par.xbar - x)/par.m;
    % pt = interp1(vV,P,vt);
    % xt = zeros(size(pt));
    % for i = 1:length(pt)
    %     xt(i) = XSteam('x_ps',pt(i),s/1e3);
    % end
    % figure();
    % plot(t,xt);
    %
    % b = linspace(5e-2,37.15e-2,100);
    % figure; plot(b,par.FdPdv((1/par.m).*par.sb.*(par.H-par.xbar-b)));title('dpdv vs x')
    % figure; plot(b,par.Fdhdp((1/par.m).*par.sb.*(par.H-par.xbar-b)));title('dhdp vs x')
    %
    
    
    [~,locs1] = findpeaks(x1);
    if length(locs1) > 1
    oscillation_period=t1(locs1(end))-t1(locs1(end-1));
    numerical_frequency(k) = 1./oscillation_period;
    fprintf('The numerical solution oscillates at a frequency of %f. \n',numerical_frequency(k));
    
    
    anfreq(k) = frequency_calculator(par);
    maxfreq(k) = analytic_frequency(par);
    fprintf(['Uthkarsh''s analytical frequency is %f. \n' ...
        'Max''s analytical frequency is %f. \n'], anfreq(k), maxfreq(k));
    else
        anfreq(k) = NaN;
        maxfreq(k) = NaN;
    end
end
hold off;
%% post processing
figure();
plot(delxys*1e2, numerical_frequency);
hold on
plot(delxys*1e2,anfreq,'r--');
plot(delxys*1e2,maxfreq,'k--')
xlabel('ybar-xbar (cm)');
ylabel('Oscillation frequency (hz)');

% compute and plot error.
err = abs((maxfreq-numerical_frequency)./numerical_frequency)*100;
figure();
plot(delxys*1e2,err);
xlabel('ybar-xbar (m)');
ylabel('Percent Error (%)');