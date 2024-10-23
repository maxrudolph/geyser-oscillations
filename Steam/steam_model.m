format long;
clear;
close all;

addpath ./XSteam_Matlab_v2.6/;
addpath ..;
reps = 3; %number of delta values to solve for
numerical_frequency = zeros(1,reps);
anfreq = zeros(1,reps);
maxfreq = zeros(1,reps);
fill = 1-.05; %fraction of chamber above lateral connector that contains water
geom=1;
switch geom
    case 0
        fprintf('Using laboratory dimensions. \n');
        par.sb = 1852e-4; % cross-section of bubble trap in m^2
        par.sc = 5.07e-4; % cross-section of (1-inch) column in m^2
        par.sl = 1; %cross-section of lateral connector in m^2
        par.L = 0; % length of lateral connector in m
        par.H = 74.7e-2; % height of bubble trap in m
        par.xbar = fill*par.H; %mean water height in bubble trap in m, must be in (0,H)
        delxys = linspace(0,par.H-par.xbar,reps);
        x0= [5e-4,1e-3,2e-3]; %initial displacement in m
    case 1
        fprintf('Using Old Faithful''s dimensions. \n ')
        par.sb = 300; % cross-section of bubble trap in m^2
        par.sc = 1; % cross-section of geyser column in m^2
        par.H = 8; % height of bubble trap in m
        par.L = 0; % length of lateral connector in m
        par.sl = 0.1; %cross-section of lateral connector in m^2
        par.xbar = fill*par.H; %mean water height in bubble trap in m, must be in (0,H)
        % par.x0= -par.sc/par.sb; %initial displacement in m
        delxys = [25 25 25];%linspace(0,20,reps);
        x0=-[1e-4,5e-4,1e-3,1e-2];
end

par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in Pa
ybars = delxys+par.xbar; %mean water height in conduit in m, must be greater than Sb/Sc*x0
numPar.v0= 0; %initial velocity in m/s
numPar.tf = 30;
figure(101); clf;

for k = 1:reps
    par.ybar = ybars(k);
    par.delxy = delxys(k);
    numPar.x0 = x0(k);
    %% Make a lookup table with pre-computed thermodynamic properties.
    par.Vol_0 = par.sb*(par.H-par.xbar); %initial (equilibrium) volume in m^3
    P_0 = par.rho*par.g*(par.delxy)+par.Pa0; %initial pressure in (Pascal)
    par.m=par.Vol_0 / XSteam('vV_p', P_0/1e5); %vapor mass in kg

    xv = 1 - 5e-1;
    sv=XSteam('sV_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
    sl=XSteam('sL_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
    s = xv*sv + (1-xv)*sl; % gives specific entropy of water+vapor mixture in J/kg/C
    V0 = XSteam('V_ps',P_0/1e5,s/1e3); % volume per unit mass
    par.m = par.Vol_0 / V0;

    P=linspace(1.01e5/1e5,(P_0+10e5)/1e5,1001); %Range of pressure for lookup table (bar)
    vV=zeros(1,length(P));  % vapor volume (m^3)
    du_dp = zeros(1,length(P)); %dU/dP in (kJ/K/Pa)
    dv_dp = zeros(1,length(P)); %dV/dP in (m^3/Pa)
    xS = zeros(1,length(P)); % vapor fraction
    T = zeros(1,length(P)); % Temperature
    for i=1:length(P)
        vV(i) = XSteam('v_ps',P(i),s/1e3)*par.m; %volume in m^3 (specific volume x mass)
        xS(i) = XSteam('x_ps',P(i),s/1e3); %vapor fraction (unitless)
        T(i) = XSteam('Tsat_p',P(i)); % temperature in degree C
        delta = sqrt(eps(P(i)));%1e-8;%bar
        du_dp(i) = (XSteam('u_ps',P(i)+delta,s/1e3)-XSteam('u_ps',P(i)-delta,s/1e3))*1e3*par.m/(2*delta*1e5);% enthalpy in mks units, pressure in Pa
        dv_dp(i) = (XSteam('v_ps',P(i)+delta,s/1e3)-XSteam('v_ps',P(i)-delta,s/1e3))*par.m/(2*delta*1e5);% m^3/Pa
    end
    
    [~,i] = sort(vV);
    par.Fdudp = griddedInterpolant(vV(i),du_dp(i),'linear','none');
    par.FdPdv = griddedInterpolant(vV(i),1./dv_dp(i),'linear','none');
    
    [t1,x1,v1] = solve_15s(0,par,numPar);

    % reconstruct pressure, temperature for the time series
    V1 = (par.H-(par.xbar+x1))*par.sb;
    P1 = interp1(vV(i),P(i),V1);
    T1 = interp1(vV(i),T(i),V1);
    X1 = interp1(vV(i),xS(i),V1);

    %Plotting
    %set(gca,'BackgroundColor,','None')
    figure(101);
    plot(t1,x1,'.');
    title('Solution x(t) of IVP (m)')
    xlabel('Time (sec)'); grid on
    ylabel('X (m)');hold on;
    fontsize(gcf,24,"points");

    figure(102);
    subplot(4,1,1);
    plot(t1,-par.sb/par.sc*x1,'.');
    title('Solution y(t) of IVP (m)')
    xlabel('Time (sec)'); grid on
    ylabel('Y (m)');hold on;
    fontsize(gcf,24,"points");
    subplot(4,1,2);
    plot(t1,P1);
    hold on
    ylabel('Pressure (bar)')
    % plot(t1,)
    subplot(4,1,3);
    plot(t1,T1);
    hold on
    ylabel('T (C)')
    subplot(4,1,4);
    plot(t1,X1);
    hold on
    ylabel('X')
%     figure();
%     plot(x1,v1);
%     title("Phase Diagram of IVP");
%     xlabel('x');ylabel('dx/dt');
    %drawnow();

    %figure();
    %y = par.ybar - par.sb/par.sc*x1;
    %plot(t1,y.*1e2);
    %title('Solution y(t)');
    %xlabel('Time (sec)');
    %ylabel('Y (cm)')
    %%compute v(t)



    %vt = par.sb*(par.H - par.xbar - x1)/par.m;
    %pt = interp1(vV,P,vt);
    %xt = zeros(size(pt));
    %for i = 1:length(pt)
    %    xt(i) = XSteam('x_ps',pt(i),s/1e3);
    %end
    %figure();
    %plot(t1,xt);

    %b = linspace(5e-2,37.15e-2,100);
    %figure; plot(b,par.FdPdv((1/par.m).*par.sb.*(par.H-par.xbar-b)));title('dpdv vs x')
    %figure; plot(b,par.Fdudp((1/par.m).*par.sb.*(par.H-par.xbar-b)));title('dhdp vs x')



    [~,locs1] = findpeaks(x1);
    if length(locs1) > 1
        oscillation_period=t1(locs1(end))-t1(locs1(end-1));
        numerical_frequency(k) = 1./oscillation_period;
        fprintf('The numerical solution oscillates at a frequency of %f. \n',numerical_frequency(k));


        anfreq(k) = frequency_calculator(par);
        fprintf('Uthkarsh''s analytical frequency is %f. \n', anfreq(k));
    else
        anfreq(k) = NaN;
    end
end
hold off;
%% post processing
% figure();
% plot(delxys*1e2, numerical_frequency);
% hold on
% plot(delxys*1e2,anfreq,'r--');
% xlabel('ybar-xbar (cm)');
% ylabel('Oscillation frequency (Hz)');
% title('Comparison of Numerical and Analytical Frequencies')
% legend('Numerical','Analytical')

% compute and plot error.
% err = ((anfreq-numerical_frequency)./numerical_frequency)*100;
% figure();
% plot(delxys*1e2,err);
% xlabel('ybar-xbar (m)');
% ylabel('Percent Error (%)');
% figure;
% reference= linspace(min(anfreq)-0.01,max(anfreq)+0.01,3);
% plot(reference,reference,'-');hold on;
% plot(numerical_frequency,anfreq, '.', 'MarkerSize', 15);  hold off;
% xlabel('Numerical Frequency (Hz)');
% ylabel('Analytical Frequency (Hz)');
% title('Comparison of Numerical and Analytical Frequencies');
figure()
% subplot(2,1,1);
plot(T,P*1e2);
xlim([95 120])
ylim([0.8 2]*1e2)
xlabel('Temperature (Celsius)');
ylabel('Pressure (kPa)');
set(gca,'Ydir','reverse')
title('Water Phase Diagram')
% subplot(2,1,2);
% plot(P,du_dp.*(1./dv_dp));
% hold on
% plot(P_0/1e5*[1 1],get(gca,'YLim'));
% xlabel('P (bar)');
% ylabel('dh/dv');
