format long;close all;clc;clear;
addpath ..
reps=10; %number of delta values
numerical_frequency = zeros(1,reps);
anfreq = zeros(1,reps);
fill = 0.8; %fraction of chamber that contains water
par.sb = 1852e-4; % cross-section of bubble trap in cm^2
par.sc = 20.26e-4;% cross-section of (1-inch) column in cm^2
par.sl = 1; %cross-section of lateral connector in cm^2
par.L = 0; % length of lateral connector in cm
par.g = 9.81;% gravitational acceleration in cm*s^-2
par.rho=1000; % density of water in gm*cm^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.H = 74.7e-2; % height of bubble trap in cm
par.alpha = 5/2;%diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in cm^−1*g*s^−2
%(rho=1)g(ybar-xbar) = Pb0-Pa0 from paper
par.xbar = fill*par.H; %mean water height in bubble trap in cm, must be within (0,H)
delxys = linspace(0,par.H-par.xbar,reps); %mean water height difference (ybar-xbar) in m
ybars = delxys+par.xbar;

numPar.h=0.001; %time step
numPar.tf=3; %time period of solution in seconds

numPar.x0= par.sc/par.sb*1e-2; %initial displacement in m
numPar.v0= 0; %initial velocity in m/s
numPar.time=0: numPar.h : numPar.tf;
numPar.j=length(numPar.time);
figure;
for k = 1:length(delxys)
    par.ybar = ybars(k);
    par.delxy = delxys(k);
    par.Pb0= par.rho*par.g*par.delxy +par.Pa0; %bubble trap pressure at equilibrium in cm^−1*g*s^−2
    [t1,x1,v1] = solve_15s(1,par,numPar);
    [x2,v2] = coldsolve(par,numPar);
    t=numPar.time;
    
    plot(t1,x1);
    hold on;
    
    [~,locs1] = findpeaks(x1);
    if length(locs1) > 1
    oscillation_period=t1(locs1(end))-t1(locs1(end-1));
    numerical_frequency(k) = 1/oscillation_period;
    fprintf('The numerical solution oscillates at a frequency of %f. \n',numerical_frequency(k));
    end
    anfreq(k) = coldfreq(par);
    fprintf('The analytical frequency is %f. \n', anfreq(k));
    % figure;
    % plot(x,v);
    % title("Phase Diagram of IVP")
    % xlabel('x');ylabel('dx/dt');
    % for i = 1:length(t)-1
    %     tdiff(i) = t(i+1)-t(i);
    % end
    % tdiff(length(t)) = tdiff(length(t)-1);
    % A = [t',x',tdiff'];
    % writematrix(A, 'eq9rk4.txt','Delimiter',' ');
end
title('Solution x(t) of IVP')
xlabel('t'); grid on;hold off;
figure;
plot(delxys*1e2, numerical_frequency,delxys*1e2,anfreq);
title('Frequency vs ybar-xbar');
xlabel('ybar-xbar (cm)');
ylabel('Frequency (Hz)');
legend('Numerical', 'Analytical')
err = abs((anfreq-numerical_frequency)./numerical_frequency)*100;
figure;
plot(delxys*1e2,err);
title('Frequency - Percent Error');
xlabel('ybar-xbar (cm)');
ylabel('Error (%)');