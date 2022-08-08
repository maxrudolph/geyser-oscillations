format long;close all;clc;clear;

reps=5; %number of delta values
nmfreq = zeros(1,reps);
anfreq = zeros(1,reps);

par.sb = 613; % cross-section of bubble trap in cm^2
par.sc = 5.07;% cross-section of (1-inch) column in cm^2
par.sl = 7.92; %cross-section of lateral connector in cm^2
par.L = 7.37; % length of lateral connector in cm
par.g = 981;% gravitational acceleration in cm*s^-2
par.rho=1; % density of water in gm*cm^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.H = 37.15; % height of bubble trap in cm
par.alpha = 5/2;%diatomic ideal gas constant - unitless
par.Pa0 =10^6; %atmospheric pressure at equilibrium in cm^−1*g*s^−2
%(rho=1)g(ybar-xbar) = Pb0-Pa0 from paper
par.xbar = 25; %mean water height in bubble trap in cm, must be within (0,H)
delxys = linspace(0,par.H-par.xbar,reps); %mean water height difference (ybar-xbar) in m
ybars = delxys+par.xbar;

numPar.h=0.005; %time step
numPar.tf=3; %time period of solution in seconds

numPar.x0= 1e-6; %initial displacement in cm
numPar.v0= 0; %initial velocity in cm*s^-1
numPar.time=0: numPar.h : numPar.tf;
numPar.j=length(numPar.time);
figure;
for k = 1:length(delxys)
    par.ybar = ybars(k);
    par.delxy = delxys(k);
    par.Pb0= par.rho*par.g*par.delxy +par.Pa0; %bubble trap pressure at equilibrium in cm^−1*g*s^−2
    [x,v] = coldsolve(par,numPar);
    t=numPar.time;
    
    plot(t,x);
    hold on;
    
    [~,locs1] = findpeaks(x);
    T=(locs1(2)-locs1(1))*numPar.h;
    nmfreq(k) = 1/T;
    fprintf('The numerical solution oscillates at a frequency of %f. \n',nmfreq(k));

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
figure;plot(delxys, nmfreq,delxys,anfreq);
err = abs((anfreq-nmfreq)./nmfreq);
figure;plot(delxys,err);