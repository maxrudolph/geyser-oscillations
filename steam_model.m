format long
clc
clear

par.delxy = 5*10^(-2); % vertical difference between xbar and ybar in m
par.sb = 613*10^(-4); % cross-section of bubble trap in m^2
par.sc = 5.07*10^(-4);% cross-section of (1-inch) column in m^2
par.sl = 7.92*10^(-4); %cross-section of lateral connector in m^2
par.L = 7.37*10^(-2); % length of lateral connector in m
par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.H = 37.15*10^(-2); % height of bubble trap in cm
par.alpha = 5/2;%diatomic ideal gas constant - unitless
par.Pa0 =1e5; %atmospheric pressure at equilibrium in Bar
%(rho=1)g(ybar-xbar) = Pb0-Pa0 from paper
fill = 0.5; %fraction of chamber above lateral connector that contains water
par.xbar = fill.*par.H;
par.ybar = par.delxy + par.xbar;


P=linspace(0.007,2,10000); %selected range of pressures in bar
vV=zeros(1,length(P));
hH = zeros(1,length(P));
dh_dp = zeros(1,length(P)); %enthalpy per bar
dv_dp = zeros(1,length(P)); %m3 per bar
Vol_0 = par.sb*(par.H-par.xbar); %initial volume in m^3
P_0 = par.rho*par.g*(par.delxy)+par.Pa0; %initial pressure in kg*s^-2*m^-1 (Pascal)
par.m=Vol_0 / XSteam('vV_p', (P_0*10^(-5))); %vapor mass in kg
s=XSteam('sV_p',(P_0*10^(-5))); %specific entropy kJ*kg^-1*degC^-1
for i=1:length(P)
    vV(i) = XSteam('v_ps', P(i),s); %specific volume in m^3*kg^-1
    hH(i) = XSteam('h_ps',P(i),s); %enthalpy in kJ*kg^-1
    delta = 1e-6;%bar
    dh_dp(i) = (XSteam('h_ps',P(i)+delta,s)-XSteam('h_ps',P(i)-delta,s))/(2*delta)/1e5;
    dv_dp(i) = (XSteam('v_ps', P(i)+delta,s)-XSteam('v_ps', P(i)-delta,s))/(2*delta)/1e5;
end

[~,i] = sort(vV);
par.Fdhdp = griddedInterpolant(vV(i),dh_dp(i));
par.FdPdv = griddedInterpolant(vV(i),1./dv_dp(i));

h=0.001; %time step
t=0:h:15; %time period of solution in seconds

x0= 0.001*1e-2; %initial displacement in cm
v0= 0; %initial velocity in cm*s^-1

assert(x0<(par.H-par.xbar)) %sense check s.t. water cannot overflow

%initialize x and v list
x= zeros(1,length(t));
v=zeros(1,length(t));
tdiff = zeros(1,length(t));
%set initial conditions
x(1)= x0;
v(1)= v0;

%pass parameters into vdot equation easily
vdot = @(t,x,v) vdot1(t,x,v,par);

%running RK4, must compute X and V slopes simultaneously
for i=1:length(t)-1
    k1x = xdot(t(i),x(i),v(i));
    k1v = vdot(t(i),x(i),v(i));
    k2x = xdot((t(i)+(.5*h)), (x(i) + .5*h*k1x),(v(i) + .5*h*k1v));
    k2v = vdot((t(i)+(.5*h)), (x(i) + .5*h*k1x),(v(i) + .5*h*k1v));
    k3x = xdot((t(i)+(.5*h)), (x(i) + .5*h*k2x),(v(i) + .5*h*k2v));
    k3v = vdot((t(i)+(.5*h)), (x(i) + .5*h*k2x),(v(i) + .5*h*k2v));
    k4x = xdot((t(i)+h), (x(i)+h*k3x), (v(i)+h*k3v));
    k4v = vdot((t(i)+h), (x(i)+h*k3x), (v(i)+h*k3v));
    x(i+1) = x(i) + (h/6)*(k1x+2*k2x+2*k3x+k4x);
    v(i+1) = v(i) + (h/6)*(k1v+2*k2v+2*k3v+k4v);
end

figure;
plot(t,x);
title('Solution x(t) of IVP')
xlabel('t'); grid on

figure;
plot(x,v);
title("Phase Diagram of IVP")
xlabel('x');ylabel('dx/dt');

function dxdt= xdot(t,x,v)
    dxdt = v;
end

function dvdt = vdot1(t,x,v,par)
    %precompute coefficients to visually declutter
    %equation 9 will look something like:  A*vdot = Bv^2 + C* O(x)
    sbsc= par.sb ./ par.sc;
    sbsl= par.sb ./ par.sl;
    A = par.xbar+x+(sbsc .* (par.ybar - (sbsc*x))) + sbsl*par.L ;
    B = (-0.5)*(1-(sbsc)^2);
    C = -par.g * (par.xbar+x);
    D = par.g * (par.ybar - (sbsc*x));
    agp = par.alpha * (1- par.gamma)*(1 / par.rho);
    vol = (1/par.m).*par.sb.*(par.H-par.xbar-x);
    E = -(par.m/par.rho).*par.Fdhdp(vol)*par.FdPdv(vol);
    F = -agp * par.Pa0*10^-5;
    
    dvdt = (1./A)*(B * v^2 +C+D+E+F);
end