function [x,v] = steam_solve(par,numPar)
%initialize x and v list
t=numPar.time; %time period of solution in seconds
h = numPar.h;
x= zeros(1,numPar.j);
v=zeros(1,numPar.j);
%% Set up parameters


%% Initial Conditions
assert(numPar.x0<(par.H-par.xbar)) %check that water cannot overflow
assert( (par.ybar - par.sb/par.sc*abs(numPar.x0)) > 0.0) % assert that water stays above lateral connector on conduit side.


%set initial conditions
x(1)= numPar.x0;
v(1)= numPar.v0;

%%Solving EOM
%pass parameters into vdot equation easily
vdot = @(t,x,v) vdot_fun(t,x,v,par);

%running RK4, must compute X and V slopes simultaneously
for i=1:numPar.j-1
    k1x = xdot_fun(t(i),x(i),v(i));
    k1v = vdot(t(i),x(i),v(i));
    k2x = xdot_fun((t(i)+(.5*h)), (x(i) + .5*h*k1x),(v(i) + .5*h*k1v));
    k2v = vdot((t(i)+(.5*h)), (x(i) + .5*h*k1x),(v(i) + .5*h*k1v));
    k3x = xdot_fun((t(i)+(.5*h)), (x(i) + .5*h*k2x),(v(i) + .5*h*k2v));
    k3v = vdot((t(i)+(.5*h)), (x(i) + .5*h*k2x),(v(i) + .5*h*k2v));
    k4x = xdot_fun((t(i)+h), (x(i)+h*k3x), (v(i)+h*k3v));
    k4v = vdot((t(i)+h), (x(i)+h*k3x), (v(i)+h*k3v));
    x(i+1) = x(i) + (h/6)*(k1x+2*k2x+2*k3x+k4x);
    v(i+1) = v(i) + (h/6)*(k1v+2*k2v+2*k3v+k4v);
end
end