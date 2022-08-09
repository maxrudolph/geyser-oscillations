function [ts,xs,vs] = solve_15s(model,par,numPar)
%initialize x and v list
t=[0, numPar.tf]; %time period of solution in seconds
j=length(t);
x= zeros(1,j);
v=zeros(1,j);
z = [x;v];
%% Initial Conditions
assert(numPar.x0<(par.H-par.xbar)) %check that water cannot overflow
assert( (par.ybar - par.sb/par.sc*abs(numPar.x0)) > 0.0) % assert that water stays above lateral connector on conduit side.

switch model
    case 0
        vdot = @(t,z) steam_RHS(t,z,par);
    case 1
        vdot = @(t,z) IG_RHS(t,z,par);
end
%set initial conditions
z0= [numPar.x0;numPar.v0];
op=odeset('RelTol',1e-10,'AbsTol',1e-12);
%%Solving EOM
%pass parameters into vdot equation easily
[ts,ys] = ode15s(vdot,t,z0,op);
xs = ys(:,1);vs = ys(:,2);

end

