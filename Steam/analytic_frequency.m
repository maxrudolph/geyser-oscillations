function freq = analytic_frequency(par)
sb_sc= par.sb ./ par.sc;
sb_sl= par.sb ./ par.sl;

% define f(x) = dH/dP * -V*dP/dV
% f(0)/(H-xbar) + (f'(0)/(H-xbar) + f(0)/(H-xbar)^2)*x
tmp1 = par.Fdhdp(par.Vol_0) * -par.Vol_0*par.FdPdv(par.Vol_0);
delta = 1e-9;
fplus  = par.Fdhdp((par.Vol_0-par.sb*delta)) * -(par.Vol_0-par.sb*delta)*par.FdPdv((par.Vol_0-par.sb*delta));
fminus = par.Fdhdp((par.Vol_0+par.sb*delta)) * -(par.Vol_0+par.sb*delta)*par.FdPdv((par.Vol_0+par.sb*delta));
dfdx = (fplus-fminus)/(2*delta);
tmp2 = dfdx/(par.H-par.xbar) - tmp1/(par.H-par.xbar)^2;

% k0 = -par.Vol_0.*par.FdPdv(par.Vol_0/par.m);

% xi0 = par.Fdhdp(par.Vol_0/par.m);
% delta = 1e-3;
% dxidx0 = (par.Fdhdp((1/par.m)*par.sb*(par.H-par.xbar-delta))-par.Fdhdp((1/par.m)*par.sb*(par.H-par.xbar+delta) )) /(2*delta);
%dkdx = dp/dvdx * v + dp/dv * dv/dx
% dpdvdx0v0  = -par.Vol_0*(par.FdPdv((1/par.m)*par.sb*(par.H-par.xbar-delta))-par.FdPdv((1/par.m)*par.sb*(par.H-par.xbar+delta) ) ) /(2*delta);
% dpdv0dvdx0 = -par.FdPdv(par.Vol_0/par.m)*par.sb;
% dkdx0 = -(dpdvdx0v0+dpdv0dvdx0);
Hx = par.H-par.xbar;
A = par.xbar+sb_sc*par.ybar+sb_sl*par.L;
B = par.g.*(1+sb_sc);
C = -1/(par.rho*par.sb);
% D=(-xi0*k0)/(Hx^2);
% E=(dxidx0*k0)/Hx;
% G =(-xi0*dkdx0)/Hx;
% F=C*(D+E+G);
F = C*tmp2;

freq = sqrt((B+F)/A)/(2*pi);



end

