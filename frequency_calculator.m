function freq = frequency_calculator(par)
sbsc= par.sb ./ par.sc;
sbsl= par.sb ./ par.sl;
k0 = -par.Vol_0.*par.FdPdv(par.Vol_0/par.m);

xi0 = par.Fdhdp(par.Vol_0/par.m);
delta = 1e-3;
dxidx0 = (par.Fdhdp((1/par.m)*par.sb*(par.H-par.xbar-delta))-par.Fdhdp((1/par.m)*par.sb*(par.H-par.xbar+delta) )) /(2*delta);
%dkdx = dp/dvdx * v + dp/dv * dv/dx
dpdvdx0v0  = -par.Vol_0*(par.FdPdv((1/par.m)*par.sb*(par.H-par.xbar-delta))-par.FdPdv((1/par.m)*par.sb*(par.H-par.xbar+delta) ) ) /(2*delta);
dpdv0dvdx0 = -par.FdPdv(par.Vol_0/par.m)*par.sb;
dkdx0 = -(dpdvdx0v0+dpdv0dvdx0);
Hx = par.H-par.xbar;
A = par.xbar+sbsc*par.ybar+sbsl*par.L;
B = par.g.*(1+sbsc);
C = -1/(par.rho*par.sb);
D=(-xi0*k0)/(Hx^2);
E=(dxidx0*k0)/Hx;
G =(-xi0*dkdx0)/Hx;
F=C*(D+E+G);
freq = sqrt((B+F)/A)/(2*pi);
end