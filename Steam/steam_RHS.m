function dvdt = steam_RHS(~,z,par)
x=z(1);
v=z(2);
%precompute coefficients to visually declutter
%equation 9 will look something like:  A*vdot = Bv^2 + C* O(x)
sb_sc= par.sb / par.sc;
sb_sl= par.sb / par.sl;
A = par.xbar+x+ sb_sc*(par.ybar - sb_sc*x) + sb_sl*par.L ;
assert( par.ybar - sb_sc*x > 0 );% assert that the liquid level in the conduit is positive

B = (-0.5)*(1-(sb_sc)^2);
C = -par.g * (par.xbar+x);
D = par.g * (par.ybar - (sb_sc*x));
agp = par.alpha * (1- par.gamma)*(1 / par.rho);

vol = par.sb*(par.H-par.xbar-x);
assert(vol>0); % assert that the bubble trap volume is positive


E = (1/par.rho).*par.Fdhdp(vol)*par.FdPdv(vol);
F = -agp * par.Pa0;

dvdt =[v; (1./A)*(B*v^2 +C+D+E+F)];
end