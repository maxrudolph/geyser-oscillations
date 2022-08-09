function dvdt = IG_RHS(~,z,par)
x=z(1);
v=z(2);

sb_sc= par.sb ./ par.sc;
sb_sl= par.sb ./ par.sl;
A = par.xbar+x+(sb_sc .* (par.ybar - (sb_sc*x))) + sb_sl*par.L ;
B = (-0.5)*(1-(sb_sc)^2);
C = -par.g * (par.xbar+x);
D = par.g * (par.ybar - (sb_sc*x));
agp = par.alpha * (1- par.gamma)*(1 / par.rho);
E = agp * par.Pb0 * ( (par.H - par.xbar) / (par.H - par.xbar - x) )^(par.gamma);
F = -agp * par.Pa0;

dvdt = [v; (1./A)*(B * v^2 +C+D+E+F)];
end