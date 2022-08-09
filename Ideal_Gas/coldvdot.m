function dvdt = coldvdot(~,x,v,par)
    
    %precompute coefficients to visually declutter
    %equation 9 will look something like:  A*vdot = Bv^2 + C* O(x)
    sbsc= par.sb ./ par.sc;
    sbsl= par.sb ./ par.sl;
    A = par.xbar+x+(sbsc .* (par.ybar - (sbsc*x))) + sbsl*par.L ;
    B = (-0.5)*(1-(sbsc)^2);
    C = -par.g * (par.xbar+x);
    D = par.g * (par.ybar - (sbsc*x));
    agp = par.alpha * (1- par.gamma)*(1 / par.rho);
    E = agp * par.Pb0 * ( (par.H - par.xbar) / (par.H - par.xbar - x) )^(par.gamma);
    F = -agp * par.Pa0;
    
    dvdt = (1./A)*(B * v^2 +C+D+E+F);
end