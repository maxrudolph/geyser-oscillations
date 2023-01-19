function [bracket,top,bottom] = bracket_val(par,model)
sbsc= par.sb ./ par.sc;
sbsl= par.sb ./ par.sl;
if model ==1
    k0 = -par.Vol_0.*par.FdPdv(par.Vol_0);
    xi0 = par.Fdudp(par.Vol_0);
    delta = 1e-3;
    dxidx0 = (par.Fdudp(par.sb*(par.H-par.xbar-delta))-par.Fdudp(par.sb*(par.H-par.xbar+delta) )) /(2*delta);
    %dkdx = dp/dvdx * v + dp/dv * dv/dx
    dpdvdx0v0  = -par.Vol_0*(par.FdPdv(par.sb*(par.H-par.xbar-delta))-par.FdPdv(par.sb*(par.H-par.xbar+delta) ) ) /(2*delta);
    dpdv0dvdx0 = -par.FdPdv(par.Vol_0)*par.sb;
    dkdx0 = -(dpdvdx0v0+dpdv0dvdx0);
    Hx = par.H-par.xbar;
    bottom = par.xbar+sbsc*par.ybar+sbsl*par.L;
    B = par.g.*(1+sbsc);
    C = 1/(par.rho*par.sb);
    D=(-xi0*k0)/(Hx^2);
    E=(dxidx0*k0)/Hx;
    G =(-xi0*dkdx0)/Hx;
    top=C*(D+E+G)+B;
    bracket = D+E+G;
end
if model==0
    sbsc= par.sb ./ par.sc;
    sbsl= par.sb ./ par.sl;
    bottom = par.xbar+sbsc*par.ybar+sbsl*par.L;
    top=1;
    bracket=1;
end
end