function freq = coldfreq(par)
sbsc= par.sb ./ par.sc;
sbsl= par.sb ./ par.sl;
A = par.xbar+sbsc*par.ybar+sbsl*par.L;
B = par.g.*(1+sbsc);
C = par.Pb0*par.gamma*(1/(par.rho*(par.H-par.xbar)));
freq = sqrt((B+C)/A)/(2*pi);
end