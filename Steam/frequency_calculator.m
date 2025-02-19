function [freq,propagated_error] = frequency_calculator(par)
% calculate the predicted oscillation frequency and its uncertainty.
% par is a data structure with the following
% sb - cross section area of bubble trap
% sc - cross section area of conduit
% H - height of tank above pickup tube
% xbar - equilibrium level of water in bubble trap
% ybar - equilibrium level of water in conduit

sbsc= par.sb ./ par.sc;
sbsl= par.sb ./ par.sl;
Hx = par.H-par.xbar;
% sig_area = 5e-4; %uncertainty in measures of area, 5 cm^2
sig_sb = 2*par.sb*5e-3/par.sb; % assume 5 mm uncertainty in radius dA/A = 2*dr/r
sig_sc = 2*par.sc*1e-4/par.sc; % assume 0.1 mm uncertainty in radius
sig_len = 5e-3; %uncertainty in measures of length, 0.5 cm (fill level quite uncertain)

k0 = -par.Vol_0.*par.FdPdv(par.Vol_0);

xi0 = par.Fdudp(par.Vol_0);
delta = sqrt(eps( (par.H-par.xbar)*par.sb) );
dxidx0 = (par.Fdudp(par.sb*(par.H-par.xbar-delta))-par.Fdudp(par.sb*(par.H-par.xbar+delta) )) /(2*delta);
%dkdx = dp/dvdx * v + dp/dv * dv/dx
dpdvdx0v0  = -par.Vol_0*(par.FdPdv(par.sb*(par.H-par.xbar-delta))-par.FdPdv(par.sb*(par.H-par.xbar+delta) ) ) /(2*delta);
dpdv0dvdx0 = -par.FdPdv(par.Vol_0)*par.sb;
dkdx0 = -(dpdvdx0v0+dpdv0dvdx0);


A = par.xbar+sbsc*par.ybar+sbsl*par.L;
B = par.g.*(1+sbsc);
C = 1/(par.rho*par.sb);
D=(-xi0*k0);
E=(dxidx0*k0);
G =(-xi0*dkdx0);
Bracket=(G/Hx+D/Hx^2+E/Hx);
F=C*Bracket;
angfreq2 = (B+F)/A;
angfreq =sqrt(angfreq2);
freq =angfreq/(2*pi);

% sig_SbSc = par.sb/par.sc *sqrt((sig_area/par.sb)^2 +(sig_area/par.sc)^2);
sig_SbSc = par.sb/par.sc *sqrt((sig_sb/par.sb)^2 +(sig_sc/par.sc)^2);
sig_gbc = par.g*sig_SbSc;
sig_bcy = par.ybar*(par.sb/par.sc)*sqrt( (sig_SbSc/(par.sb/par.sc))^2 + sig_len^2*par.ybar^-2 );
sig_A = sqrt(sig_len^2+sig_bcy^2);
sig_hx = sqrt(2*sig_len^2);
sig_invhx = sig_hx/(Hx)^2;
sig_invhx2 = sqrt(8)*sig_len/Hx^3;
sig_C = sig_sb/(par.rho*par.sb^2);
sig_bracket = sqrt(2*sig_len^2*G^2*Hx^-4+8*sig_len^2*D^2*Hx^-6+2*sig_len^2*E^2*Hx^-4);
sig_F = abs(F)*sqrt((sig_C/C)^2+(sig_bracket/Bracket)^2);
sig_top = sqrt(sig_gbc^2+sig_F^2);
sig_angfreq2 = angfreq2*sqrt((sig_top/(B+F))^2+(sig_A/A)^2);
sig_angfreq = 0.5*sig_angfreq2/angfreq2;
propagated_error = sig_angfreq/(2*pi);
end