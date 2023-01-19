function propagated_error = propagation(par,model)
if model==1 %steam error propagation
    hx = par.H - par.xbar;
    sig_area = 5e-6; %uncertainty in measures of area, 1 cm^2
    sig_len = 1e-3; %uncertainty in measures of length, 1 mm
    sig_SbSc = par.sb/par.sc *sqrt((sig_area/par.sb)^2 +(sig_area/par.sc)^2);
    sig_gbc = par.g*sig_SbSc;
    sig_hx = sqrt(2*sig_len^2);
    sig_invhx = sig_hx/(hx)^2;
    sig_invhx2 = 2/hx^3;
    sig_bracket = sqrt(2*sig_invhx^2+sig_invhx2^2);
    sig_rs = par.rho;
    sig_invrs = sig_rs*(par.rho*par.sb)^-2;
    [bracket,top,bottom] = bracket_val(par,model);
    sig_bigterm = (par.rho*par.sb)^-1 *bracket*sqrt((sig_invrs*par.rho*par.sb)^2 + (sig_bracket/bracket)^2);
    sig_top = sqrt(sig_gbc^2+sig_bigterm^2);
    sig_bcy = par.ybar*(par.sb/par.sc)*sqrt( (sig_SbSc/(par.sb/par.sc))^2 + par.ybar^-2 );
    sig_bottom = sqrt(1+sig_bcy^2);
    sig_omega2 = (frequency_calculator(par))^2*sqrt((sig_top/top)^2+(sig_bottom/bottom)^2);
    propagated_error = 0.5*sig_omega2*(2*pi)^-1;
end
if model==0
    hx = par.H - par.xbar;
    sig_area = 5e-6; %uncertainty in measures of area, 1 cm^2
    sig_len = 1e-3; %uncertainty in measures of length, 1 mm
    sig_SbSc = par.sb/par.sc *sqrt((sig_area/par.sb)^2 +(sig_area/par.sc)^2);
    sig_gbc = par.g*sig_SbSc;
    sig_hx = sqrt(2*sig_len^2);
    sig_invhx = sig_hx/(hx)^2;
    sig_bracket = sig_invhx^2;
    sig_rs = par.rho;
    sig_invrs = sig_rs*(par.rho*par.sb)^-2;
    [~,~,bottom] = bracket_val(par,model);
    sbsc= par.sb ./ par.sc;
    top = par.g.*(1+sbsc)+1/(par.rho*par.sb)*(par.gamma/hx);
    bracket = par.Pb0*par.gamma;
    sig_bigterm = (par.rho*par.sb)^-1 *bracket*sqrt((sig_invrs*par.rho*par.sb)^2 + (sig_bracket/bracket)^2);
    sig_top = sqrt(sig_gbc^2+sig_bigterm^2);
    sig_bcy = par.ybar*(par.sb/par.sc)*sqrt( (sig_SbSc/(par.sb/par.sc))^2 + par.ybar^-2 );
    sig_bottom = sqrt(1+sig_bcy^2);
    sig_omega2 = (coldfreq(par))^2*sqrt((sig_top/top)^2+(sig_bottom/bottom)^2);
    propagated_error = 0.35*sig_omega2*(2*pi)^-1;
end
end