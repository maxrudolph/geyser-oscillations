% Application of geyser model to old faithful
clear;
% close all;
addpath XSteam_Matlab_v2.6;
addpath Ideal_Gas;
addpath Steam
%closed = true;

% OFG parameters
par.sb = 40; % 300 m^2 
par.sc = 0.5;   % 1 m^2
par.g = 9.81;% gravitational acceleration in m*s^-2
par.rho = 1000; % water density in kg*m^-3
par.gamma=7/5; % adiabatic exponent for diatomic ideal gas - unitless
par.alpha = 5/2; %diatomic ideal gas constant - unitless
par.Pa0 = 1e5; %atmospheric pressure at equilibrium in Pa
par.H = 8;
par.sl=par.sc;
par.L=5;

% EL Jefe Parameters:
par.sb = 4;
par.sc = 0.07;
par.sl = par.sc;
par.H = 1;
par.L = 1;
% par.Pa0 = 6.1e4;


nx=50;
% delta = 1e-6;%bar
%xx = linspace(4,6.999,nx);
xx = logspace(-3,log10(par.H),nx);% actually H-xbar
ny = 51;
yy = linspace(0,30,ny);

freq = zeros(ny,nx);
for a=1:nx
    for j=1:ny
        if( yy(j) >= par.H- xx(a) )
            par.xbar = par.H-xx(a);
            par.ybar = yy(j);
            par.delxy = par.ybar-par.xbar;
            par.Vol_0 = par.sb*(par.H-par.xbar); %initial (equilibrium) volume in m^3
            P_0 = par.rho*par.g*(par.delxy) + par.Pa0; %initial pressure in (Pascal)
            xv = 1 - 1e-2;
            % xv=0.5;
            sv=XSteam('sV_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
            sl=XSteam('sL_p',P_0/1e5)*1e3; %specific entropy J*kg^-1*degC^-1
            s = xv*sv + (1-xv)*sl; % gives specific entropy of water+vapor mixture in J/kg/C
            V0 = XSteam('V_ps',P_0/1e5,s/1e3); % volume per unit mass
            par.m = par.Vol_0 / V0; %vapor mass in kg
            
            P=linspace((P_0-9e4)/1e5,(P_0+9e4)/1e5,51); %Range of pressure for lookup table (bar)
            vV=zeros(1,length(P));  % vapor volume (m^3)
            du_dp = zeros(1,length(P)); %dU/dP in (kJ/K/Pa)
            dv_dp = zeros(1,length(P)); %dV/dP in (m^3/Pa)
            
            for i=1:length(P)
                delta = sqrt(eps(P(i)));
                vV(i) = XSteam('v_ps',P(i),s/1e3)*par.m; %volume in m^3
                xx1 = XSteam('x_ps',P(i),s/1e3);
                assert(xx1>0);
                du_dp(i) = (XSteam('u_ps',P(i)+delta,s/1e3)-XSteam('u_ps',P(i)-delta,s/1e3))*1e3*par.m/(2*delta*1e5);% enthalpy in mks units, pressure in Pa
                dv_dp(i) = (XSteam('v_ps',P(i)+delta,s/1e3)-XSteam('v_ps',P(i)-delta,s/1e3))*par.m/(2*delta*1e5);% m^3/Pa
            end
            [~,i] = sort(vV);
            %vV = vV(~isnan(vV)); du_dp = du_dp(~isnan(du_dp)); dv_dp=dv_dp(~isnan(dv_dp));
            par.Fdudp = griddedInterpolant(vV(i),du_dp(i),'linear','none');
            par.FdPdv = griddedInterpolant(vV(i),1./dv_dp(i),'linear','none');
            [AA,~] =frequency_calculator(par);
            freq(j,a) = AA;
        else
            freq(j,a) = NaN;
        end
    end
end

%% compute the u-tube freq
Futube = zeros(ny,nx);
for a=1:nx
    for j=1:ny
        par.xbar = par.H-xx(a);
        par.ybar = yy(j);
        par.delxy = par.ybar-par.xbar;
        LL = par.xbar + par.ybar + par.L;
        Futube(j,a) = 1/(2*pi)*sqrt(2*9.81/LL);
    end
end
figure, 
contourf(xx,yy,Futube,1024,'Color','none');
hold on
contour_levels = [0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 1.0 1.2 1.4 1.6 1.8 2 4 6 8 10]
contour(xx,yy,(freq),contour_levels,'Color','k','ShowText','on');

set(gca,'ColorScale','log');
hold on
hcb=colorbar;

%% 
figure;
contourf(xx,yy,(freq),1024,'Color','none');
set(gca,'ColorScale','log');
hold on
hcb=colorbar;
% 1 Hz contour
contour_levels=logspace(-1,2,25);
contour_levels = [0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 1.0 1.2 1.4 1.6 1.8 2 4 6 8 10]
contour(xx,yy,(freq),contour_levels,'Color','k','ShowText','on');
fn = get(gca,'FontName')
set(gca,'FontSize',16);
xlabel('Trapped Vapor Thickness (${H}-\bar{x}$) (m)','Interpreter','latex');
ylabel('Conduit Fill Level ($\bar{y}$) (m)','Interpreter','latex','FontName','Helvetica');
hcb.Label.String = 'Frequency (Hz)';
hcb.Label.FontSize=16;
set(gcf,'Color','w');
set(gca,'XScale','log');
set(gca,'XTick',10.^[-3 -2 -1 0 1])
%saveas(gcf,'ofg_application.svg');
% export_fig('ofg_application.eps');
exportgraphics(gcf,'Figure8-ofg_application_steam.pdf','ContentType','vector')


%% Now, solve the ODE directly for ofg.
p.Sb=par.sb;
p.Sc=par.sc;
p.ybar = 20;
p.H = 7;
p.xbar = p.H-.01;
p.rho = 1000;
p.g = 10;
p.patm=1e5;
p.p0 = p.patm+p.rho*p.g*(p.ybar-p.xbar);
p.Sl = 1;
p.Ll = 0;
p.gamma=7/5;
p.alpha=5/2;
p.Fa=-1.0;

time=     [0 20]; % starting and ending times for time integration

% define initial conditions
xx0=zeros(2,1);
xx0(1)=p.Fa * p.Sc/p.Sb;% x
xx0(2)=0.0; % dx/dt

options=odeset('MaxStep',1e-3);
odefun = @(t,x) sohn_model3_dimensional(t,x,p); % use anonymous function to pass parameters as argument
% solve the problem
[t,y] = ode45( odefun, time, xx0, options);
h=figure;
h.Position(3:4) = [1000 500];
subplot(1,2,1);
plot(t,-y(:,1)*p.Sb/p.Sc);
xlabel('Time (s)');

ylabel('y (m)');
axis square

subplot(1,2,2);
plot(y(:,1)*p.Sb/p.Sc,y(:,2)*p.Sb/p.Sc);
xlabel('y (m)');
ylabel("dy/dt (m/s)");
hold on
axis square;
