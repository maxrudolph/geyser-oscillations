function  [Pxx,Exx,pX,f]= mt_cspek_phs(x,Fs,unc);
%
%   Function "mt_cspekk" computes the power spectral densities of the 
% vectors "x" and "y", the cross-spectrum, transfer-function, and coherence
% using the multitaper method. 
% Compute "e" and "v" using the MATLAB function "dpss".
% 	The normalization is such that the PSD is defined from 0 to the 
% Nyquist frequency, and its integral over
% this range is the function mean squared amplitude."
% i.e. sum(Pxx*df) = var(x), where df = 1/(nfft*dt);
%
%	Pxx  = X-vector power spectral density
%	Pyy  = Y-vector power spectral density
%	Pxy  = Cross spectral density
%	Txy  = Complex transfer function from X to Y = Pxy./Pxx
%	Cxy  = Coherence function between X nd Y = (abs(Pxy).^2)./(Pxx.*Pyy)
%   pX  = phase of Sx(f), radians
%   pY  = phase of Sy(f), radians
%   pXY  = phase of Sxy(f), radians
%
%
%                                                         j.a. collins
%************************************************************************
%
% modifications by r.a. sohn, 6/06
% to include uncertainty estimates and phase outputs
% input param now includes unc, flag for reporting uncertainty
%
% 
%  DPSS   Discrete prolate spheroidal sequences (Slepian sequences).
%     [E,V] = DPSS(N,NW) are the first 2*NW discrete prolate spheroidal sequences
%     (DPSSs, or Slepian sequences) of length N (in the columns of E) and 
%     their corresponding concentrations (in vector V) in the frequency band 
%     |w|<=(2*pi*W) (where  W = NW/N is the half-bandwidth and w is in 
%     radians/sample).  E(:,1) is the length N signal most concentrated in the 
%     frequency band |w|<=(2*pi*W) radians/sample, E(:,2) is the signal 
%     orthogonal to E(:,1) which is most concentrated in this band, E(:,3) is the
%     signal orthogonal to both E(:,1) and E(:,2) which is most concentrated in 
%     this band, etc.  
%  
%     For multi-taper spectral analysis, typical choices for NW are 2, 5/2, 3, 
%     7/2, or 4.
%  
%     [E,V] = DPSS(N,NW,K) are the K most band-limited discrete prolate spheroidal
%     sequences.  [E,V] = DPSS(N,NW,[K1 K2]) returns the K1-th through the 
%     K2-th sequences.
%  
%     [E,V] = DPSS(N,NW,'spline') uses spline interpolation to compute the DPSSs 
%     from existing DPSSs in the DPSS database with length closest to N.
%     [E,V] = DPSS(N,NW,'spline',Ni) interpolates from existing length Ni DPSSs.
%     DPSS(N,NW,'linear') and DPSS(N,NW,'linear',Ni) use linear interpolation, 
%     which is much faster but less accurate than spline interpolation.  
%     'linear' requires Ni > N. [E,V] = DPSS(N,NW,'calc') uses the direct 
%     algorithm (default).
%  
%     Use a trailing 'trace' argument to find out which method DPSS uses, e.g.,
%     DPSS(...,'trace').
%
%************************************************************************
%
% modifications by t. barreyre, 03/2015
% to include jackknifed error estimates @ all frequencies for spectra, coherence, 
% transfer function and phase as well as correction for phase ambiguity (unwrapping)
% output param now includes: 
%    Er_p: jackknifed error estimates for phase 
%    Er_c: jackknifed error estimates for coherence 
%    Er_Tf: jackknifed error estimates for transfer function
%

if (nargin < 3) 
    error ('Not enough input arguments!');
end

x = x(:);		       % Make sure x is a column vector
x = detrend(x,0);      % remove mean
% y = y(:);		       % Make sure x is a column vector
% y = detrend(y,0);      % remove mean

N = length(x);
[e,v] = dpss(N,3); % from Rob's code in Rudolph and Sohn 2018

nfft = max(256,2^nextpow2(N));
% compute PSD for x
k = length(v);


%%%%% Compute the windowed dfts and spectral estimates for x %%%%%%%%%%%%%
if N<=nfft
    Sk_xx = abs(fft(e(:,1:k).*x(:,ones(1,k)),nfft)).^2;
    Sk_x  =     fft(e(:,1:k).*x(:,ones(1,k)),nfft);
else
    % use CZT to compute DFT on nfft evenly spaced samples around the
    % unit circle:
    Sk_xx = abs(czt(e(:,1:k).*x(:,ones(1,k)),nfft)).^2;
    Sk_x  =     czt(e(:,1:k).*x(:,ones(1,k)),nfft);
end

% Select the proper points from fft:
if isreal(x)
    if rem(nfft,2)==0, M=nfft/2+1; else M=(nfft+1)/2; end
else
    M = nfft;
end
Sk_xx = Sk_xx(1:M,:);
Sk_x  = Sk_x(1:M,:);

%%%%% Compute the windowed dfts and spectral estimates for y %%%%%%%%%%%%%
% Compute the windowed dfts and the
%   corresponding spectral estimates:
% if N<=nfft
%     Sk_yy = abs(fft(e(:,1:k).*y(:,ones(1,k)),nfft)).^2;
%     Sk_y  =     fft(e(:,1:k).*y(:,ones(1,k)),nfft);
% else
%     % use CZT to compute DFT on nfft evenly spaced samples around the
%     % unit circle:
%     Sk_yy = abs(czt(e(:,1:k).*y(:,ones(1,k)),nfft)).^2;
%     Sk_y  =     czt(e(:,1:k).*y(:,ones(1,k)),nfft);
% end
% 
% % Select the proper points from fft:
% if isreal(y)
%     if rem(nfft,2)==0, M=nfft/2+1; else M=(nfft+1)/2; end
% else
%     M = nfft;
% end
% Sk_yy = Sk_yy(1:M,:);
% Sk_y  = Sk_y(1:M,:);


% Find the adaptive weights via iteration.  The algorithm converges 
% so fast that results are usually 'indistinguishable' after about 
% three iterations.  This version uses the equations from P&W pp 368-370

%%%%%%%% Get adaptive weights for x %%%%%%%%%%%%%%%%%%%
sig2=x'*x/N;                    % Power
S=(Sk_xx(:,1)+Sk_xx(:,2))/2;    % Initial spectrum estimate
Stemp=zeros(M,1);
S1=zeros(M,1);

% Set tolerance for acceptance of spectral estimate:
tol=.0005*sig2/M;
i=0;
a=sig2*(1-v);


% Do the iteration:
while sum(abs(S-S1)/M)>tol
    i=i+1;
    % calculate weights
    b=(S*ones(1,k))./(S*v'+ones(M,1)*a');
    % calculate new spectral estimate
    wk=(b.^2).*(ones(M,1)*v');
    S1=sum(wk'.*Sk_xx')./ sum(wk');
    S1=S1';
    Stemp=S1; S1=S; S=Stemp;  % swap S and S1
end
S_xx = S;
b_x = b; wk_x = wk;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Get adaptive weights for y %%%%%%%%%%%%%%%%%%%
% sig2=y'*y/N;                    % Power
% S=(Sk_yy(:,1)+Sk_yy(:,2))/2;    % Initial spectrum estimate
% Stemp=zeros(M,1);
% S1=zeros(M,1);
% 
% % Set tolerance for acceptance of spectral estimate:
% tol=.0005*sig2/M;
% i=0;
% a=sig2*(1-v);
% 
% % Do the iteration:
% while sum(abs(S-S1)/M)>tol
%     i=i+1;
%     % calculate weights
%     b=(S*ones(1,k))./(S*v'+ones(M,1)*a');
%     % calculate new spectral estimate
%     wk=(b.^2).*(ones(M,1)*v');
%     S1=sum(wk'.*Sk_yy')./ sum(wk');
%     S1=S1';
%     Stemp=S1; S1=S; S=Stemp;  % swap S and S1
% end
% S_yy = S;
% b_y = b; wk_y = wk;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Cross Spectrum %%%%%%%%% 
%wk_xy = (b_x.*b_y).*(ones(M,1)*v');
%S_tmp = conj(Sk_x).*Sk_y; 
%S_xy = sum(wk_xy'.*S_tmp')./ sum(wk_xy');
%S_xy = S_xy';

%cross_num = zeros(length(Sk_x(:,1)),1);
%cross_denom = zeros(length(Sk_x(:,1)),1);
%for n = 1:k
%    cross_num   = v(n)*b_x(:,n).*b_y(:,n).*conj(Sk_x(:,n)).*Sk_y(:,n) + cross_num;
%    cross_denom = v(n)*b_x(:,n).*b_y(:,n) + cross_denom;
%end
%S_xy = cross_num./cross_denom; 

 

%%%%% Get Chave's defn of adaptive weights  %%%%% 
B = var(x)*(1-v);
for j = 1:k
    wt(:,j) = sqrt(v(j))*S_xx(:)./(v(j)*S_xx(:) + B(j));
end
for i = 1:M
    sumwt(i) = sum (wt(i,:).^2);
end
sumwt = sqrt(sumwt);
for j = 1:k
    wt(i,j) = wt(i,j)/sumwt(j);
end
wt_x = wt;

% B = var(y)*(1-v);
% for j = 1:k
%     wt(:,j) = sqrt(v(j))*S_yy(:)./(v(j)*S_yy(:) + B(j));
% end
% for i = 1:M
%     sumwt(i) = sum (wt(i,:).^2);
% end
% sumwt = sqrt(sumwt);
% for j = 1:k
%     wt(i,j) = wt(i,j)/sumwt(j);
% end
% wt_y = wt;

%for i = 1:M
%    for j = 1:k
%        eigcr(j)   = v(j)*wt_x(i,j)*wt_y(i,j)*conj(Sk_x(i,j))*Sk_y(i,j);
%        eigsp_x(j) = v(j)*(wt_x(i,j)*abs(Sk_x(i,j)))^2;
%        eigsp_y(j) = v(j)*(wt_y(i,j)*abs(Sk_y(i,j)))^2;
%    end
%    cross(i) = sum(eigcr)/k;
%    spec_x(i) = sum(eigsp_x)/k;
%    spec_y(i) = sum(eigsp_y)/k;
%    c2(i) = abs(cross(i))^2/(spec_x*spec_y);
%    phase(i) = atan2(imag(cross(i)),real(cross(i)));
%end
%S_xx = spec_x(:);
%S_yy = spec_y(:);
%S_xy = cross(:);


for j = 1:k
%     eigcr(:,j)   = v(j).*wt_x(:,j).*wt_y(:,j).*conj(Sk_x(:,j)).*Sk_y(:,j);
    eigsp_x(:,j) = v(j).*(wt_x(:,j).*abs(Sk_x(:,j))).^2;
%     eigsp_y(:,j) = v(j).*(wt_y(:,j).*abs(Sk_y(:,j))).^2;
    Sx(:,j) = v(j).*wt_x(:,j).*Sk_x(:,j);
%     Sy(:,j) = v(j).*wt_y(:,j).*Sk_y(:,j);
end
% cross = sum(eigcr,2)/k;
spec_x = sum(eigsp_x,2)/k;
% spec_y = sum(eigsp_y,2)/k;
S_x = sum(Sx,2)/k;
% S_y = sum(Sy,2)/k;
S_xx = spec_x(:);
% S_yy = spec_y(:);
% S_xy = cross(:);
pX=angle(S_x);
% pY=angle(S_y);
% pXY  = angle(cross);

dt = 1/Fs;

% For real signals return the single-sided spectrum with full power.
if isreal(x),    
   S_xx = [S_xx(1); 2*S_xx(2:end-1); S_xx(end)];
end
% if isreal(y),    
%    S_yy = [S_xx(1); 2*S_yy(2:end-1); S_yy(end)];
% end
% if (isreal(x) & isreal(y)),    
%    S_xy = [S_xx(1); 2*S_xy(2:end-1); S_xy(end)];
% end
Pxx = dt*S_xx;
% Pyy = dt*S_yy;
% Pxy = dt*S_xy;
% Cxy = (abs(Pxy).^2)./(Pxx.*Pyy);
% Txy = Pxy./Pxx;

% confidence interval estimation, power spectra DOF from P&W, p. 370
% 
alp=0.95;
alp1=(1-alp)/2;
Exx=zeros(M,2);
Eyy=zeros(M,2);
ETxy=zeros(M,2);

DOFxys=NaN(1,M);
if unc==1	% if uncertainty flag is on...
  for i=1:M
    DOFx(i)=round(2*(b_x(i,:).^2*v).^2/(b_x(i,:).^4*v.^2));	%P&W p.370
%     DOFy(i)=round(2*(b_y(i,:).^2*v).^2/(b_y(i,:).^4*v.^2));    
%     DOFxy=mean([DOFx(i) DOFy(i)]);
%     DOFxys(i)=DOFxy; %Let's output the DOFs at each point
    Exx(i,2)=Pxx(i)*chi2inv(1-alp1,DOFx(i))/DOFx(i);
    Exx(i,1)=Pxx(i)*chi2inv(alp1,DOFx(i))/DOFx(i);
%     Eyy(i,2)=Pyy(i)*chi2inv(1-alp1,DOFy(i))/DOFy(i);
%     Eyy(i,1)=Pyy(i)*chi2inv(alp1,DOFy(i))/DOFy(i);
%     varT=abs(Txy(i))^2/(Cxy(i)*DOFxy);	% Bendat and Piersol, p. 320
%     ETxy(i,1)=norminv(alp1,abs(Txy(i)),sqrt(varT));
%     ETxy(i,2)=norminv(1-alp1,abs(Txy(i)),sqrt(varT));
    
    % Zero Significant Level of coherency (zsl)
    % where alpha is desired confidence interval, 0< alpha ? 1.0
    alpha=0.9;
%     zsl(i)=1-(1-alpha)^(2/(DOFxy-2));
    % zsl=1-(1-alpha).^(2./(DOFxys-2)); if we calculate it outside the loop
    
    %lines uncommented, cwp 1/11
%     r=(4/(DOFxy-4))*finv(alp1,4,DOFxy-4)*(1-Cxy(i))*Pyy(i)/Pxx(i);	% chave notes, 
    %ETxy(i,1)=abs(Txy(i))-r;
    %ETxy(i,2)=abs(Txy(i))+r;
  end
end


df = Fs/nfft;
F_nyq = Fs/2;
f = [0:df:F_nyq];
f = f(:);


% remove zero-frequency component
Pxx(1) = [];
% Pxy(1) = [];
% Pyy(1) = [];
% Cxy(1) = [];
% Txy(1) = [];
% ETxy(1,:)=[];
Exx(1,:)=[];
% Eyy(1,:)=[];
pX(1)=[];
% pY(1)=[];
% pXY(1) = [];
f(1) = [];
% DOFxys(1) = [];
% zsl(1) = [];


% normalize
%df = f(2)-f(1);
%Pxx = Pxx/df;
% The variable check should be ~= 1
check = sum(Pxx*df)/var(x)

% normalize

df = f(2)-f(1);
scale_x = var(x)./(sum(Pxx*df));
Pxx = Pxx*scale_x;
% scale_y = var(y)./(sum(Pyy*df));
% Pyy = Pyy*scale_y;


% %% Error calculation - bootstraping on the independent fourier estimates (eigsp_x/_y & eigcr)
% %modifications by t. barreyre, 03/2015 
% 
% % Calculation of the phase estimate error: Er_p
% jackstatp = jackknife(@mean,eigcr')';
% pXYj  = angle(jackstatp);
% pXYj(1,:) = []; 
% 
% %filter: not estimating phase lag error if coherence Cxy is < zsl
% %pXYj(Cxy<zsl',:) = NaN;
% 
% %solving phase ambiguity (i.e., correction) to compute correct average of jackknifed phase vector
% pXYj_pi = pXYj; %pXYj_pi phase values are in [-pi pi]
% pXYj(pXYj<0) = pXYj(pXYj<0)+2*pi; %put all phases in [0 2pi[
% idw = find(max(pXYj,[],2)-min(pXYj,[],2) > pi); %index for unwarping
% for ib=1:length(idw)
%     pXYj(idw(ib),pXYj(idw(ib),:)>pi) = pXYj(idw(ib),pXYj(idw(ib),:)>pi)-2*pi;
%     pXYj_pi(idw(ib),pXYj_pi(idw(ib),:)<0) = pXYj_pi(idw(ib),pXYj_pi(idw(ib),:)<0)+2*pi;
% end
% idwf = find(std(pXYj,[],2) > std(pXYj_pi,[],2)); %flag on index for unwarping -> 60deg error in phase average due to ambiguity
% pXYj(idwf,:) = pXYj_pi(idwf,:);
% 
% mpXYj = mean(pXYj,2); MpXYj=repmat(mpXYj,1,k); %compute means of the jackknifed phase vector for each freq
% Er_p = sqrt(((k-1)./k).*sum(abs(pXYj-MpXYj).^2,2)); %compute jackknifed error estimates for phase
% 
% 
% % Calculation of the coherency and transfer function errors: Er_c & Er_Tf
% jackstatSxx = jackknife(@mean,eigsp_x')';
% jackstatSyy = jackknife(@mean,eigsp_y')';
% jackstatSxy = jackknife(@mean,eigcr')';
% 
% % For real signals return the single-sided spectrum with full power.
% if isreal(x),    
%    jackstatSxx = [jackstatSxx(1,:); 2*jackstatSxx(2:end-1,:); jackstatSxx(end,:)];
% end
% if isreal(y),    
%    jackstatSyy = [jackstatSxx(1,:); 2*jackstatSyy(2:end-1,:); jackstatSyy(end,:)];
% end
% if (isreal(x) && isreal(y)),    
%    jackstatSxy = [jackstatSxx(1,:); 2*jackstatSxy(2:end-1,:); jackstatSxy(end,:)];
% end
% 
% Pxxj = dt.*jackstatSxx;
% Pyyj = dt.*jackstatSyy;
% Pxyj = dt.*jackstatSxy;
% Cxyj = (abs(Pxyj).^2)./(Pxxj.*Pyyj);
% Txyj = Pxyj./Pxxj;
% Cxyj(1,:) = []; Txyj(1,:) = [];
% 
% mCxyj = mean(Cxyj,2); MCxyj=repmat(mCxyj,1,k); 
% mTxyj = mean(Txyj,2); MTxyj=repmat(mTxyj,1,k);
% Er_c = sqrt(((k-1)./k).*sum(abs(Cxyj-MCxyj).^2,2)); %compute jackknifed error estimates for coherence
% Er_Tf = sqrt(((k-1)./k).*sum(abs(Txyj-MTxyj).^2,2)); %compute jackknifed error estimates for transfer function

return;
