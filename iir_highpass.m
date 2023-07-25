clc;close all; clear all;
k=0:255;
w=(pi/255).*k   % to chosse 256 points between 0 and pi
alpha=0.767326
const=(1+alpha)/2
 
H=const.*((1-exp(-1*j*w))./(1-(alpha.*exp(-1*j*w))))
figure;
subplot(2,1,1)
plot(w/pi,abs(H)); grid on
title('Magnitude Spectrum |H(e^{j\omega})|')
xlabel('\omega /\pi');
ylabel('Amplitude');
subplot(2,1,2)
plot(w/pi,angle(H)); grid on
title('Phase Spectrum arg[H(e^{j\omega})]')
xlabel('\omega /\pi');
ylabel('Phase, radians')
 
b=const.*[1 1] % co-efficient vectors of numerator and denominater i.e poles and zero
a=[1 alpha]

figure;
[gdelay,w1]=grpdelay(b,a) % callculating group deay for filter
plot(w1,gdelay)
grid on
title('group delay')
xlabel('w')
ylabel('group delay')
 
figure;
zplane(b,a) % ploting poles and zeros on pzmap
grid on
title('pzmap')
