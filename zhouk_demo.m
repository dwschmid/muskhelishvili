%==============================================================================================================
% ZHOUK_DEMO
%
% Demonstration of the Joukowski transform
%
% 2002, Dani Schmid
%
% DISCLAIMER OF WARRANTY: 
% Since the Software is provided free of charge, the Software is provided on an AS IS basis,
% without warranty of any kind, including without limitation the warranties of merchantability,
% fitness for a particular purpose and non-infringement. The entire risk as to the quality and performance 
% of the Software is borne by you. Should the Software prove defective, 
% you assume the entire cost of any service and repair. 
%
% LIMITATION OF LIABILITY: 
% UNDER NO CIRCUMSTANCES AND UNDER NO LEGAL THEORY, TORT, CONTRACT, OR OTHERWISE, 
% SHALL THE AUTHORS BE LIABLE TO YOU OR ANY OTHER PERSON FOR ANY INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES OF ANY CHARACTER INCLUDING, WITHOUT LIMITATION, DAMAGES FOR LOSS OF GOODWILL, 
% WORK STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, OR ANY AND ALL OTHER COMMERCIAL DAMAGES OR LOSSES
%==============================================================================================================

%RESOLUTION
nt          = 200;
Theta       = 0:2*pi/nt:2*pi;

%SETUP THREE DIFFERENT CIRCLES
SLIT        =    exp(i*Theta);
ELLE        =  2*exp(i*Theta);
JOUK        = 2*exp(i*Theta)-0.9696+i*0.3473;

%PLOT IN ZETA
figure(1)
clf
subplot(121)
hold on;
plot(real(SLIT), imag(SLIT), '-k');
plot(real(ELLE), imag(ELLE), '--k');
plot(real(JOUK), imag(JOUK), ':k');
axis equal
grid on
title('\zeta-Plane');

%TRANSFORM 
SLIT    = SLIT + 1./SLIT;
ELLE    = ELLE + 1./ELLE;
JOUK    = JOUK + 1./JOUK;

%PLOT IN Z
subplot(122)
hold on;
plot(real(SLIT), imag(SLIT), '-k');
plot(real(ELLE), imag(ELLE), '--k');
plot(real(JOUK), imag(JOUK), ':k');
axis equal;
grid on
title('z-Plane');

