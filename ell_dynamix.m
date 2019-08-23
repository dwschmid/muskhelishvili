%==============================================================================================================
% ELL_DYNAMIX.M
%
% 2D pressure, stress and maximum shear stress caused by an elliptical inclusion
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


%COMPLEX NUMBER DEFINITION
I       = sqrt(-1);
i       = sqrt(-1);

%VISCOSITY CONTRAST BETWEEN CLAST AND MATRIX
mc      = 1000;

%FAR FIELD FLOW
er      = 0;
gr      = 1;
alpha   = -30/180*pi;

%ASPECT RATIO, t CANNOT BE 1 or SMALLER, USE t=1.001 FOR CIRCULAR INCLUSION APPROXIMATION
t       = 2;
rc      = sqrt((t-1)*(t+1))/(t-1);

%SOLUTION CONSTANTS
BC 	    = (2.*er-I.*gr).*exp(+2.*I.*alpha);
B1      = rc.^4.*mc+rc.^4-1+mc;
B2      = rc.^4.*mc+rc.^4-mc+1;
B3      = rc.^4.*mc-mc-rc.^4+1;
B4      = -rc.^4.*mc-mc-rc.^4+1;
B5      = rc.^8.*mc-mc-rc.^8+1;

%RESOLUTION
rs      = 100;
ts      = 200;

%CLAST GRID IS FROM 1..rc
[rho, theta]    = meshgrid(1:(rc-1)/rs:rc, 0:2*pi/ts:2*pi);
zeta_clast      = rho.*exp(i*theta);
p_clast         = real(-I.*mc.*B4./B1.*gr+2.*rc.^2.*(mc-1).*(I.*mc.*imag(BC)./B1-real(BC)./B2));
tau_clast       = -2.*mc.*rc.^4.*(I.*imag(BC)./B1+real(BC)./B2);
tau_clast       = sqrt((real(tau_clast)).^2 + (imag(tau_clast)).^2);
%Correct size of arrays
p_clast         = ones(size(rho))*p_clast;
tau_clast       = ones(size(rho))*tau_clast;

%MATRIX
[rho, theta]    = meshgrid(rc:2*rc/rs:3*rc, 0:2*pi/ts:2*pi);
zeta_mat        = rho.*exp(i*theta);
p_mat           = -2.*rc.^2.*real(B3.*(-I.*imag(BC).*B2+real(BC).*B1)./(zeta_mat.^2-1)./B1./B2);
str_mat         = conj(zeta_mat+1./zeta_mat).*((-I.*gr./zeta_mat.^3+2.*B3.*rc.^2.*(I.*imag(BC)./B1-real(BC)./B2)./zeta_mat.^3)./(1-1./(zeta_mat.^2)).^2-2.*(-1./2.*I.*gr.*(1-1./(zeta_mat.^2))-B3.*rc.^2.*(I.*imag(BC)./B1-real(BC)./B2)./zeta_mat.^2)./(1-1./(zeta_mat.^2)).^3./zeta_mat.^3)+(-(real(BC)+I.*imag(BC)).*(1-1./(zeta_mat.^2))-B5.*(I.*imag(BC)./B1-real(BC)./B2)./(zeta_mat.^3-zeta_mat).^2.*(3.*zeta_mat.^2-1))./(1-1./(zeta_mat.^2));
tau_mat         = sqrt((real(str_mat)).^2 + (imag(str_mat)).^2);

%TRANSLATE zeta -> z
z_clast         = zeta_clast+1./zeta_clast;
z_mat           = zeta_mat+1./zeta_mat;

%PLOT IN ZETA (IMAGE) DOMAIN
figure(1);
clf
subplot(211)
pcolor(real(zeta_clast), imag(zeta_clast), p_clast);
hold on;
pcolor(real(zeta_mat),   imag(zeta_mat),   p_mat);
shading interp;
axis image
axis off
title('Pressure')
colorbar
colormap(jet)

subplot(212)
pcolor(real(zeta_clast), imag(zeta_clast), tau_clast);
hold on;
pcolor(real(zeta_mat),   imag(zeta_mat),   tau_mat);
shading interp;
axis image
axis off
title('\tau')
colorbar
colormap(jet)

%PLOT IN Z (PHYSICAL) DOMAIN
figure(2);
clf
subplot(211)
pcolor(real(z_clast), imag(z_clast), p_clast);
hold on;
pcolor(real(z_mat),   imag(z_mat),   p_mat);
shading interp;
axis image
axis off
title('Pressure')
colorbar
colormap(jet)

subplot(212)
pcolor(real(z_clast), imag(z_clast), tau_clast);
hold on;
pcolor(real(z_mat),   imag(z_mat),   tau_mat);
shading interp;
axis image
axis off
title('\tau')
colorbar
colormap(jet)
