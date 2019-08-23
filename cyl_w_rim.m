%==============================================================================================================
% CYL_W_RIM
%
% Complete two-dimensional field of pressure, maximum shear stress and stream function
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

%DEFINE i
i   = sqrt(-1);
I   = sqrt(-1);

%INPUT PARAMETERS
er  = -0;   %Negative values indicate horizontal compression
gr  = 1;    %Positive value indicate top to the left shear
rl  = 1.2;
ml  = 1e+3;
mc  = 1e+3;

%DESIRED BOX SIZE
bs  = 3; %Times rl

%K's
K1 = ml*mc;
K2 = (ml-mc)*(ml-1);
K3 = (mc+ml)*(ml+1);
K4 = ml/(ml-1);
K5 = ml/(ml+1);
K6 = (ml-mc)*(ml+1);
K7 = (mc+ml)*(ml-1);
K8 = (ml-mc)*(mc+ml)*(1-ml+ml^2);
K9 = (ml-mc)*(mc+ml)*(3-8*ml+3*ml^2);

%Q's
Q0 = K2^2+K2*K3*(-4*rl^2+6*rl^4-4*rl^6)+K3^2*rl^8;
Q1 = 4*K1*K2*(rl^2-rl^4)/Q0;
Q2 = (-16*K1*K2*rl^2+12*K1*K2*rl^4+4*K1*K3*rl^8)/Q0;
Q3 = K2*K4*(-2*K2*rl^2+2*K3*rl^8)/Q0;
Q4 = K2*(ml^2+K1)*(2*rl^2-2*rl^4)/Q0;
Q5 = K2*K4*(-2*K2*rl^4+2*K3*rl^8)/Q0;
Q6 = K3*(K2*K5*(-8*rl^2+6*rl^4)+2*(ml^2+K1)*rl^8)/Q0;
Q7 = (-K2*K6*rl^2+4*K2*K7*rl^4-6*K2*K7*rl^6+4*K8*rl^8-K3*K7*rl^10)/Q0;
Q8 = (-K2*K6*rl^4+4*K2*K7*rl^6-2*K9*rl^8+4*K2*K7*rl^10-K3*K7*rl^12)/Q0;

%RESOLUTION
nr          = 100;
nt          = 200;
Theta       = 0:2*pi/nt:2*pi;

%CLAST
[R, THETA]  = meshgrid(0:1/nr:1, Theta);
z           = R.*exp(i*THETA);
x           = real(z);
y           = imag(z);
Z_CLAST     = z; %Save z's for later contour plots
%Pressure
figure(1);
clf
PRES_CLAST  = -6*Q1*real(z.^2*(-2*er+i*gr));
pcolor(real(z), imag(z), PRES_CLAST);
hold on;
%Streamfun
figure(2);
clf
STREAM_FUN_CLAST = er.*(-2./mc.*Q1.*y.^3.*x-1./mc.*Q2.*y.*x-2./mc.*Q1.*x.^3.*y)+(-1./4.*y.^2-1./4.*1./mc.*Q2.*y.^2-1./2.*1./mc.*y.^4.*Q1).*gr;
pcolor(real(z), imag(z), STREAM_FUN_CLAST);
hold on;
%TAU
figure(3);
clf
TAU_CLAST = sqrt((-2.*Q2.*er+6.*Q1.*real(conj(z).*z.*(I.*gr-2.*er))).^2+(Q2.*gr+6.*Q1.*imag(conj(z).*z.*(I.*gr-2.*er))).^2);
pcolor(real(z), imag(z), TAU_CLAST);
hold on;


%LUBR
[R, THETA]  = meshgrid(1:(rl-1)/nr:rl, Theta);
z           = R.*exp(i*THETA);
x           = real(z);
y           = imag(z);
Z_LUBR      = z;
%Pressure
figure(1)
PRES_LUBR   = 2*real((i*gr+2*er)*Q3./z.^2-3*(-2*er+i*gr)*Q4*z.^2);
pcolor(real(z), imag(z), PRES_LUBR);
%Streamfun
figure(2)
STREAM_FUN_LUBR = -1./2./ml.*er.*(4.*x.*y.^3.*Q4+2.*y.*Q3./(x.^2+y.^2).*x-2.*Q5.*y.*x./(x.^2+y.^2).^2-x.*(-2.*Q3./(x.^2+y.^2).*y-6.*Q4.*(x.^2.*y-1./3.*y.^3))-2.*Q4.*(x.^3.*y-y.^3.*x)+2.*Q6.*y.*x)-1./2.*1./ml.*gr.*(1./2.*y.^2.*ml-x.*(Q3./(x.^2+y.^2).*x-3.*Q4.*x.*y.^2)-3./2.*Q4.*x.^2.*y.^2+3./4.*Q4.*y.^4-x.^2.*Q3./(x.^2+y.^2)-Q4.*(3./2.*x.^2.*y.^2-1./4.*y.^4)-Q5.*(-x.^2./(x.^2+y.^2).^2+1./2./(x.^2+y.^2))+1./2.*Q6.*y.^2);
pcolor(real(z), imag(z), STREAM_FUN_LUBR);
%TAU
figure(3);
TAU_LUBR    = sqrt(real(-conj(z).*(2.*(I.*gr+2.*er).*Q3./z.^3+6.*(I.*gr-2.*er).*Q4.*z)+3.*(I.*gr+2.*er).*Q5./z.^4-(I.*gr-2.*er).*Q6).^2+imag(-conj(z).*(2.*(I.*gr+2.*er).*Q3./z.^3+6.*(I.*gr-2.*er).*Q4.*z)+3.*(I.*gr+2.*er).*Q5./z.^4-(I.*gr-2.*er).*Q6).^2);
pcolor(real(z), imag(z), TAU_LUBR);

%MAT
[R, THETA]  = meshgrid(rl:(sqrt(2*bs^2*rl^2)-rl)/nr:sqrt(2*bs^2*rl^2), Theta);
z           = R.*exp(i*THETA);
x           = real(z);
y           = imag(z);
Z_MAT       = z;
%Pressure
figure(1)
PRES_MAT   = 2*real((i*gr+2*er)*Q7./z.^2);
pcolor(real(z), imag(z), PRES_MAT);
%Streamfun
figure(2)
STREAM_FUN_MAT = er.*(-2.*Q7.*y.*x./(x.^2+y.^2)+Q8.*y.*x./(x.^2+y.^2).^2-x.*y)+(-1./2.*y.^2+Q7.*x.^2./(x.^2+y.^2)+1./2.*Q8.*(-x.^2./(x.^2+y.^2).^2+1./2./(x.^2+y.^2))).*gr;
pcolor(real(z), imag(z), STREAM_FUN_MAT);
%TAU
figure(3);
TAU_MAT     = sqrt((-2.*er+real(2.*conj(z).*(I.*gr+2.*er).*Q7./z.^3-3.*(I.*gr+2.*er).*Q8./z.^4)).^2+(gr+imag(2.*conj(z).*(I.*gr+2.*er).*Q7./z.^3-3.*(I.*gr+2.*er).*Q8./z.^4)).^2);
pcolor(real(z), imag(z), TAU_MAT);


%FINALIZE PLOTS========================================
%Complex coordinates of clast and layer
Clast   =    exp(i*Theta);
Layer   = rl*exp(i*Theta);

figure(1)
axis off
axis equal
axis([-bs*rl bs*rl -bs*rl bs*rl]);
shading interp;
plot(real(Clast), imag(Clast), '--k');
plot(real(Layer), imag(Layer), '--k');
title(['P \epsilon:', num2str(er),' \gamma:', num2str(gr),' \mu_c: ',num2str(mc),' \mu_l: ',num2str(ml), ' r_l: ', num2str(rl)]);
cb_h        = colorbar('horiz.');
pos_cb      = get(cb_h, 'pos');
set(cb_h, 'pos', [.26666, pos_cb(2), .504, pos_cb(4)])
set(get(cb_h,'Title'),'String','P');

figure(2)
axis equal
axis([-bs*rl bs*rl -bs*rl bs*rl]);
shading interp;
plot(real(Clast), imag(Clast), '--k');
plot(real(Layer), imag(Layer), '--k');
colorbar('horiz.');
%Build up complete 2D STREAM_FUN
%Z           = [Z_CLAST;          Z_LUBR;          Z_MAT         ];
%STREAM_FUN  = [STREAM_FUN_CLAST; STREAM_FUN_LUBR; STREAM_FUN_MAT];
Z           = [ Z_MAT          ];
STREAM_FUN  = [ STREAM_FUN_MAT ];
contour(real(Z), imag(Z), STREAM_FUN, 10, 'k')
axis off;
axis equal
axis([-bs*rl bs*rl -bs*rl bs*rl]);
shading interp;
plot(real(Clast), imag(Clast), '--k');
plot(real(Layer), imag(Layer), '--k');
title(['\Theta \epsilon:', num2str(er),' \gamma:', num2str(gr),' \mu_c: ',num2str(mc),' \mu_l: ',num2str(ml), ' r_l: ', num2str(rl)]);
cb_h        = colorbar('horiz.');
pos_cb      = get(cb_h, 'pos');
set(cb_h, 'pos', [.26666, pos_cb(2), .504, pos_cb(4)])
set(get(cb_h,'Title'),'String','\Theta');

figure(3)
axis off
axis equal
axis([-bs*rl bs*rl -bs*rl bs*rl]);
shading interp;
plot(real(Clast), imag(Clast), '--k');
plot(real(Layer), imag(Layer), '--k');
title(['\tau \epsilon:', num2str(er),' \gamma:', num2str(gr),' \mu_c: ',num2str(mc),' \mu_l: ',num2str(ml), ' r_l: ', num2str(rl)]);
cb_h        = colorbar('horiz.');
pos_cb      = get(cb_h, 'pos');
set(cb_h, 'pos', [.26666, pos_cb(2), .504, pos_cb(4)])
set(get(cb_h,'Title'),'String','\tau');
