%==============================================================================================================
% CYL_P_INTERF.M
%
% Pressure in the matrix at the circular clast-matrix interface. 
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

%SETUP ANGLES
theta   = 0:2*pi/359:2*pi;

%VISCOSITIES
mm      = 1;
mcs     = [1,2,10,1e6];

%FAR FIELD FLOW
er      = 0;
gr      = 1;

%CALCULATE AND PLOT
Styles  = {':k', '-.k', '--k', '-k'};
figure(1);
clf
counter=1;
for mc=mcs;
    Pressure    = 2.*mm.*(-2.*mc.*cos(2.*theta).*er-mc.*gr.*sin(2.*theta)+2.*mm.*cos(2.*theta).*er+mm.*gr.*sin(2.*theta))./(mc+mm);    
    plot(theta/pi*180, Pressure, Styles{counter});
    hold on;
    counter     = counter+1;
end

axis([0 360 -2 2]);
set(gca, 'XTick', [0:45:360]);
xlabel('\theta');
ylabel('p/(\mu_m\gamma)', 'Rotation', 0);
title('Pressure Around Cylindrical Inclusion')
legend({'\mu_c/\mu_m=1','\mu_c/\mu_m=2','\mu_c/\mu_m=10','\mu_c/\mu_m=\infty'}, 'location', 'NorthEastOutside');
	 