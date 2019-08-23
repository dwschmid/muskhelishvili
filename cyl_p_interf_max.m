%==============================================================================================================
% CYL_P_INTERF_MAX.M
%
% Maximum pressure in the matrix at the circular clast-matrix interface
% as a function of the viscosity contrast between clast and matrix. 
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

%LOGARITHMIC RANGE OF VISCOSITIES
mm          = 1;
mc          = logspace(-3,3);

%BOUNDARY CONDITION
gr          = 1;

%MAXIMUM PRESSURE EXPRESSION
Pressure    = 2.*mm.*(mc-mm)./(mc+mm).*gr;

%PLOT PRESSURE vs. THETA
figure(1);
clf
plot(mc, Pressure, '-k');
set(gca, 'XScale', 'log');
grid on;
set(gca, 'XMinorgrid', 'off');

xlabel('\mu_c/\mu_m');
ylabel('p/(\mu_m\gamma)', 'Rotation', 0);
title('Max. Pressure as f(\mu_c/\mu_m)')
