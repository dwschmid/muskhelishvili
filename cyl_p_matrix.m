%==============================================================================================================
% CYL_P_MATRIX.M
%
% Pressure around circular inclusion
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

%FAR FIELD FLOW - VISCOSITIES - GEOMETRY
gr	    = 1;
er      = 0;
mm      = 1;
mc      = 1e6;
rc	    = 1;

%PRESSURE CALCULATION IN THE Z-PLANE
[X,Y]	= meshgrid(-2:.01:2);
Z	    = X+i*Y;
P	    = -2.*mm.*(mc-mm)./(mc+mm).*real(rc^2./Z.^2.*(i*gr+2*er));

%PRESSURE IS ONLY FOR THE OUTSIDE OF THE CLAST
P(abs(Z)<rc) = NaN;

%PLOTTING
pcolor(X,Y,P)
axis image;
shading interp;
hold on;
contour(X,Y,P, [-1.5,-1,-.5,0,.5,1,1.5], 'k');
colorbar('horiz');
