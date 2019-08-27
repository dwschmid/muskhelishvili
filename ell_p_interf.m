%==============================================================================================================
% ELL_P_INTERF.M
%
% Pressure in the matrix at the elliptical clast-matrix interface. 
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


%CLEAR FIGURE
figure(1);
clf

%COMPLEX NUMBER DEFINITION
I       = sqrt(-1);
i       = sqrt(-1);

%VISCOSITY CONTRAST BETWEEN CLAST AND MATRIX
mc      = 1e+3;

%FAR FIELD FLOW
er      = -1;
gr      = 0;
alpha   = 0/180*pi;

%ASPECT RATIO, t CANNOT BE 1 or SMALLER, USE t=1.001 FOR CIRCULAR INCLUSION APPROXIMATION
t       = [1.0001, 2, 10, 20];
Styles  = {':k', '-.k', '--k', '-k'};

%CIRCUMFERENCE
theta   = 0:2*pi/360:2*pi;

for m=1:length(t)
    %TRANSLATE ASPECT RATIO INTO RADIUS
    rc      = sqrt((t(m)-1)*(t(m)+1))/(t(m)-1);
    
    %PRESSURE ON RADIUS
    press   = (2.*mc.*rc.^4-2.*mc-2.*rc.^4+2).*rc.^2./(mc.*rc.^4+mc-1+rc.^4)./(1+rc.^4+mc.*rc.^4-mc).*((rc.^2.*(cos(theta).^2-sin(theta).^2)-1)./((rc.^2.*(cos(theta).^2-sin(theta).^2)-1).^2+4.*rc.^4.*sin(theta).^2.*cos(theta).^2).*(-2.*er.*cos(2.*alpha)-gr.*sin(2.*alpha)).*(mc.*rc.^4+mc-1+rc.^4)+2.*rc.^2.*sin(theta).*cos(theta)./((rc.^2.*(cos(theta).^2-sin(theta).^2)-1).^2+4.*rc.^4.*sin(theta).^2.*cos(theta).^2).*(-gr.*cos(2.*alpha)+2.*er.*sin(2.*alpha)).*(1+rc.^4+mc.*rc.^4-mc));
    
    %PLOTTING
    plot(theta/pi*180, press, Styles{m});
    hold on;
end
achsen  = axis;
axis([0 360 achsen(3:4)]);
set(gca, 'Xtick', [0:45:360])
grid on;

title('Pressure around elliptical inclusion');
xlabel('\theta');
ylabel('Pressure');
legend({'t=1','t=2','t=10','t=20'}, 'Location', 'NorthEastOutside');
