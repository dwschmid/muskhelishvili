%==============================================================================================================
% ELL_ROT_RATE
%
% Analytical formula for the rotation rate of an ellipse in combined, inclined  simple & pure shear
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
clf;

%FAR FIELD FLOW
er      = 0;
gr      = 1;
alpha   = -pi/2:(pi/2)/100:pi/2;

%ELLIPSE ASPECT RATIO
t       = 6;

%VISCOSITIES
mc      = [1e6,      1,  1/10, 1/100];
Styles  = {':k', '-.k', '--k', '-k'};

for m=1:length(mc)
    %ROTATAION RATE
	rot_rate = (-1./2.*(t.^2-mc(m).*t.^2+mc(m)-1)./ ...
                   (mc(m).*t.^2+mc(m)+2.*t).*cos(2.*alpha)-1./2).*gr ...
                   -1./2.*(2.*mc(m).*t.^2-2.*t.^2-2.*mc(m)+2)./ ...
                   (mc(m).*t.^2+mc(m)+2.*t).*sin(2.*alpha).*er;
    
    %PLOT
    plot(rot_rate, alpha/pi*180, Styles{m});
    hold on;
end
grid on;
axis tight;
xlabel('Rotation Rate');
ylabel('\alpha');
legend({'\mu_c/\mu_m=\infty','\mu_c/\mu_m=1','\mu_c/\mu_m=1/10','\mu_c/\mu_m=1/100'}, 'Location', 'NorthEastOutside')
