function [maxx,maxy,xfit,yfit] = maxfit(X,Y,dk,N)
% gives back the maximum of a curve using a parabolic fit
% data: X, Y
% dk: number data points used (for the fit) on one side
% xfit, yfit: parabola curve (for illustrative purposes)
if (nargin < 3)
    dk = 3;
end
dk1 = dk;
dk2 = dk;   
[~,maxk] = max(Y);
params = polyfit(X(maxk-dk1:maxk+dk2),Y(maxk-dk1:maxk+dk2),2);
maxx = -params(2)/params(1)/2;
maxx0 = X(maxk);
dx = X(2)-X(1);
if abs(abs(maxx-maxx0)-dx/2) < dx/6
    dk1 = dk-(sign(maxx-maxx0)+1)/2;
    dk2 = dk-1+(sign(maxx-maxx0)+1)/2;
    params = polyfit(X(maxk-dk1:maxk+dk2),Y(maxk-dk1:maxk+dk2),2);
    maxx = -params(2)/params(1)/2;
end
maxy = params(3)-params(2)^2/params(1)/4;
if (nargin < 4)
    N = 100;
end
xfit = linspace(X(maxk-dk1),X(maxk+dk2),N);
yfit = params(1)*xfit.^2+params(2)*xfit+params(3);
end