function y = nonlin(x,threshold,slope,exponent)
% nonlinearity (x: MP, y: FR)
if nargin < 4, exponent = 1; end
if nargin < 3, slope = 1; end
if nargin < 2, threshold = 0; end
y = slope*(x-threshold).^exponent;
y(x<threshold) = 0;
end