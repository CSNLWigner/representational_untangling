function Z = circular_gabor(FOV,gabor_params,grid_size)
% discretized circular Gabor filter on a Cartesian grid
% /distances and angles are measured in degrees/
% FOV (field of vision) = [left,right,top,bottom]
% left,right,top,bottom: angles measured from the line of sight
% simplified parameter passing is also accepted for FOV:
% FOV = [x,y] -> left = right = x and top = bottom = y
% FOV = [x] or simply x -> all angles are the same 
% FOV = [] -> default value (90) is used for all angles
% gabor_params = [x0,y0,lambda,theta,phi,sigma,rho]
% x0,y0: center of the filter measured from the line of sight
% lambda: wavelength of the cosine wave component
% theta: orientation of the parallel stripes
% phi: phase of the cosine at the center of the filter
% increasing (positive) phi -> upward shift
% sigma: standard deviation of the Gaussian envelope
% rho: used for the definition of the receptive field boundary  
% meaning of rho in the coordinate system aligned to the filter: 
% radius = rho*sigma at the boundary
% rho = Inf can be used for creating non truncated filters
% simplified parameter passing is also accepted for gabor_params:
% gabor_params = [x0,y0,lambda,theta,phi,sigma]
% -> rho = Inf is considered as default value
% gabor_params = [lambda,theta,phi,sigma]
% -> x0 = y0, rho = Inf
% grid_size = [Nx,Ny] or simply N = Nx = Ny
% Example: gabor([],[0,0,2,45,0,1,3.5],101)
if length(FOV) == 4
    left = FOV(1); right = FOV(2); top = FOV(3); bottom = FOV(4);
elseif length(FOV) == 2
    left = FOV(1); right = left; top = FOV(2); bottom = top;
elseif length(FOV) == 1
    left = FOV(1); right = left; top = left; bottom = left;
else
    left = 90; right = 90; top = 90; bottom = 90; 
end
gp = num2cell(gabor_params);
if length(gabor_params) == 7
    [x0,y0,lambda,theta,phi,sigma,rho] = gp{:};
elseif length(gabor_params) == 6
    [x0,y0,lambda,theta,phi,sigma] = gp{:};
    rho = Inf;
elseif length(gabor_params) == 4
    [lambda,theta,phi,sigma] = gp{:};
    x0 = 0; y0 = 0; rho = Inf;
else
    error('Wrong number of filter parameters!');
end
if length(grid_size) == 2
    Nx = grid_size(1);
    Ny = grid_size(2);
else
    Nx = grid_size(1); Ny = Nx;
end
X0 = 1:Nx;
Y0 = 1:Ny;
X0 = (left+right)*(X0-0.5)/Nx-left;
Y0 = (top+bottom)*(Y0-0.5)/Ny-top;
[X Y] = meshgrid(X0,Y0);
X = X-x0;
Y = Y-y0;
theta = pi*theta/180;
YY = cos(theta)*Y-sin(theta)*X;
F = cos(2*pi*(-YY/lambda+phi/360));
% Gaussian normalization:
G = exp(-((X/sigma).^2+(Y/sigma).^2)/2)/sigma/sigma/pi/2;
g0 = exp(-rho^2/2)/sigma/sigma/pi/2;
G(G<g0) = 0;
Z = G.*F;
end