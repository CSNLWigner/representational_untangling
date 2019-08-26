function rho = mixed_rho(r,mu,sigma,nu,A,kappa)
% Gaussian-rectification model
% potential /u/ distribution: Gauss(mu,sigma)
% parameters of the nonlinearity u -> r :
% threshold (nu), prefactor (A), exponent (kappa)
% mixed type firing rate /r/ distribution:
% discrete component: p0(r=0|theta,phi)  
% continuous component: rho(r|theta,phi)
if kappa == 1
    rho = exp(-(r/A+nu-mu).^2/sigma^2/2)/sigma/A/sqrt(2*pi);
else
    rho = (A*r.^(kappa-1)).^(-1/kappa).*exp(-((r/A).^(1/kappa)+nu-mu).^2/sigma^2/2)/sigma/kappa/sqrt(2*pi);
end