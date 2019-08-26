function p0 = mixed_p0(mu,sigma,nu)
% Gaussian-rectification model
% potential /u/ distribution: Gauss(mu,sigma)
% parameters of the nonlinearity u -> r :
% threshold (nu), prefactor (A), exponent (kappa)
% mixed type firing rate /r/ distribution:
% discrete component: p0(r=0|theta,phi)  
% continuous component: rho(r|theta,phi)
p0 = .5+.5*erf((nu-mu)/sigma/sqrt(2));
end