function I = cresponse(GPARAMS,WPARAMS)
% analytical responses of circular Gabor filters to cosine wave stimuli
% I: M x N matrix containing circular Gabor filter responses
% filter response range: [-1,1]
% N: number of circular Gabor filters in the population
% M: number of stimuli in the stimulus bank
% GPARAMS: N x G matrix containing the parameters of all filters
% G = 6: number of parameters defining a single circular Gabor filter
% x0, y0, lambda_G, theta_G, phi_G, sigma
% /distances and angles are measured in degrees/ 
% x0, y0: center of the filter measured from the line of sight
% glambda: wavelength of the cosine wave component
% gtheta: orientation of the parallel stripes in the wave
% gphi: phase of the cosine at the center of the filter
% sigma: standard deviations of the circular Gaussian envelope
% WPARAMS: M x W matrix containing the parameters of all stimuli
% W: number of parameters defining a single cosine wave stimulus
% W = 5 is used in the general case:
% slambda, stheta, sphi, DC, AC
% slambda: wavelength of the cosine wave
% stheta: orientation of the grating
% sphi: phase of the cosine at the origo (the line of sight)
% DC: phase-averaged mean stimulus intensity
% AC: amplitude of the phase modulated stimulus intensity 
% simplified parameter passing is also accepted (W = 3):
% given parameters: slambda, stheta, sphi
% -> restricted intensity range: [0,1]
[N,G] = size(GPARAMS);
[M,W] = size(WPARAMS);
if G == 6
    X0 = repmat(GPARAMS(:,1)',M,1);
    Y0 = repmat(GPARAMS(:,2)',M,1);
    GLAMBDA = repmat(GPARAMS(:,3)',M,1);
    GTHETA = repmat(GPARAMS(:,4)',M,1)*pi/180;
    GPHI = repmat(GPARAMS(:,5)',M,1)*pi/180;
    SIGMA = repmat(GPARAMS(:,6)',M,1);
else
    error('Wrong number of filter parameters!');
end
if W == 5 || W == 3
    SLAMBDA = repmat(WPARAMS(:,1),1,N);
    STHETA = repmat(WPARAMS(:,2),1,N)*pi/180;
    SPHI = repmat(WPARAMS(:,3),1,N)*pi/180;
else
    error('Wrong number of wave parameters!');
end
if W == 5 
    DC = repmat(WPARAMS(:,4),1,N);
    AC = repmat(WPARAMS(:,5),1,N);
end
if W == 3
    DC = .5*ones(M,N);
    AC = .5*ones(M,N);
end
SIG2 = 2*pi*pi*SIGMA.*SIGMA;
J0 = cos(GPHI).*exp(-2*(pi*SIGMA./GLAMBDA).^2);
PHI0 = 2*pi*(sin(STHETA).*X0-cos(STHETA).*Y0)./SLAMBDA+SPHI;
K0 = exp(-SIG2.*(1./GLAMBDA./GLAMBDA+1./SLAMBDA./SLAMBDA));
KAPPA = 2*SIG2.*cos(GTHETA-STHETA)./GLAMBDA./SLAMBDA;
I1 = (AC/2).*K0.*(cos(PHI0-GPHI).*exp(KAPPA)+cos(PHI0+GPHI).*exp(-KAPPA));
I = 4*(DC.*J0+I1);
end