function phi = nu2phi(nu)
%NU2PHI convert nu [Poisson parameter] to azimuthal coordinate on the lune [phi]
% 
% INPUT
%   nu      Poisson parameter (-infinity,infinity)
%
% OUTPUT
%   phi     azimuthal angle on eigenvalue lune [0,180]
%   
% Note: Here nu is a mathematical parameter, and there are no restrictions
%       on it, such as the allowable Poisson values of -1 <= nu <= 0.5.
%
% See TapeTape2013"The classical model for moment tensors"
% 
% See examples in phi2nu.m
%
% Carl Tape, 2022-08-30
%

tanphi = (1 - 2*nu)./(sqrt(2)*(1+nu));
% given nu, you can only get phi or phi+180
phi = 180/pi * atan(tanphi);
% force phi values to [0,180]
ineg = find(phi < 0);
phi(ineg) = phi(ineg) + 180;

bfigure = false;    
if bfigure
    zeta = 90*ones(size(phi));
    lam = phizeta2lam(phi,zeta);
    [gamma,delta] = lam2lune(lam);
    figure; hold on;
    plot(gamma,delta,'.');
    for ii=1:length(gamma)
        text(gamma(ii),delta(ii),sprintf('(\\nu=%.2f, \\phi=%.0f)',nu(ii),phi(ii)));
    end
    axis([-30 30 -90 90]); grid on
end
