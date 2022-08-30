function nu = phi2nu(phi)
%PHI2NU convert azimuthal coordinate on the lune [phi] to nu [Poisson parameter]
% 
% INPUT
%   phi     azimuthal angle on eigenvalue lune [-180,180]
%
% OUTPUT
%   nu      Poisson parameter [-infty,infty]
%   
% Note: Here nu is a mathematical parameter, and there are no restrictions
%       on it, such as the allowable Poisson values of -1 <= nu <= 0.5.
%
% See Tape and Tape (2013), "The classical model for moment tensors"
%
% See examples below.
%
% Carl Tape, 2022-08-30
%

tanphi = tan(phi*pi/180);
% solve TT2013 Eq 33b for nu
nu = (1 - sqrt(2)*tanphi)./(2 + sqrt(2)*tanphi);

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
    axis([-40 40 -90 90]); grid on
end

%==========================================================================
% EXAMPLES

if 0==1
    % values from Fig. 5b of TT2013 (set bfigure = true)
    clear, clc, close all, format short
    nu = [-1:0.25:0.5 -1000]';
    % get phi values from nu
    phi = nu2phi(nu);
    % get nu values from phi
    nucheck = phi2nu(phi);
    [nu nucheck]
    
    % values from Fig. 6 of TT2013 (set bfigure = true)
    clear, clc, close all, format short
    phi = [15:30:165 -165:30:-15]';
    % get nu values from phi
    nu = phi2nu(phi);
    % get phi values from nu
    % note: the recovered phi will be correct to +/-180
    phicheck = nu2phi(nu);
    [phi phicheck]
end

%==========================================================================