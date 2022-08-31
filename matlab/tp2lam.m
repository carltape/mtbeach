function lam = tp2lam(thetadc,phi,rho)
%TP2LAM convert theta-phi angles on the lune to eigenvalues
%
% INPUT
%   thetadc n x 1 vector of polar angles (from DC), degrees [0,90]
%   phi     n x 1 vector of azimuthal angles (from top/north), degrees [-180,180]
% optional
%   rho     1 x n vector of magnitudes
%
% OUTPUT
%   lam     3 x n set of eigenvalue triples
%
% note: rho = sqrt(2)*M0  -----  M0 = rho/sqrt(2)
%
% See TapeTape2013 "The classical model for moment tensors".
% Reverse function of lam2tp.m
%
% Carl Tape, 2022-08-30
%

deg = 180/pi;

% convert to row vectors in radians
thetadc = thetadc(:)' / deg;
phi = phi(:)' / deg;
n = length(thetadc);

% default: unit sphere
if nargin==2, rho = ones(size(thetadc)); end

% Eq 14 of TT2013
V = 1/sqrt(6) * [ sqrt(3) 0 -sqrt(3) ;
                  -1 2 -1 ;
                  sqrt(2) sqrt(2) sqrt(2) ];

% Eq 17 of TT2013
F = NaN(3,n);
F(1,:) = rho.*sin(thetadc).*cos(phi);
F(2,:) = rho.*sin(thetadc).*sin(phi);
F(3,:) = rho.*cos(thetadc);

% Eq 21 of TT2013
Yrot = rotmat(90,2);
Xrot = rotmat(180,1);
% note that for this V, the inverse is the transpose
lam  = V' * Xrot * Yrot * F;

%==========================================================================
% EXAMPLE

if 0==1
    % generate grid of thetadc-phi values -- out to 30 deg only
    % note: need to check cases where the input dots (Lambdas) are OUTSIDE the fundamental lune
    dphi = 45; dtheta = 10; THETADCMAX = 30;
    phivec = [(-180+dphi):dphi:180];
    thetadcvec = [dtheta:dtheta:THETADCMAX];
    [X,Y] = meshgrid(phivec,thetadcvec);
    % add the DC and +/- ISO
    phi = [X(:) ; 0 ; 0 ; 180];
    thetadc = [Y(:) ; 0 ; 90 ; 90];

    lam = tp2lam(thetadc,phi);
    [thetadc_check,phi_check,rho_check] = lam2tp(lam);
    % check
    [thetadc thetadc_check phi phi_check]

    % plot
    [gamma,delta] = lam2lune(lam);
    figure; hold on;
    plot(gamma,delta,'.');
    for ii=1:length(phi), text(gamma(ii),delta(ii),sprintf('(%.1f, %.1f)',thetadc(ii),phi(ii))); end
    axis([-30 30 -90 90]); grid on;
end

%==========================================================================