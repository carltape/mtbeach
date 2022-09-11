function [thetadc,phi,rho] = lam2tp(lam)
%LAM2TP convert eigenvalues to theta-phi angles on the lune
%
% INPUT
%   lam     3 x n set of eigenvalue triples
%
% OUTPUT
%   thetadc n x 1 vector of polar angles (from DC), degrees [0,90]
%   phi     n x 1 vector of azimuthal angles (from top/north), degrees [-180,180]
% optional
%   rho     1 x n vector of magnitudes
%
% note: rho = sqrt(2)*M0  -----  M0 = rho/sqrt(2)
%
% See example in tp2lam.m
%
% See TapeTape2013 "The classical model for moment tensors".
% Reverse function of tp2lam.m
%
% Carl Tape, 2022-08-30
%

deg = 180/pi;

[lam,n] = lamsort(lam);
% row vectors
lam1 = lam(1,:);
lam2 = lam(2,:);
lam3 = lam(3,:);

% magnitude of lambda vector (rho of TT2012beach -- see p. 490 within text)
rho = sqrt(lam1.^2 + lam2.^2 + lam3.^2);

% TT2013 Eq. S1 (from Eq. 12)
% thetadc: the angle between the DC and the lune point
% (This arc distance will be equal to the angular distance in matrix space
% between a full moment tensor and a double couple moment tensor having the
% same frame and same magnitude.)
%thetadc = acos( cos(delta/deg) .* cos(gamma/deg) ) * deg;
thetadc = acos( (lam1 - lam3) ./ (sqrt(2)*rho) ) * deg;

% TT2013 Eq. 24
phi = atan2( lam1-2*lam2+lam3, sqrt(2)*(lam1+lam2+lam3) ) * deg;

% special cases
XEPS = 1e-6;
phi(abs(thetadc) < XEPS) = 0;
phi(abs(phi + 180) < XEPS) = 180;
phi(abs(phi) < XEPS) = 0;

% column vectors
thetadc = thetadc(:);
phi = phi(:);
