function [Klam,thetaK] = Kphi(phi)
%KPHI given phi angle on the lune, return eigenvalues of crack tensor K
%
% A crack tensor is a moment tensor with normal vector and slip vector
% parallel (opening/tensional crack) or opposite (closing/compressional crack).
% This function  returns the eigenvalues of a crack tensor
% given an input parameter phi. This implements Eq 39 of TapeTape2013
% "The classical model for moment tensors". (See Figure 6.)
% 
% INPUT
%   phi     azimuthal angles on the lune (phi = 0 points toward +ISO), degrees
%
% OUTPUT 
%   Klam    3 x n set of normalized eigenvalues of the crack tensors
% optional:
%   thetaK  angular distance from DC to crack tensor
%
% See TapeTape2013 "The classical model for moment tensors"
%
% See example below.
%
% Carl Tape, 2013-12-18
%

deg = 180/pi;

% ensure row vector
phi = phi(:)';

sinphi = sin(phi/deg);
cosphi = cos(phi/deg);

% TT2013 Eq 14
V = 1/sqrt(6) * [ sqrt(3) 0 -sqrt(3) ;
                  -1 2 -1 ;
                  sqrt(2) sqrt(2) sqrt(2) ];

% TT2013 Eq 39 [see Figure 6]
% rphi is the r coordinate of Lambda^K(phi) in the cylindrical coordinate
% system attached to the vwu cartesian system.
rphi = (4*sinphi.^2 + cosphi.^2).^(-1/2);
ivec = [sqrt(3)*abs(sinphi) ; -sinphi ; cosphi];
% note: V-transpose = V-inverse for this V
Klam = repmat(rphi,3,1) .* (V' * ivec);

% the arc distance from the DC to Lambda^K (see note above about rphi)
thetaK = asin(rphi) * deg;

%==========================================================================
% EXAMPLE

if 0==1
    % values from Fig. 6 of TT2013
    phi = [15:30:165 -165:30:-15]';
    Klam = Kphi(phi);
    % check
    [thetadc,phi_check] = lam2tp(Klam);
    [gamma,delta] = lam2lune(Klam);
    [phi Klam' gamma phi_check]
    
    % alternatively we could convert phi to nu, then use nualpha2lam.m
    nu = phi2nu(phi);
    alpha = [0*ones(1,6) 180*ones(1,6)];    % half on each side of lune
    Klam_check = nualpha2lam(nu,alpha);
    norm(Klam - Klam_check)
end

%==========================================================================
