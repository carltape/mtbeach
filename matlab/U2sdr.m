function [kappa,theta,sigma,K,N,S] = U2sdr(U,bdisplay)
%U2SDR converts basis U to strike-dip-rake angles
%
% INPUT
%   U           3 x 3 x n set of bases in SOUTH-EAST-UP convention
%   bdisplay    OPTIONAL (if present, will display details)
%               
% OUTPUT
%   kappa       strike angle, degrees: [0,360]
%   theta       dip angle, degrees: [0,90]
%   sigma       rake (or slip) angle, degrees: [-90,90]
% optional:
%   K           strike vector (SOUTH-EAST-UP)
%   N           normal vector (SOUTH-EAST-UP)
%   S           slip vector (SOUTH-EAST-UP)
%
% See TapeTape2012beach "A geometric setting for moment tensors"
%
% calls NS2sdr.m
% called by CMT2TT.m
%
% Carl Tape, 2012-12-01
%

global EPSVAL
EPSVAL = 1e-6;

if nargin==1, bdisplay=false; end

% U is assumed to be 3 x 3 x n
[~,~,n] = size(U);

% moment tensor orientation
% TT2012beach Section 6.3

Yrot = rotmat(45,2);

% compute candidate fault vectors
S = zeros(3,n);
N = zeros(3,n);
for ii=1:n
   V = U(:,:,ii) * Yrot;    % V = U * Yrot (TT2012beach p. 487)
   S(:,ii) = V(:,1);        % fault slip vector
   N(:,ii) = V(:,3);        % fault normal vector
end

[kappa,theta,sigma,K,N,S] = NS2sdr(N,S,bdisplay);

%==========================================================================
% EXAMPLES

if 0==1
    % single set of angles
    kappa = 320;
    theta = 10;
    sigma = 20;
    U = sdr2U(kappa,theta,sigma)
    [kappacheck,thetacheck,sigmacheck] = U2sdr(U);
    kappa,kappacheck,theta,thetacheck,sigma,sigmacheck
    
    % large set of angles
    deg = 180/pi;
    n = 1e4;
    kappa = randomvec(0,360,n);
    h = randomvec(0,1,n);
    theta = deg*acos(h);
    sigma = randomvec(-90,90,n);
    
    U = sdr2U(kappa,theta,sigma);
    [kappacheck,thetacheck,sigmacheck] = U2sdr(U);
    norm(kappa(:) - kappacheck(:)) / norm(kappa)
    norm(theta(:) - thetacheck(:)) / norm(theta)
    norm(sigma(:) - sigmacheck(:)) / norm(sigma)
end
    
%==========================================================================
