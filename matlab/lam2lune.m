function [gamma,delta,M0,thetadc,lamdev,lamiso] = lam2lune(lam)
%LAM2LUNE convert eigenvalues to lune coordinates (gamma, delta, M0)
%
% INPUT
%   lam         3 x n set of eigenvalues for a set of moment tensors
%               
% OUTPUT
%   gamma       angle from DC meridian to lune point (-30 <= gamma <= 30)
%   delta       angle from deviatoric plane to lune point (-90 <= delta <= 90)
%   M0          seismic moment, M0 = ||lam|| / sqrt(2)
%   thetadc     angle from DC to lune point (0 <= thetadc <= 90)
%                  (see also zeta angle in lam2phizeta.m)
%   lamdev      eigenvalues of deviatoric component
%   lamiso      eigenvalues of isotropic component
%
% Reverse program for lune2lam.m
% See also CMT2all.m
%
% See TapeTape2012beach "A geometric setting for moment tensors".
%
% Carl Tape, 2011-04-01
%

deg = 180/pi;

[lam,n] = lamsort(lam);
% row vectors
lam1 = lam(1,:);
lam2 = lam(2,:);
lam3 = lam(3,:);

% magnitude of lambda vector (rho of TT2012beach -- see p. 490 within text)
rho = sqrt(lam1.^2 + lam2.^2 + lam3.^2);

% seismic moment
M0 = rho / sqrt(2);

% TT2012beach Eq. 21a (and 23)
% numerical safety 1: if trace(M) = 0, delta = 0
% numerical safety 2: is abs(bdot) > 1, adjust bdot to +1 or -1
delta = zeros(1,n);         % initialized to delta=0
idev = find(sum(lam) ~= 0);
bdot = (lam1 + lam2 + lam3) ./ (sqrt(3)*rho);
bdot(bdot > 1) = 1; bdot(bdot <-1) = -1;
delta(idev) = 90 - acos(bdot(idev)) * deg;

% TT2012beach Eq. 21a
gamma = atan((-lam1 + 2*lam2 - lam3) ./ (sqrt(3)*(lam1 - lam3))) * deg;
% note: we set gamma=0 for (1,1,1) and (-1,-1,-1)
XEPS = 1e-6;
biso = find(abs(lam1-lam3) < XEPS);
%biso = lam1==lam3;
gamma(biso) = 0;

% TT2013 Eq. S1 (from Eq. 12)
% thetadc: the angle between the DC and the lune point
% (This arc distance will be equal to the angular distance in matrix space
% between a full moment tensor and a double couple moment tensor having the
% same frame and same magnitude.)
%thetadc = acos( cos(delta/deg) .* cos(gamma/deg) ) * deg;
thetadc = acos( (lam1 - lam3) ./ (sqrt(2)*rho) ) * deg;

% extra output
trM = sum(lam); 
lamiso = repmat(1/3*trM,3,1);
lamdev = lam - lamiso;

% column vectors
delta = delta(:);
gamma = gamma(:);
M0 = M0(:);
thetadc = thetadc(:);

%==========================================================================
% EXAMPLE

if 0==1
    clear, close all, clc
    % generate moment tensor source types (on the lune) having fixed magnitude
    gvec = linspace(-30,30,100);
    bvec = linspace(-89,89,100);
    [G,B] = meshgrid(gvec,bvec);
    gamma0 = G(:);
    delta0 = B(:);
    M00 = 1e16*ones(length(gamma0),1);
    lam = lune2lam(gamma0,delta0,M00);
    
    [gamma,delta,M0,thetadc] = lam2lune(lam);
    
    figure; nr=3; nc=1;
    subplot(nr,nc,1); plot(gamma-gamma0,'.'); title('gamma residual');
    subplot(nr,nc,2); plot(delta-delta0,'.'); title('delta residual');
    subplot(nr,nc,3); plot(M0-M00,'.'); title('M0 residual');
    
    % thetadc on the lune
    figure; scatter(gamma,delta,4^2,thetadc,'filled');
    caxis([0 90]); colorbar;
    xlabel('gamma, deg'); ylabel('delta, deg');
    title('angle between DC and MT point');
end

%==========================================================================

