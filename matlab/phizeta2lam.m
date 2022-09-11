function lam = phizeta2lam(phi,zeta)
%PHIZETA2LAM converts phi and zeta to unit lambda vector
%
% INPUT
%   phi     n x 1 vector of azimuthal angles on the lune, degrees [-180,180]
%   zeta    n x 1 vector of crack fraction within CDC model [0,90]
%
%  OUTPUT
%   lam     3 x n set of normalized eigenvalue triples, sorted lam1 >= lam2 >= lam3
%
% Reverse function is lam2phizeta.m
% See TapeTape2013 "The classical model for moment tensors"
% 
% Carl Tape, 2013-01-08
%

bdisplay = false;

deg = 180/pi;

% row vectors
phi = phi(:)';
zeta = zeta(:)';

cosp = cos(phi/deg);
sinp = sin(phi/deg);
cosz = cos(zeta/deg);
sinz = sin(zeta/deg);

% TapeTape2013
D  = 1/sqrt(2)*[0 0 1 ; 0 0 0 ; 1 0 0];         % Eq 42
rk = (4*sinp.^2 + cosp.^2).^(-1/2);             % Eq 40

% calculate moment tensor (with fixed but arbitrary orientation)
n = length(phi);
Mh = NaN(3,3,n);
for ii=1:n
   % Eq 41
   K = rk(ii)/sqrt(3)*[cosp(ii) - sqrt(2)*sinp(ii) 0 0
                        0 cosp(ii) - sqrt(2)*sinp(ii) 0
                        0 0 cosp(ii) + 2*sqrt(2)*sinp(ii)];
   % Eq 43
   Mh(:,:,ii) = cosz(ii)*D + sinz(ii)*K;
   if bdisplay
       disp(sprintf('phizeta2lam.m: event %i/%i',ii,n));
       phi, zeta
       cosp, sinp, cosz, sinz
       D, rk, K
       Mh(:,:,ii)
   end
end
% convert to 6 x n array
Mhvec = Mvec2Mmat(Mh,0);
% calculate eigenvalues (note: basis U is irrelevant)
lam = CMTdecom(Mhvec);

% unit norm
lam1 = lam(1,:);
lam2 = lam(2,:);
lam3 = lam(3,:);
lam = [lam1 ; lam2 ; lam3];
rho = sqrt( lam1.^2 + lam2.^2 + lam3.^2 );
lam = lam ./ repmat(rho,3,1);

%==========================================================================
% EXAMPLE

if 0==1
    % grid of phi-zeta values
    pvec = -175:5:175;
    zvec = 5:5:85;
    [phi,zeta] = meshgrid(pvec,zvec);
    phi = phi(:); zeta = zeta(:);
    lam = phizeta2lam(phi,zeta);

    % convert back to phi-zeta, then check
    [phicheck,zetacheck] = lam2phizeta(lam);
    norm(phi-phicheck)
    norm(zeta-zetacheck)
end

%==========================================================================
