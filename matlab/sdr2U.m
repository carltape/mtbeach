function U = sdr2U(kappa,theta,sigma)
%SDR2U strike-dip-rake to rotation matrix U 
%
% INPUT
%   kappa   strike angle, degrees
%   theta   dip angle, degrees
%   sigma   slip (or rake) angle, degrees
%
% OUTPUT
%   U       3 x 3 x n set of bases in SOUTH-EAST-UP convention
%
% See TapeTape2012 "A geometric setting for moment tensors"
%
% See examples in U2sdr.m
%
% calls sd2NK.m
% called by TT2CMT.m
%
% Carl Tape, 2012-12-01
%

nk = length(kappa);
nt = length(theta);
ns = length(sigma);
if length(unique([nk nt ns])) > 2
    whos kappa theta sigma
    error('input must have same dimension');
end
n = nk;

for ii=1:n
    if theta(ii)==0
        warning('%i/%i input fault is horizontal, so strike angle (%.1f) is undefined',ii,n,kappa(ii));
        disp(sprintf('resetting slip angle (%.1f) to 0',sigma(ii)));
        sigma(ii) = 0;
    end
end

[N,K] = sd2NK(kappa,theta);

% TT2012beach Eq 27c
S = NaN(3,n);
for ii=1:n
    S(:,ii) = rotmat_gen(N(:,ii),sigma(ii)) * K(:,ii);
end

% TT2012beach Eq 28 (or Proposition 2)
U = NaN(3,3,n);
Yrot = rotmat(-45,2);
for ii=1:n
    V = [S(:,ii) cross(N(:,ii),S(:,ii)) N(:,ii)];
    Ux = V*Yrot;
    U(:,:,ii) = Ux;
end

%==========================================================================
