function [N,K] = sd2NK(kappa,theta)
%SD2NK strike and dip angles to fault normal vector and strike vector
%
% INPUT
%   kappa   strike angle, degrees
%   theta   dip angle, degrees
%
% OUTPUT
%   N       3 x n set of fault normal vectors in SOUTH-EAST-UP convention
%   K       3 x n set of stike vectors in SOUTH-EAST-UP convention
%
% See WTape and CTape (2012) "A geometric setting for moment tensors" (TT2012).
%
% See examples in U2sdr.m
% called by sdr2U.m
%
% Carl Tape, 2012-12-01
%

nk = length(kappa);
nt = length(theta);
%ns = length(sigma);
if length(unique([nk nt])) > 2
    whos kappa theta
    error('input must have same dimension');
end
n = nk;

% NOTE: Algorithmically, it would be simpler to compute V directly from the
% expression in TT2012 Proposition 2, since this requires fewer calculations.
% (In the case here, we do not need the fault vectors at all.)
% The implementaion below is more conceptual.
% The basis is specified from the components of the north and zenith vectors.

% for north-west-up basis (TT2012)
%north = [1 0 0]'; zenith = [0 0 1]';

% for south-east-up basis (TT2013)
north = [-1 0 0]'; zenith = [0 0 1]';

% TT2012, p. 485
phi = -kappa;

% TT2012, Eq 27ab
N = NaN(3,n);
S = NaN(3,n);
for ii=1:n
    K(:,ii) = rotmat(phi(ii),3) * north;
    N(:,ii) = rotmat_gen(K(:,ii),theta(ii)) * zenith;
end

%==========================================================================
