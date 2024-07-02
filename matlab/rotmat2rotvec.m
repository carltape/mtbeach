%
% function [u,umag] = rotmat2rotvec(U)
% Carl Tape, 11-March-2011
%
% Convert a rotation matrix to a rotation vector.
%
% See also v2U.m for reverse program.
% 
% calls Wperp.m, fangle.m
% called by Ueigvec.m
%

function [u,umag] = rotmat2rotvec(U,EPSVAL)

% for checking GCMT U computed from trend and plunge angles accurate to
% single digits, we need to use EPSVAL = 1e-2 (instead of 1e-6)
if nargin==1, EPSVAL = 1e-6; end

deg = 180/pi;

% for a rotation matrix there will be one eigenvector with
% eigenvalue =1, and two eigenvectors with complex eigenvalues
[V,D] = eig(U);

ipick = find( abs(diag(D) - 1) <= EPSVAL);

if length(ipick) == 1
    v = V(:,ipick);     % rotation axis candidate
    
elseif length(ipick) == 3
    u = zeros(3,1);
    umag = 0;
    disp('rotmat2rotvec.m: U is the identity matrix');
    return
    
else
    % typically this error can be avoided by lowering EPSVAL
    ipick, V, D, U
    error('rotmat2rotvec.m: U must have exactly one eigenvector with eigenvalue =1');
end

% wperp for candidate rotation axis
w = Wperp(v);

% candidate rotation angle
rotangle_candidate = fangle_signed(w, U*w, v);

% rotation axis
u = sin( rotangle_candidate/2 / deg) * v;

% rotation angle
umag = fangle(w, U*w);
        
%==========================================================================
