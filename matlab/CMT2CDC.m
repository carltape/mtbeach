function [K1,D1,K2,D2,KU1,DU1,KU2,DU2,rho,zeta,U] = CMT2CDC(Mmat,bdisplay)
%CMT2CDC moment tensor to basic crack-plus-double couple (two possibilities)
%
% This function decomposes an input moment tensor into the sum of a 
% double couple and a crack tensor with its plane coinciding with a DC plane.
%   M = K1 + D1 = K2 + D2
%
% NOTE: ALL INPUT AND OUTPUT MATICES ARE 3 X 3 SYMMETRIC MATRICES.
% NOTE: THE BASIS FOR THE OUTPUT IS THE SAME AS THE BASIS FOR THE INPUT.
%
% INPUT
%   Mmat    moment tensor
%
% OUTPUT
%   K1      crack tensor for 1st decomposition
%   D1      double couple tensor for 1st decomposition
%   K2      crack tensor for 2nd decomposition
%   D2      double couple tensor for 2nd decomposition
%   KU1     
%   DU2     
%   KU2     
%   DU2     
%
% References:
%   Tape and Tape (2013 GJI), Appendix A (see also TT2013AppA.m)
%   Minson et al. (2007 JGR)
%   Aki and Richards (1980)
% Applications:
%   Alvizuri and Tape (2018 SRL)
%

if nargin==1, bdisplay=false; end
if bdisplay, disp('entering CMT2DC.m'); end

deg = 180/pi;

M = Mvec2Mmat(Mmat,0);

% eigenvalues and basis
% note: U has same basis as input M
[lams,U] = CMTdecom(M,1);
lam1 = lams(1);
lam2 = lams(2);
lam3 = lams(3);

% normalized eigenvalues
rho = sqrt(sum(lams.^2));
lamh = lams / rho;
lam1h = lamh(1);
lam2h = lamh(2);
lam3h = lamh(3);

% unit double couple tensor
D = 1/sqrt(2) * [0 0 1 ; 0 0 0 ; 1 0 0];

if bdisplay
    Mmat,U,det(U),lams,D
end

% unit crack tensor
K0 = [lam2 0 0 ; 0 lam2 0 ; 0 0 lam1-lam2+lam3 ];
Knorm = sqrt(sum(diag(K0).^2));     % norm(K0(:))
K = K0 / Knorm;

alpha = acos( (lam1-2*lam2+lam3) / (lam1-lam3) )*deg;

zeta = acos( sqrt(2*(lam1-lam2)*(lam2-lam3)) / rho )*deg;
cosz = cos(zeta/deg);
sinz = sin(zeta/deg);

% note: Eqs 36 and 56 are the same
N1 = U*rotmat(-alpha/2,2)*[1 0 0]';
N2 = U*rotmat(alpha/2,2)*[1 0 0]';

% unit double couple tensors
chi = (180 - alpha)/2;       % Eq 51
U1 = U*rotmat(chi,2);
U2 = U*rotmat(180,1)*rotmat(chi,2);
DU1 = U1*D*U1';
DU2 = U2*D*U2';

% unit crack tensors
KU1 = U1*K*U1';
KU2 = U2*K*U2';

D1 = rho*cosz*DU1;
K1 = rho*sinz*KU1;
D2 = rho*cosz*DU2;
K2 = rho*sinz*KU2;

% check
if bdisplay
    disp('check from CMT2CDC.m:');
    Mmat
    Mcheck1 = D1 + K1
    Mcheck2 = D2 + K2
    norm(Mmat(:) - Mcheck1(:)) ./ norm(Mmat(:))
    norm(Mmat(:) - Mcheck2(:)) ./ norm(Mmat(:))
end

%==========================================================================

if 0==1
    % EXAMPLE 1: see TT2013AppA.m
    
    % EXAMPLE 2: Minson et al. (2007), Eq. A13
    % note 1: This is a contrived example in the sense that one can see
    % from inspection that one of the decompositions will be the diagonal
    % (3,1,1) [crack] and the off-diagonal matrix [DC]. Also the crack
    % tensor is selected such that nu = 0.25, which is the constraint
    % preferred by Minson ("In our inversions, we assume a Poisson solid
    % (nu = 0.25)...")
    % note 2: The basis is assumed to be north-east-down, but this choice
    % does not matter for the sake of the CDC decomposition.
    clear, close all, clc
    Mmat = [3 1 0 ; 1 1 0 ; 0 0 1];
    [K1,D1,K2,D2,KU1,DU1,KU2,DU2,rho,zeta,U] = CMT2CDC(Mmat,true);
    % display in same order as in Minson
    K2,D2,K1,D1
end

%==========================================================================
