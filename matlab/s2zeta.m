function zeta = s2zeta(s)
%S2ZETA convert s to zeta angle
%
% INPUT
%   s       parameter [0,1]
%
% OUTPUT
%   zeta    zeta angle, in RADIANS [0,pi/2]
%
% This is the reverse function of zeta2s.m
%
% Carl Tape, 2022-09-09
%

s = round2one(s);
exp1 = 2 * cos( 1/3 * (pi + acos(s)) );
exp1 = round2one(exp1);
zeta = asin(exp1);

%--------------------------------------------------------------------------
function x = round2one(x)

% round values near 1 to avoid returning complex numbers
ETRSH = 1e-6;
bround = and(x > 1, x < 1+ETRSH);
x(bround) = 1;

%==========================================================================
% EXAMPLE

if 0
   %% check that these functions are inverses
   s = linspace(0,1);
   zeta = s2zeta(s);
   s_check = zeta2s(zeta);
   figure; plot(s,s_check,'.'); grid on;
   norm(s - s_check) / norm(s)
   
    %% uniform random MTs for a fixed nu/phi
    n = 1e6;
    s = rand(n,1);
    zeta = s2zeta(s);
    
    deg = 180/pi;
    dedge = [0:2:90];
    figure; hold on;
    plot_histo(zeta,dedge/deg,3); xlabel('\zeta, radians');
    % analytical curves
    [zetac,fzeta] = pdf_zeta(true);      % for fixed phi/nu
    plot(zetac,fzeta,'r','linewidth',2);
    [zetac,fzeta] = pdf_zeta(false);     % for all moment tensors
    plot(zetac,fzeta,'r--','linewidth',2);
    %figure; plot_histo(zeta*deg,dedge); xlabel('\zeta, degrees');
    
    %% compare with a swath of uniform MTs on the lune
    n = 1e6;
    [M,v,w,kappa,sigma,h,lam] = uniformMT(n);
    [phi,zeta] = lam2phizeta(lam);
    nu = phi2nu(phi);
    % subset of nu
    NUTARGET = 0.25;
    NUEPSILON = 0.01;
    isub = find( abs(nu-NUTARGET) < NUEPSILON );
    nu1 = NUTARGET - NUEPSILON;
    nu2 = NUTARGET + NUEPSILON;
    zetaX = zeta(isub);
    % plot
    figure; plot_histo(zetaX,[0:2:90]); xlabel('\zeta, degrees');
    title(sprintf('subset for \\nu between %.2f and %.2f',nu1,nu2));
    
end

%==========================================================================
