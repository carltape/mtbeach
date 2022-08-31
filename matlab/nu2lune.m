function [gamma,delta,thetadc,phi] = nu2lune(nu,n,blune,brandom)
%NU2LUNE get a set of lune points for a fixed nu value
%
% INPUT
%   nu      Poisson parameter [-infty,infty]
%   n       number of points
% optional:
%   blune   =true for uniform spacing on the lune
%   brandom =false for regularly space points; =true for random
%
% OUTPUT
%   gamma   lune longitude
%   delta   lune latitude
% optional:
%   thetadc polar angle on the lune [0,90]
%   phi     azimuthal angle on the lune [-180,180]
%
% EXAMPLE: [gamma,delta] = nu2lune(0.25,11); figure; plot(gamma,delta,'.');
% See examples below.
%
% TapeTape2013 "The classical model for moment tensors"
%
% calls Kphi.m, uniform_zeta.m, phizeta2lam.m, tp2lam.m
% 
% Carl Tape, 2022-08-30
%

if rem(n,2)==0
    warning('input number of points must be odd');
    warning(sprintf('increasing from %i to %i',n,n+1));
    n = n+1;
end

if nargin==2
   blune = true;
   brandom = false;
end
if nargin==3, brandom = false; end

deg = 180/pi;

% number of points on the arc between the DC and the lune boundary
nh = ceil(n/2);

% get the two phi angles; construct vector of phi values
phitar  = nu2phi(nu);
phitar2 = wrapTo180(phitar + 180);
phi1vec = phitar  * ones(nh,1);
phi2vec = phitar2 * ones(nh-1,1);
phi     = [phi1vec ; phi2vec];

if blune
    % angular distance from DC to crack tensor
    [~,thetaK] = Kphi(phitar);
    
    % get a set of uniform thetadc values
    if brandom
        thetadc = thetaK * rand(2*nh-1,1);
    else
        theta1vec = linspace(0,thetaK,nh)';
        theta2vec = theta1vec(2:end);
        thetadc   = [theta1vec ; theta2vec];
    end

    % eigenvalues
    lam = tp2lam(thetadc,phi);
    
else
    % get a set of zeta values
    if brandom
        zeta = uniform_zeta(2*nh-1);
    else
        zeta1vec = uniform_zeta(nh);
        zeta2vec = zeta1vec(2:end);
        zeta = [zeta1vec ; zeta2vec];
    end
    
    % eigenvalues
    lam = phizeta2lam(phi,zeta);
    % calculate thetadc (we already have phi)
    thetadc = lam2tp(lam);
end

% lune coordinates
[gamma,delta] = lam2lune(lam);
%[gamma,delta,~,thetadc] = lam2lune(lam);

% sort by increasing lune longitude
%[~,isort] = sortrows([gamma delta],[1 2],{'ascend' 'descend'});
[~,isort] = sortrows(gamma,'ascend');
if abs(gamma) < 1e-6, [~,isort] = sortrows(delta,'descend'); end
gamma   = gamma(isort);
delta   = delta(isort);
thetadc = thetadc(isort);
phi     = phi(isort);

%==========================================================================

if 0
    %% three choices of nu
    phi = [0 nu2phi(0.25) 90];
    nu = phi2nu(phi);
    n = 21;
    figure; hold on;
    for ii=1:length(nu)
        % default settings: regular spacing on lune
        [gamma,delta,thetadc,phi] = nu2lune(nu(ii),n);
        % random points (still regular spacing on lune)
        %[gamma,delta] = nu2lune(nu(ii),n,true,true);
        % spacing on the lune
        dtheta = abs(thetadc(2)-thetadc(1));

        plot(gamma,delta,'.','markersize',12);
        text(gamma(1),delta(1)-2,sprintf('(\\nu=%.2f, \\phi=%.1f, \\Delta\\theta=%.1f)',nu(ii),phi(ii),dtheta));
    end
    axis equal, axis([-30 30 -90 90]); grid on;
    xlabel('Lune longitude \gamma'); ylabel('Lune latitude \delta');
    title(sprintf('nu2lam.m with n=%i and %i \\nu values',n,length(nu)));
    
    %% highlighting the two different spacing options (blune)
    nu = 0.25;
    n = 21;
    for kk=1:2
        if kk==1, blune = false; stlab = 'based on zeta'; else blune = true; stlab = 'based on theta'; end
        [gamma,delta,thetadc,phi] = nu2lune(nu,n,blune);
        % plot
        xnu = phi2nu(phi);
        figure; hold on;
        plot(gamma,delta,'.'); title(stlab);
        axis([-30 30 -90 90]); grid on
        for ii=1:length(gamma)
            text(gamma(ii),delta(ii),sprintf('(\\theta=%.1f, \\phi=%.1f, \\nu=%.2f)',thetadc(ii),phi(ii),xnu(ii)));
        end
    end

    %% thetaK = thetadc_crack as a function of phi
    deg = 180/pi;
    phi = [-180:180];
    % eigenvalues of crack tensor
    [Klam,thetaK] = Kphi(phi);
    % OPTION A: analytical (thetaK from Kphi.m)
    % OPTION B: numerical; from eigenvalues
    lam1 = Klam(1,:); lam2 = Klam(2,:); lam3 = Klam(3,:);
    lammag = sqrt(lam1.^2 + lam2.^2 + lam3.^2);
    % see lam2lune.m
    thetaK_num = acos( (lam1 - lam3) ./ (sqrt(2)*lammag) ) * deg;
    % plot
    figure; hold on;
    plot(phi,thetaK_num,'b');
    plot(phi,thetaK,'r--');
    axis([-180 180 0 90]); set(gca,'xtick',[-180:30:180]); grid on;
    xlabel('phi, degrees'); ylabel('thetadc, degrees');

    %% what does a sector of dphi look like on the vw-rectangle?
    n = 99;
    nuvec = [0.1 0.25 0.4];     % TT2013 Fig. 13
    figure; hold on;
    for ii=1:length(nuvec)
        [gamma,delta] = nu2lune(nuvec(ii),n);
        [v,w] = lune2rect(gamma,delta);
        plot(v,w);
    end
    axis equal, axis([-1/3 1/3 -3*pi/8 3*pi/8]); grid on
    xlabel('v'); ylabel('w');

end

%==========================================================================