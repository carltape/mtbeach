function [gamma,delta,thetadc,phi] = nu2lune(nu,n,blune,brandom)
%NU2LUNE get a set of lune points for a fixed nu value
%
% INPUT
%   nu      Poisson parameter (-infinity,infinity)
%   n       number of points
% optional:
%   blune   =true for uniform spacing on the lune
%   brandom =false for regularly space points; =true for random
%
% OUTPUT
%   gamma   lune longitude [-30,30]
%   delta   lune latitude [-90,90]
% optional:
%   thetadc polar angle on the lune [0,90]
%   phi     azimuthal angle on the lune [-180,180]
%
% See examples below.
% TapeTape2013 "The classical model for moment tensors"
%
% calls phi2lamK.m, s2zeta.m, phizeta2lam.m, tp2lam.m
% 
% Carl Tape, 2022-08-30
%

if nargin==2
   blune = true;
   brandom = false;
end
if nargin==3, brandom = false; end

deg = 180/pi;

% get the two phi angles; construct vector of phi values
phitar  = nu2phi(nu);
phitar2 = wrapTo180(phitar + 180);
disp(sprintf('nu = %.2f phi1 = %.2f phi2 = %.2f',nu,phitar,phitar2));

if brandom
    ivec = round(1 + rand(n,1));
    phi = NaN(n,1);
    phi(ivec==1) = phitar;
    phi(ivec==2) = phitar2;
    
else
    % number of points on the arc between the DC and the lune boundary
    nh = floor(n/2);
    phi1vec = phitar  * ones(nh,1);
    phi2vec = phitar2 * ones(nh,1);
    % if an odd number of points is specified, insert the DC
    if rem(n,2)==0
        warning('for regular spacing, specify an odd number of points so that the DC is included');
        phi = [phi1vec ; phi2vec];
    else
        phi = [phi1vec ; phitar ; phi2vec];
    end
end

if blune    % uniform spacing on the lune
    % angular distance from DC to crack tensor
    [~,thetaK] = phi2lamK(phitar);
    
    % get a set of uniform thetadc values
    if brandom
        thetadc = thetaK * rand(n,1);
    else
        thetadc = abs( linspace(0,2*thetaK,n)' - thetaK );
    end
    
    % eigenvalues
    lam = tp2lam(thetadc,phi);
    
else        % uniform moment tensors
    % get a set of zeta values
    if brandom
        s = rand(n,1);
        zeta = s2zeta(s);
    else
        s = abs(linspace(-1,1,n)');
        zeta = s2zeta(s);
    end
    zeta = zeta*deg;    % radians to degrees
    
    % eigenvalues
    lam = phizeta2lam(phi,zeta);
    % calculate thetadc (we already have phi)
    thetadc = lam2tp(lam);
end

% lune coordinates
[gamma,delta] = lam2lune(lam);
%[gamma,delta,~,thetadc] = lam2lune(lam);

if ~brandom
    % sort by increasing lune longitude
    %[~,isort] = sortrows([gamma delta],[1 2],{'ascend' 'descend'});
    [~,isort] = sortrows(gamma,'ascend');
    if abs(gamma) < 1e-6, [~,isort] = sortrows(delta,'descend'); end
    gamma   = gamma(isort);
    delta   = delta(isort);
    thetadc = thetadc(isort);
    phi     = phi(isort);
end

%==========================================================================

if 0
    %% uniform spacing on the lune
    nu = 0.25;
    n = 7;          % an even number will return a warning
    [gamma,delta,thetadc,phi] = nu2lune(nu,n,blune);
    xnu = phi2nu(phi);
    figure; hold on; plot(gamma,delta,'.-');
    axis([-30 30 -90 90]); grid on
    for ii=1:length(gamma)
        text(gamma(ii),delta(ii),sprintf('(\\theta=%.1f, \\phi=%.1f, \\nu=%.2f)',thetadc(ii),phi(ii),xnu(ii)));
    end
    
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
        % plot
        plot(gamma,delta,'.-','markersize',12);
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
        plot(gamma,delta,'.-'); title(stlab);
        axis([-30 30 -90 90]); grid on
        for ii=1:length(gamma)
            text(gamma(ii),delta(ii),sprintf('(\\theta=%.1f, \\phi=%.1f, \\nu=%.2f)',thetadc(ii),phi(ii),xnu(ii)));
        end
    end
    
    %%
    nu0 = 0.25;
    blune = false;
    brandom = true;
    n = 1e5;
    [gamma,delta,thetadc,phi] = nu2lune(nu0,n,blune,brandom);
    lam = lune2lam(gamma,delta);
    [nu,alpha] = lam2nualpha(lam);
    % plot
    figure; hold on; plot_histo(alpha*pi/180,linspace(0,pi,38),3);
    [xalpha,falpha] = pdf_alpha(nu0);
    plot(xalpha,falpha,'r','linewidth',2);

    %% thetaK = thetadc_crack as a function of phi
    deg = 180/pi;
    phi = [-180:180];
    % eigenvalues of crack tensor
    % OPTION A: analytical (thetaK from phi2lamK.m)
    % OPTION B: numerical; from eigenvalues (see lam2lune.m)
    [lamK,thetaK] = phi2lamK(phi);
    lam1 = lamK(1,:); lam2 = lamK(2,:); lam3 = lamK(3,:);
    lammag = sqrt(lam1.^2 + lam2.^2 + lam3.^2);
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
    xlabel('v'); ylabel('w'); title({'nu2lune.m',sprintf('\\nu = (%.2f, %.2f, %.2f)',nuvec)});

end

%==========================================================================