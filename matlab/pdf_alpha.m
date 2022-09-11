function [alpha,falpha] = pdf_alpha(nu,n)
%PDF_ALPHA probability density function for alpha angle (radians) for uniform moment tensors
%
% INPUT
%   nu      =[] to indicate that all moment tensors are considered
%           =scalar to indicate that nu is fixed
% optional
%   n       number of points
%
% OUTPUT
%   alpha   angle between slip and normal vector, in RADIANS [0,pi]
%   falpha  f(alpha)
%
% EXAMPLE: [alpha,falpha] = pdf_alpha(0.25); plot(alpha,falpha);
%
% See TapeTape2013 "The classical model for moment tensors"
% called by plotMT_nualpha.m
%
% Carl Tape, 2022-09-09
%

if nargin==1, n = 501; end

if isempty(nu)  % all moment tensors
    
    alpha = linspace(0,pi,n)';
    sina = sin(alpha);
    cosa = cos(alpha);
    % AlphaZetaThetaForCarloTOSS.pdf
    falpha = 27 * sina.^3 ./ ( 2*(3 + cosa.^2).^(5/2) );
    
else            % all moment tensors with fixed nu
    
    % the equation is defined for alpha = [0,90]
    nh = ceil(n/2);
    alpha = linspace(0,pi/2,nh)';
    sina = sin(alpha);
    cos2a = cos(2*alpha);
    
    phi_deg = nu2phi(nu);       % degrees
    phi = phi_deg * pi/180;     % radians
    disp(sprintf('input nu = %.2f, phi = %.2f deg = %.4f rad',nu,phi_deg,phi));
    % PDFforAlphaTOSS.pdf
    cos2p = cos(2*phi);
    sinp  = sin(phi);
    exp1 = 5 - 3*cos2p;
    exp2 = 4 + cos2a - 3*cos2p;
    falpha = 54 * sqrt(exp1) .* sina.^3 .* sinp.^4 ./ exp2.^(5/2);
    
    % double it to cover alpha = [90,180]
    a1 = alpha;
    a2 = linspace(pi/2,pi,nh)'; a2(1) = [];
    f1 = falpha;
    f2 = flipud(falpha(1:end-1));
    falpha = [f1 ; f2];
    alpha = [a1 ; a2];

    % since we doubled the range of the PDF, we need to cut the amplitude in half
    falpha = falpha/2;
end

%==========================================================================
% EXAMPLES

if 0
    % all moment tensors
    figure; hold on;
    nu = [];
    [alpha,falpha] = pdf_alpha(nu); plot(alpha,falpha,'r');
    nu = 0.25;
    [alpha,falpha] = pdf_alpha(nu); plot(alpha,falpha,'r--');
    legend('all MTs',sprintf('MTs with \\nu = %.2f',nu));
    axis equal, axis([0 pi -0.1 3]); grid on; xlabel('\alpha, radians');
    
    %% vary nu/phi
    phivec = [7 10 15 20 30 45 90];
    nuvec = phi2nu(phivec);
    nplot = length(nuvec);
    xplot = linspace(pi/2,2,nplot);
    figure; hold on;
    % show the PDF for all MTs
    [alpha,falpha] = pdf_alpha([]); plot(alpha,falpha,'k','linewidth',4);
    % loop over different nu/phi
    for ii=1:nplot
        [alpha,falpha] = pdf_alpha(nuvec(ii));
        plot(alpha,falpha);
        stlab = sprintf('\\phi = %.1f, \\nu = %.2f',phivec(ii),nuvec(ii));
        yplot = interp1(alpha,falpha,xplot(ii));
        plot(pi/2,max(falpha),'ko','markerfacecolor','k');
        plot(xplot(ii),yplot,'ko');
        text(xplot(ii) + 0.1,yplot,stlab);
    end
    axis equal, axis([0 pi -0.1 4]); grid on;
    xlabel('\alpha, radians');
    
end

%==========================================================================
