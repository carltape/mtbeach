function plotMT_cdc(nu,alpha,phi,zeta)
%PLOTMT_CDC plot histograms of moment tensor parameters nu,alpha,phi,zeta
%
% CDC = crack-plus-double couple model for seismic moment tensors
% See TapeTape2013 "The classical model for moment tensors"
%
% INPUT
%   nu      Poisson parameter from classical model (-infinity, infinity)
%   alpha   angle between fault normal and slip vector, degrees [0,180]
%   phi     azimuthal angle on the lune, degrees [-180,180]
%   zeta    crack fraction within CDC model [0,90]
% 
% See also lam2nualpha.m, lam2phizeta.m
% See run_uniformMT.m for examples.
%
% Carl Tape, 2018-08-06
%

PLOT_UNIFORM_CURVES = true;
PLOT_AS_PDFS = true;
buse_greek = false;             % =true unless you have an issue rendering latex fonts

if PLOT_AS_PDFS, itype = 3; else itype = 2; end

edges_alpha = [0:5:180];
edges_phi   = [-180:10:180];
edges_zeta  = [0:5:90];

% edges for nu are tricky, since nu can range from -infinity to infinity
% Note: for real materials, nu varies from -1 to 0.5
NUMAX = 5;
NUMIN = -NUMAX;
nbinnu = 51;
edges_nu = linspace(NUMIN,NUMAX,nbinnu);
%edges_nu = linspace(-1.001,0.501,nbinnu);
%edges_nu = [-1:0.1:0.5];

edges_alpha = [0:5:180];
edges_phi   = [-180:10:180];
edges_zeta  = [0:2:90];

% angles in radians, since we are dealing with PDFs
deg = 180/pi;
if PLOT_AS_PDFS
   alpha = alpha/deg;
   phi = phi/deg;
   zeta = zeta/deg;
   edges_alpha = edges_alpha/deg;
   edges_phi = edges_phi/deg;
   edges_zeta = edges_zeta/deg;
   xtag = 'radians';
else
   xtag = 'degrees';
end

%fsizex = 14;

figure; nr=2; nc=2;

% plots for NU (Poisson parameter)
subplot(nr,nc,1); hold on;
plot_histo(nu,edges_nu,itype);
if PLOT_UNIFORM_CURVES
    npt = 500;
    xplot = linspace(NUMIN,NUMAX,npt);
    if 1==1
        % Untitled-3.nb
        qplot = 1 - 2*xplot + 3*xplot.^2;
        K = sqrt(2) / (atan((1-3*NUMIN)/sqrt(2)) +  atan((-1+3*NUMAX)/sqrt(2)));
        uplot = K ./ qplot;
        umax = sqrt(2) ./ (2/3 .* (atan((1-3*NUMIN)/sqrt(2)) +  atan((-1+3*NUMAX)/sqrt(2))) );
    else
        % NuForCarloTOSS.pdf
        qplot = 1./ ((xplot - 1/3).^2 + 2/9);
        dx = xplot(2) - xplot(1);
        C = 1 / (sum(qplot) * dx);    % crude integration (should be analytical)
        uplot = C*qplot;
        umax = C*9/2;
    end
    plot(xplot,uplot,'r','linewidth',2);
    % maximum of distribution
    plot(1/3,umax,'ro','markerfacecolor','k','markersize',6);
end
if buse_greek
    xlabel('\nu, Poisson ratio');
else
    xlabel('nu, Poisson ratio');
end
box on;    % note sure why this is needed, but it is
%title('(a)','Units','normalized','Position',[-0.1 1.05],'HorizontalAlignment','left','fontsize',fsizex-2);

% plots for ALPHA (angle between slip and normal vector)
subplot(nr,nc,2); hold on;
plot_histo(alpha,edges_alpha,itype);
if PLOT_UNIFORM_CURVES
    xplot = linspace(0,pi,npt);
    sina = sin(xplot);
    cosa = cos(xplot);
    % AlphaZetaThetaForCarloTOSS.pdf
    uplot = 27*sina.^3 ./ (2*(3 + cosa.^2).^(5/2));
    plot(xplot,uplot,'r','linewidth',2);
    % maximum of distribution
    umax = sqrt(3)/2;
    plot(pi/2,umax,'ro','markerfacecolor','k','markersize',6);
end
if buse_greek
    xlabel(sprintf('\\alpha = \\angle (N, S), %s',xtag));
else
    xlabel(sprintf('alpha = angle(N, S), %s',xtag));
end
%set(gca,'xtick',[0:30:180]);
%title('(b)','Units','normalized','Position',[-0.1 1.05],'HorizontalAlignment','left','fontsize',fsizex-2);

% plots for PHI (azimuth spoke on the lune; similar to nu)
subplot(nr,nc,3); hold on;
plot_histo(phi,edges_phi,itype);
if PLOT_UNIFORM_CURVES
    xplot = linspace(-pi,pi,npt);
    uplot = 2 ./ (5*pi - 3*pi*cos(2*xplot));
    plot(xplot,uplot,'r','linewidth',2);
end
if buse_greek
    xlabel(sprintf('\\phi, azimuth on lune, %s',xtag));
else
    xlabel(sprintf('phi, azimuth on lune, %s',xtag));
end
%set(gca,'xtick',[-180:60:180]);
%title('(c)','Units','normalized','Position',[-0.1 1.05],'HorizontalAlignment','left','fontsize',fsizex-2);

% plots for ZETA (crack fraction for the crack-plus-double-couple model)
subplot(nr,nc,4); hold on;
plot_histo(zeta,edges_zeta,itype);
if PLOT_UNIFORM_CURVES
    xplot = linspace(0,pi/2,npt);
    cosz = cos(xplot);
    sinz = sin(xplot);
    % PDFsForZetaAndTheta.pdf
    uplot = 4*cosz.^3.*sinz;
    plot(xplot,uplot,'r','linewidth',2);
    % maximum of distribution
    umax = 3*sqrt(3)/4;
    plot(pi/6,umax,'ro','markerfacecolor','k','markersize',6);
end
if buse_greek
    xlabel(sprintf('\\zeta, crack fraction, %s',xtag));
else
    xlabel(sprintf('zeta, crack fraction, %s',xtag));
end
%title('(d)','Units','normalized','Position',[-0.1 1.05],'HorizontalAlignment','left','fontsize',fsizex-2);

%==========================================================================
