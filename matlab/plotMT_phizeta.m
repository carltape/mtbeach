function plotMT_phizeta(phi,zeta)
%PLOTMT_PHIZETA plot histograms of moment tensor parameters phi,zeta
%
% INPUT
%   phi     azimuthal angle on the lune, degrees [-180,180]
%   zeta    crack fraction within CDC model [0,90]
% 
% See TapeTape2013 "The classical model for moment tensors"
% See also plotMT_nualpha.m, lam2phizeta.m, phizeta2lam.m
% See run_uniformMT.m for examples.
%
% Carl Tape, 2018-08-06
%

PLOT_UNIFORM_CURVES = true;
PLOT_AS_PDFS = true;
buse_greek = true;  % =true unless you have an issue rendering latex fonts

if PLOT_AS_PDFS, itype = 3; else itype = 2; end

edges_phi   = [-180:10:180];
edges_zeta  = [0:2:90];

% angles in radians, since we are dealing with PDFs
deg = 180/pi;
if PLOT_AS_PDFS
   phi = phi/deg;
   zeta = zeta/deg;
   edges_phi = edges_phi/deg;
   edges_zeta = edges_zeta/deg;
   xtag = 'radians';
else
   xtag = 'degrees';
end

figure; nr=2; nc=1;

% plots for PHI (azimuthal angle on the lune; similar to nu)
subplot(nr,nc,1); hold on;
plot_histo(phi,edges_phi,itype);
if PLOT_UNIFORM_CURVES
    [xplot,uplot] = pdf_phi;
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
subplot(nr,nc,2); hold on;
plot_histo(zeta,edges_zeta,itype);
if PLOT_UNIFORM_CURVES
    [xplot,uplot] = pdf_zeta(false);
    plot(xplot,uplot,'r','linewidth',2);
    [xplot,uplot] = pdf_zeta(true);
    plot(xplot,uplot,'r--','linewidth',2);
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
