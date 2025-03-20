function plotMT_nualpha(nu,alpha)
%PLOTMT_NUALPHA plot histograms of moment tensor parameters nu,alpha
%
% INPUT
%   nu      Poisson parameter from classical model (-infinity, infinity)
%   alpha   angle between fault normal and slip vector, degrees [0,180]
%
% See TapeTape2013 "The classical model for moment tensors"
% See also plotMT_nualpha.m, lam2nualpha.m, nualpha2lam.m
% See run_uniformMT.m for examples.
%
% Carl Tape, 2018-08-06
%

PLOT_UNIFORM_CURVES = true;
PLOT_AS_PDFS = true;
buse_greek = true;             % =true unless you have an issue rendering latex fonts

if PLOT_AS_PDFS, itype = 3; else itype = 2; end

edges_alpha = [0:5:180];
% get the nu edges from the plotting limits of the analytical function
[xplot,uplot] = pdf_nu;
NUMIN = min(xplot);
NUMAX = max(xplot);
nbinnu = 51;
edges_nu = linspace(NUMIN,NUMAX,nbinnu);
%edges_nu = linspace(-1.001,0.501,nbinnu);
%edges_nu = [-1:0.1:0.5];

% angles in radians, since we are dealing with PDFs
deg = 180/pi;
if PLOT_AS_PDFS
   alpha = alpha/deg;
   edges_alpha = edges_alpha/deg;
   xtag = 'radians';
else
   xtag = 'degrees';
end

%fsizex = 14;

figure; nr=2; nc=1;

% plots for NU (Poisson parameter)
subplot(nr,nc,1); hold on;
plot_histo(nu,edges_nu,itype);
if PLOT_UNIFORM_CURVES
    plot(xplot,uplot,'r','linewidth',2);
    % maximum of the plotted distribution
    umax = sqrt(2) ./ (2/3 .* (atan((1-3*NUMIN)/sqrt(2)) +  atan((-1+3*NUMAX)/sqrt(2))) );
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
    nu0 = [];
    [xplot,uplot] = pdf_alpha(nu0);
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

%==========================================================================
