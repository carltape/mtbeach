function [p,t,pcum] = omegapdf(x,bfigure)
%OMEGAPDF probability density function for random moment tensors
%
% INPUT
%   x       option A: an integer number of points for discretization OR
%           option B: an input vector of points at which to evaluate the PDF
%               INPUT IS IN RADIANS, NOT DEGREES
%
% OUTPUT
%   p       PDF discretized based on input choice
%   t       discretization of omega (IN RADIANS)
%   pcum    CDF discretized based on input choice
%
% See also plot_omega.m, kaganpdf.m, omegadcpdf.m, run_unformMT.m
%
% EXAMPLES (see longer example below):
%   x = linspace(0,pi,181); p = omegapdf(x,true);
%   [p,t] = omegapdf(10,true);
%   [p,~,pcum] = omegapdf(x); figure; plot(x,pcum);
%

if nargin==1, bfigure=false; end

% key points on curve
t0 = 0;
t3 = pi;

% evaluate the PDF at the bin centers
if length(x)==1
   n = round(x); 
   dt = t3/n;
   t = [dt/2 : dt : t3-dt/2];   % centers of bins
   tx1 = t - dt/2;              % left edge of bin
   tx2 = t + dt/2;              % right edge of bin
else
   t = x;
   n = length(t);
   dt = t(2)-t(1);
   tx1 = [];
   tx2 = [];
end

p = omegafun(t);
pcum = omegafuncum(t);

if length(x)==1
    % numerically integrate the PDF for each bin
    pint = zeros(n,1);
    xtol = 1e-12;
    for ii=1:n
        pint(ii) = integral(@(x) omegafun(x), tx1(ii), tx2(ii),'AbsTol',xtol,'RelTol',xtol);
    end
    
    % integration check
    disp(sprintf('integration check using %i bars: %.14f',n,sum(p)*dt));
    disp(sprintf('integration check using matlab: %.14f',sum(pint)));

    p = pint/dt;
end

% % cumulative density function
% cint = cumsum(pint);    % pint will be very accurate
% %cint = zeros(n,1);
% %for ii=1:n
% %    cint(ii) = integral(fun, t0, tx2(ii),'AbsTol',xtol,'RelTol',xtol);
% %end

if bfigure
   tsmooth = linspace(t0,t3,400);
   psmooth = omegafun(tsmooth);
   %pt1 = 4/pi;
   %pt2 = 4/pi*(-8/3 + sqrt(8));
   
   figure; nr=2; nc=1; ymax = 1.5;
   
   subplot(nr,nc,1); hold on;
   plot(tsmooth,psmooth,'b');
   plot(t,p,'r.--');
   %plot([t1 t1],[0 pt1],'k--',[t2 t2],[0 pt2],'k--');
   %plot([t0 t1 t2 t3],[0 pt1 pt2 0],'bo','markersize',14,'markerfacecolor','r','linewidth',1);
   axis equal; axis([0 t3 0 ymax]);
   xlabel('t, radians');
   ylabel('p_{\omega}(t)');
   title(sprintf('PDF for random moment tensors (n = %i)',n));
   
   subplot(nr,nc,2); hold on;
   bar(t,p,1,'c');
   plot(t,p,'r.');
   plot(tsmooth,psmooth,'b');
   %plot([t1 t1],[0 pt1],'k--',[t2 t2],[0 pt2],'k--');
   %plot([0 t1 t2 t3],[0 pt1 pt2 0],'bo','markersize',14,'markerfacecolor','r','linewidth',1);
   axis equal; axis([t0 t3 0 ymax]);
   xlabel('t, radians');
   ylabel('p_{\omega}(t)');
   title(sprintf('PDF for random moment tensors (n = %i)',n));
end

%--------------------------------------------------------------------------

function p = omegafun(t)

K = 3*pi/8;     % integral of sin^4(x)
p = sin(t).^4 / K;

%--------------------------------------------------------------------------

function pcum = omegafuncum(t)

% The cumulative distribution is the integral of the probability density sin(t).^4 / K
pcum = (12*t - 8*sin(2*t) + sin(4*t)) / (12*pi);

%==========================================================================
