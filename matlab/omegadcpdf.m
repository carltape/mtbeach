function [p,t] = omegadcpdf(x,bfigure)
%OMEGADCPDF probability density function for random orientations double couple moment tensor
%
% This function has been dubbed the 'mesa' based on its appearance.
%
% INPUT
%   x   A: an integer number of points for discretization OR
%       B: an input vector of points with for which to evaluate the PDF
%           INPUT IS IN RADIANS, NOT DEGREES
%
% OUTPUT
%   p   PDF discretized based on input choice
%   t   discretization of omega (IN RADIANS)
%
% See also omegapdf.m and kaganpdf.m
%
% EXAMPLES (see longer example below):
%   t = linspace(0,pi,181); p = omegadcpdf(t,true);
%   [p,t] = omegadcpdf(10,true);
%

if nargin==1, bfigure=false; end

% key points on curve
t0 = 0;
t1 = pi/3;
t2 = 2*pi/3;
t3 = pi;

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

% evaluate the PDF at the bin centers
p = omegafun(t);

% check area under curve
if 0==1
    xtol = 1e-12;
    deg = 180/pi;
    ta = t0; tb = t1; pint = integral(@(x) omegafun(x), ta, tb,'AbsTol',xtol,'RelTol',xtol);
    disp(sprintf('area from %6.1f to %6.1f = %14.12f',ta*deg,tb*deg,pint));
    ta = t1; tb = t2; pint = integral(@(x) omegafun(x), ta, tb,'AbsTol',xtol,'RelTol',xtol);
    disp(sprintf('area from %6.1f to %6.1f = %14.12f',ta*deg,tb*deg,pint));
    ta = t2; tb = t3; pint = integral(@(x) omegafun(x), ta, tb,'AbsTol',xtol,'RelTol',xtol);
    disp(sprintf('area from %6.1f to %6.1f = %14.12f',ta*deg,tb*deg,pint));
    ta = t0; tb = t3; pint = integral(@(x) omegafun(x), ta, tb,'AbsTol',xtol,'RelTol',xtol);
    disp(sprintf('area from %6.1f to %6.1f = %14.12f',ta*deg,tb*deg,pint));
end

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
   
   figure; nr=2; nc=1; 
   
   subplot(nr,nc,1); hold on;
   plot(tsmooth,psmooth,'b');
   plot(t,p,'r.--');
   %plot([t1 t1],[0 pt1],'k--',[t2 t2],[0 pt2],'k--');
   %plot([t0 t1 t2 t3],[0 pt1 pt2 0],'bo','markersize',14,'markerfacecolor','r','linewidth',1);
   axis equal; axis([0 t3 0 1.5]);
   xlabel('t, radians');
   ylabel('p_{\omega}(t)');
   title(sprintf('PDF for random orientations (n = %i)',n));
   
   subplot(nr,nc,2); hold on;
   % bars will be centered on the input points (which you may have intended
   % to be edges); therefore some bins may go outside the x-plotting range
   bar(t,p,1,'c');
   plot(t,p,'r.');
   plot(tsmooth,psmooth,'b');
   %plot([t1 t1],[0 pt1],'k--',[t2 t2],[0 pt2],'k--');
   %plot([0 t1 t2 t3],[0 pt1 pt2 0],'bo','markersize',14,'markerfacecolor','r','linewidth',1);
   axis equal; axis([t0 t3 0 1.5]);
   xlabel('t, radians');
   ylabel('p_{\omega}(t)');
   title(sprintf('PDF for random orientations (n = %i)',n));
end

%--------------------------------------------------------------------------

function p = omegafun(t)
% the mesa function -- note that this is APPROXIMATED by a pre-evaluated
% function; in the future we could port the mathematica equations here

t0 = 0;
t1 = pi/3;
t2 = 2*pi/3;
t3 = pi;

i1 = find(and(t >= t0, t <= t1));
i2 = find(and(t > t1, t <= t2));
i3 = find(and(t > t2, t <= t3));

% % see estimation below (based on 10 million values)
% P2poly  = [-0.2752 0 0.6960];
% P13poly = [0.4855 0 -0.6828 0 0.3798 0 0.2766 0 0];
% 
% pshift1 = 0;
% pshift2 = pi/2;
% pshift3 = pi;
% 
% p(i1) = polyval(P13poly,t(i1) - pshift1);
% p(i2) = polyval( P2poly,t(i2) - pshift2);
% p(i3) = polyval(P13poly,t(i3) - pshift3);

% from Mathematica calculations
%P1omega = [0 5 10 15 20 25 30 35 40 45 50 55 57 59 60]'*pi/180;
%P1pts = [0 0.0024 0.0097 0.0221 0.0396 0.0627 0.0919 0.1280 0.1722 0.2267 ...
%    0.2959 0.3906 0.4431 0.5192 0.6191]';
P1omega = [0:0.5:60]*pi/180;
P1pts = [
     0.
0.000024241
0.0000969677
0.000218191
0.00038793
0.00060621
0.000873064
0.00118853
0.00155267
0.00196552
0.00242716
0.00293765
0.00349707
0.00410552
0.00476308
0.00546987
0.00622598
0.00703155
0.00788669
0.00879156
0.00974629
0.010751
0.011806
0.0129113
0.0140671
0.0152736
0.0165312
0.0178398
0.0191998
0.0206115
0.022075
0.0235906
0.0251587
0.0267794
0.0284531
0.03018
0.0319606
0.0337951
0.035684
0.0376275
0.039626
0.04168
0.0437899
0.0459561
0.0481791
0.0504593
0.0527972
0.0551934
0.0576484
0.0601627
0.0627369
0.0653717
0.0680677
0.0708256
0.073646
0.0765297
0.0794773
0.0824898
0.0855679
0.0887126
0.0919247
0.0952051
0.0985549
0.101975
0.105467
0.109031
0.112669
0.116382
0.120172
0.124039
0.127986
0.132013
0.136123
0.140316
0.144596
0.148964
0.153421
0.15797
0.162614
0.167354
0.172193
0.177134
0.18218
0.187334
0.192599
0.197979
0.203476
0.209097
0.214844
0.220723
0.226738
0.232896
0.239201
0.245661
0.252283
0.259075
0.266044
0.273201
0.280556
0.28812
0.295907
0.303932
0.31221
0.320761
0.329606
0.33877
0.348282
0.358177
0.368495
0.379286
0.39061
0.402544
0.415183
0.428658
0.443143
0.458889
0.476281
0.495959
0.519169
0.549105
0.619061
    ]';
itype = 'linear';
%itype = 'pchip';
p(i1) = interp1(P1omega,P1pts,t(i1),itype);
p(i3) = interp1(-P1omega+pi,P1pts,t(i3),itype);

% middle parabola (note: we assume that the curve is a parabola,
% but it may not even be a polynomial for that matter)
%P2poly  = [-0.2827601209636421 0.6965810255067267];   % Walt Mathematica
P2poly  = [-0.28276 0 0.696581];
%P2poly  = [-0.282112 0 0.696403];
pshift2 = pi/2;
p(i2) = polyval( P2poly,t(i2) - pshift2);

%==========================================================================
% EXAMPLE

if 0==1
    % plot the probability density function
    %% option 1: provide input points in vector t
    n = 181; t = linspace(0,pi,n); p = omegadcpdf(t,true);
    %% option 2: specify the number of discretization points
    [p,t] = omegadcpdf(10,true);
    %t = linspace(0,pi,11); p = omegadcpdf(t,true);
    
    %% cumulative density function
    dt = 3 * pi/180;        % try 0.5, 1, 3, 10
    t = [0:dt:pi];
    p = omegadcpdf(t,true);
    pcum = cumsum(p)*dt;    % CRUDE integration
    pcum(end)               % check (ideally should be 1)
                            % 1.00024 (0.5), 1.00069 (1), 1.0035 (3), 1.0216 (10)
    figure; plot(t,pcum,'.'); axis equal; axis([0 pi 0 1.01]);
    
    %% compare omegaDC PDF with omega PDF
    nz = 1000;
    t = linspace(0,pi,nz);
    f = omegapdf(t);
    fdc = omegadcpdf(t);
    % note: the sum between two random variables is convolution;
    %       the difference between two random variables is correlation
    %figure; plot(xcorr(f,fdc));
    %figure; plot(conv(f,fdc));
    fndc = (f - fdc) / sum(f - fdc);
    figure; hold on;
    plot(t,f,'g','linewidth',2);
    plot(t,fdc,'r','linewidth',2);
    plot(t,fndc,'k','linewidth',2);
    legend('PDF \omega','PDF \omegaDC','\omega - \omegaDC');
    plot([0 pi],[0 0],'k--');
    axis([0 pi -1 1.5]);
    xlabel('\omega, radians');
    
    %% plot a pre-saved distribution
    n = 1e4;
    %n = 1e7;
    load(sprintf('/home/carltape/dwrite/mt/omega_xi0_%9.9i',n));
    plot_omegadc(omega);
    plot_xi0(xi0);
end

%==========================================================================
% ESTIMATING THE PDF FROM A SET OF omegas
% EVENTUALLY THIS CAN BE DELETED

if 0==1
    %%
    clear, close all, clc
    bsaveomega = false;

    if bsaveomega
        deg = 180/pi;
        n = 1e7;
        % generate random orientations
        kappa = randomvec(0,360,n);
        sigma = randomvec(-90,90,n);
        h = randomvec(0,1,n);
        theta = deg * acos(h);
        % calculate DC with random oreientations
        gamma = zeros(n,1);
        delta = zeros(n,1);
        M0 = sqrt(2)*ones(n,1);
        M1 = TT2CMT(gamma,delta,M0,kappa,theta,sigma);
        figure; plotMT_TT(gamma,delta,M0,kappa,theta,sigma);

        % here we just fix the reference moment tensor -- in principle we
        % could take any moment tensor from within the set and then
        % calculate the omega to all other moment tensors
        Mref = 1/sqrt(2)*[1 0 -1 0 0 0]';
        M2   = repmat(Mref,1,n);
        iortho = 0;
        % note: probably safter to use CMT2omega.m
        [omega,xi0,U] = CMT2omegadc_xi0(M1,M2,iortho);
        save(sprintf('/home/carltape/dwrite/mt/omega_xi0_%9.9i',n),'omega','xi0');
        return
    else
       n = 1e7;
       load(sprintf('/home/carltape/dwrite/mt/omega_xi0_%9.9i',n));
    end

    %plot_xi0(xi0);
    %%
    nbin = 180;
    [N,Nplot,centers] = plot_omegadc(omega,nbin);
    return

    figure; nr=2; nc=1;
    subplot(nr,nc,1); hold on;
    plot(centers,Nplot,'k.-');
    %axis equal, axis tight

    % three sections of the mesa
    % include one extra point from inside the middle region
    x = centers;
    i1 = find(x <= pi/3);
    i3 = find(x >= 2*pi/3);
    %i1 = [i1 ; i1(end)+1];
    %i3 = [i3(1)-1 ; i3];
    i13 = [i1 ; i3];
    %i13 = find(or(x <= pi/3, x >= 2*pi/3));
    i2 = find(and(x >= pi/3, x <= 2*pi/3));

    %Norder = 4;
    %P1 = polyfit(centers(i1),Nplot(i1),Norder);
    x = centers(i1); d = Nplot(i1);
    P = [x.^8 x.^6 x.^4 x.^2]\d;
    P1 = [P(1) 0 P(2) 0 P(3) 0 P(4) 0 0]
    %P = [x.^6 x.^4 x.^2]\d;
    %P1 = [P(1) 0 P(2) 0 P(3) 0 0]
    Nfit1 = polyval(P1,centers(i1));
    plot(centers(i1),Nfit1,'r--');

    pshift2 = pi/2;
    %P2 = polyfit(centers(i2)-pshift,Nplot(i2),Norder);
    x = centers(i2)-pshift2; d = Nplot(i2);
    %P = [x.^4 x.^2 x.^0]\d;
    %P2 = [P(1) 0 P(2) 0 P(3)]
    P = [x.^2 x.^0]\d;
    P2 = [P(1) 0 P(2)]
    Nfit2 = polyval(P2,centers(i2)-pshift2);
    plot(centers(i2),Nfit2,'r--');

    pshift3 = pi;
    %P3 = polyfit(centers(i3)-pshift,Nplot(i3),Norder);
    x = centers(i3)-pshift3; d = Nplot(i3);
    P = [x.^8 x.^6 x.^4 x.^2]\d;
    P3 = [P(1) 0 P(2) 0 P(3) 0 P(4) 0 0]
    Nfit3 = polyval(P3,centers(i3)-pshift3);
    plot(centers(i3),Nfit3,'r--');

    P13poly = mean([P1 ; P3])
    
    subplot(nr,nc,2); hold on;
    plot(centers(i1),Nfit1-Nplot(i1),'k.-');
    plot(centers(i2),Nfit2-Nplot(i2),'k.-');
    plot(centers(i3),Nfit3-Nplot(i3),'k.-');

    % SUMMARY HERE:
    % 1. The sin^n(x) function does not fit at all
    % 2. The exponential function is pretty good, but not as good as the
    %    polynomial fit.
    for itype = 1:2
    nx = 100;
    switch itype
        case 1
            % sine fit
            pvec = linspace(1,6,nx);
            avec = exp(linspace(-1.5,1,nx));
        case 2
            % exponential
            pvec = exp(linspace(-1.5,-0.5,nx));
            avec = exp(linspace(-5,-2,nx));
    end

    [P,A] = meshgrid(pvec,avec);
    p = P(:);
    a = A(:);
    n = length(a);
    fall = NaN(n,1);
    x = centers;
    for ii=1:n
        switch itype
            case 1
                yfit = a(ii)*sin(x).^p(ii);
                fall(ii) = norm( yfit(i13) - Nplot(i13) );
            case 2
                yfit = a(ii)*(exp(x/p(ii)) - 1);
                fall(ii) = norm( yfit(i1) - Nplot(i1) );
        end
    end
    figure;
    if itype==1, scatter(p,log(a),5^2,fall,'fill'); end
    if itype==2, scatter(log(p),log(a),5^2,fall,'fill'); end
    xlabel('p'); ylabel('a');
    caxis([0 1]);
    colorbar

    [~,imin] = min(fall);
    figure; hold on; plot(x,Nplot,'k.-');
    if itype==1, plot(x,a(imin)*sin(x).^p(imin),'r--'); end
    if itype==2, plot(x,a(imin)*(exp(x/p(imin))-1),'r--'); end
    ylim([0 1]);
    end

end
    
%==========================================================================
