%
% run_uniformMT.m
%
% run script for uniformMT.m
% Examples to generate subsets of uniform moment tensors
% See TapeTape2015 "A uniform parameterization of moment tensors"
%
% Carl Tape, 2015-07-26
%

clear, clc, close all

deg = 180/pi;
bprint = false;

% GET USER INPUT
stlabs = {  'random full moment tensor',
            'random moment tensors for fixed nu',
            'random deviatoric moment tensors',
            'random moment tensors with fixed eigenvalues',
            'random double couple moment tensors',
            'regular grid of full moment tensors',
            'regular grid of double couple moment tensors',
            'will generate error (intentionally)',
            'insights into the sin^4(omega) distribution',
            'insights into random vs uniform grid'};
nex = length(stlabs);
disp('run_uniformMT.m examples:');
for ii=1:nex, disp(sprintf('   %2i : %s',ii,stlabs{ii})); end
iex = input(sprintf('Select your example (1-%i), then hit ENTER: ',nex));  
if ~any(iex==[1:nex])
    error('iex must be 1-%i',nex);
else
    stlab1x = stlabs{iex};
    disp(sprintf('run_uniformMT.m example %i: %s',iex,stlab1x));
end

% reference moment tensor for setting omega=0 in the distributions
% (alternatively you can choose a moment tensor that is in the set)
Mref = 1/sqrt(2)*[1 0 -1 0 0 0]';

% default number of moment tensors to generate
n = 1e5;

switch iex
    case 1
    % randomly generated uniform full moment tensors
    [M,v,w,kappa,sigma,h] = uniformMT(n);
    %iref = randi(n); Mref = M(:,iref);
    omega = CMT2omega(Mref,M);
    plot_omega(omega);  

    case 2
    % randomly generated uniform moment tensors for fixed nu (or phi)
    nu0 = 0.25;
    %nu0 = -1;   % deviatoric moment tensors
    blune = false;
    brandom = true;
    [gamma,delta] = nu2lune(nu0,n,blune,brandom);
    % random uniform orientations (see uniformMT.m)
    kappa = randomvec(0,360,n);
    sigma = randomvec(-90,90,n);
    h     = randomvec(0,1,n);
    theta = acos(h)*deg;
    M0 = 1/sqrt(2);     % ensure that ||M|| = rho = 1
    % moment tensors
    M = TT2CMT(gamma,delta,M0,kappa,theta,sigma);
    [v,w] = lune2rect(gamma,delta);
    
    case 3
    % randomly generated uniform deviatoric moment tensors
    % This is equivalent to case 2 with nu0 = -1.
    % WARNING: One cannot simply draw random points on a line segment
    %          between two points (v1,w1) and (v2,w2) in the vw plane, and
    %          then use random orientations. But this does work for the case
    %          of w = constant [line of latitude].
    w0 = 0;  % could be any permissible w value +/-3pi/8 = +/- 1.1781
             % here we choose it for the deviatoric MTs
    [~,v,w,kappa,sigma,h] = uniformMT(n);
    w = w0*ones(size(w));    % reassign all w values to w0
    rho = sqrt(2)*ones(n,1);
    M = TT152CMT(rho,v,w,kappa,sigma,h);

    case 4
    % randomly generated uniform moment tensors with fixed lune point
    gamma0 = -25; delta0 = 60;     % double couple
    %gamma0 = -75; delta0 = 60;     % will generate error exit
    [M,v,w,kappa,sigma,h] = uniformMT(n,gamma0,delta0);
    %iref = randi(n); Mref = M(:,iref);
    omega = CMT2omega(Mref,M);

    case 5
    % randomly generated uniform double couple moment tensors
    gamma0 = 0; delta0 = 0;     % double couple
    [M,v,w,kappa,sigma,h] = uniformMT(n,gamma0,delta0);
    %iref = randi(n); Mref = M(:,iref);
    omega = CMT2omega(Mref,M);
    plot_omegadc(omega);  

    case 6
    % regular grid of uniform full moment tensors
    %n = [2 3 5 4 2];
    %n = [6 18 18 9 5];    % 5-degree increments in strike/dip/rake
    %n = [6 18 36 18 10];  % 10-degree increments in strike/dip/rake
    %n = [7 7 7 7 7];
    n = [13 35 11 11 11];   % lune points match Alvizuri et al. (2018), Figure 5d
                            % IF you set UNIFORM_GAMMAB = true in uniformMT.m
    [M,v,w,kappa,sigma,h] = uniformMT(n);
    ntotal = length(v);

    case 7
    % regular grid of uniform double couple moment tensors
    ntotal = 1e5;
    %ntotal = 46656;                 % number of points for n = [0 0 72 36 18]
    npart = round(ntotal^(1/3));
    %n = [0 0 5 4 2];
    %n = [0 0 36 18 9];              % 10-degree increments
    n = [0 0 72 36 18];             % 5-degree increments
    %n = [0 0 npart npart npart];    % equal number of points for each interval
    %n = [0 0 8 180 100];            % PROBLEMATIC!
    gamma0 = 0; delta0 = 0;     % double couple
    [M,v,w,kappa,sigma,h] = uniformMT(n,gamma0,delta0);
    ntotal = length(v);
    %iref = randi(ntotal); Mref = M(:,iref);
    omega = CMT2omega(Mref,M);
    plot_omegadc(omega); 

    case 8
    %--> THIS WILL GENERATE AN ERROR, since the 0 0 indicates a fixed lune point,
    % but no lune point is specified as the 2nd and 3rd arguments of uniformMT
    n = [0 0 5 4 2];    % error
    %n = [1 1 5 4 2];    % error
    %n = [2 2 5 4 2];    % no error
    [M,v,w,kappa,sigma,h] = uniformMT(n);

    case 9
    % insights into the sin^4(omega) distribution
    figure; nr=2; nc=1;

    % angular distance from a point on a circle (1-sphere) to all other points on a circle
    n = 1e6;
    theta = 2*pi*rand(n,1);
    [rx,ry] = pol2cart(theta,1);
    ipick = randi(n);
    rx0 = rx(ipick)*ones(n,1);
    ry0 = ry(ipick)*ones(n,1);
    theta = acos( rx0.*rx + ry0.*ry );
    subplot(nr,nc,1); plot_histo(theta*deg,[0:2:180],3); ylim([0 1e-2]);
    title('angular distance from any point on the circle to uniformly distributed points');
    %figure; hold on; plot(rx,ry,'.'); plot(rx(ipick),ry(ipick),'ro');

    % angular distance from a point on a sphere (2-sphere) to all other points
    n = 1e6;
    z = randomvec(-1,1,n);
    ylat = asin(z)*deg;         % delta = 90 - acos(z)*deg;
    xlon = randomvec(0,360,n);
    ipick = randi(n);
    xlon0 = xlon(ipick)*ones(n,1);
    ylat0 = ylat(ipick)*ones(n,1);
    % this would be cleaner if written as a dot product
    theta = distance(ylat0,xlon0,ylat,xlon,'degrees');

    subplot(nr,nc,2); hold on;
    [N,Nplot,centers] = plot_histo(theta,[0:2:180],3); ylim([0 1e-2]);
    x = centers;
    f = sin(x/deg);
    plot(x,f/2/deg,'r.-');
    title('distance from any point on the sphere to uniformly distributed points');
    xlabel('angular distance, deg');
    %figure; hold on; plot(xlon,ylat,'.'); plot(xlon(ipick),ylat(ipick),'ro');
    
    error
    
    case 10
    % insights into random vs regular grid
    
    nb = 16;
    a1=0; a2=2; b1=0; b2=1;
    aran = a2-a1; bran = b2-b1;
    n = nb*(2*nb);
    ar = 2*rand(n,1); br = rand(n,1);
    da = sqrt((aran*bran)/n);
    avec = [(a1+da/2) : da : (2-da/2)]';
    db = da; bvec = [(b1+db/2) : da : (1-db/2)]';
    [A,B] = meshgrid(avec,bvec); ag = A(:); bg = B(:);
    
    figure; nr=2; nc=1;
    subplot(nr,nc,1); plot(ar,br,'.'); axis equal, axis([a1 a2 b1 b2]);
    title(sprintf('%i uniform random points',length(ar)));
    subplot(nr,nc,2); plot(ag,bg,'.'); axis equal, axis([a1 a2 b1 b2]);
    title(sprintf('%i uniform grid points',length(ag)));
    
    error
        
end

%==========================================================================

% Question: How many fewer points are needed to calculate a uniform
% distribution of moment tensors, in comparison to using assuming the
% eigenvalue triples are uniformly distributed on the lune?
%
% Answer (see TapeTape2015):
% Suppose person H takes eigenvalue triples uniformly distributed on the lune and
% then takes the same number of random orientations at each of his lune points.
% And suppose H wants to achieve the same density of moment tensors that we
% have at the double couple. Then from our Eq 44b H wants density 4/pi
% everywhere on the lune. The integral of that over the lune is
% 4/pi x lune area = 4/pi x (4 pi^2)/6 = 8 pi/3.
% The integral of our density is 1, so H is using about 8 times as many points.
%
%disp(sprintf('savings in number of gridpoints is the factor %.2f',8*pi/3));

%==========================================================================
% plot various moment tensor quantities

disp('plotting various quantities...');

n = length(kappa);
M0 = CMT2m0(1,M);
theta = acos(h)*deg;
[gamma,delta] = rect2lune(v,w);

% if 0
%     % check sigmaproj vs sigma
%     bdisplay = false;
%     [nu,alpha,N1,N2,lam] = CMT2faultpar(M,bdisplay);
%     N = N1; S = N2;
%     bSproj = true;
%     [kappax,thetax,sigmax,K,N,S,sigmaprojx] = NS2sdr(N,S,bdisplay,bSproj);
%     
%     figure; nr=3; nc=1;
%     subplot(nr,nc,1); plot_histo(sigma,[-90:4:90]);
%     subplot(nr,nc,2); plot_histo(sigmax,[-90:4:90]);
%     subplot(nr,nc,3); plot_histo(sigmaprojx,[-90:4:90]);
%     
%     figure; plot(sigmax,sigmaprojx,'.');
%     grid on; axis([-90 90 -90 90]);
%     
%     figure; plot(sigmax,sigma,'.');
%     grid on; axis([-90 90 -90 90]);
% 
%     % fault angles for the two DCs that make up the two CDCs
%     % (this is slow because CMT2CDC.m is not vectorized)
%     kappa1 = NaN(n,1); theta1 = NaN(n,1); sigma1 = NaN(n,1);
%     kappa2 = NaN(n,1); theta2 = NaN(n,1); sigma2 = NaN(n,1);
%     for ii=1:n
%         if mod(ii,1000)==0, disp(sprintf('==== %i out of %i ====',ii,n)); end
%         Mtar = Mvec2Mmat(M(:,ii),1);
%         [K1,D1,K2,D2] = CMT2CDC(Mtar);
%         [~,~,~,kappa1(ii),theta1(ii),sigma1(ii)] = CMT2TT(Mvec2Mmat(D1,0),bdisplay);
%         [~,~,~,kappa2(ii),theta2(ii),sigma2(ii)] = CMT2TT(Mvec2Mmat(D2,0),bdisplay);
%     end
% end

% lune longitude, lune latitude, strike, dip, slip
plotMT_TT(gamma,delta,M0,kappa,theta,sigma);
if bprint, orient tall; print(gcf,'-dpng',sprintf('run_uniformMT_iex%i_hist1',iex)); end
% v, w, strike, dip, slip
plotMT_TT(v,w,M0,kappa,theta,sigma,true);
if bprint, orient tall; print(gcf,'-dpng',sprintf('run_uniformMT_iex%i_hist2',iex)); end

% plot correlation matrices
figure; nr=2; nc=1; clims = [-1 1];
% moment tensor entries
R = corrcoef(M');
Mlab = {'Mrr','Mtt','Mpp','Mrt','Mrp','Mtp'};
subplot(nr,nc,1); imagesc(R); axis equal; axis tight; caxis(clims); colorbar
set(gca,'xtick',[1:6],'xticklabel',Mlab,'ytick',[1:6],'yticklabel',Mlab);
title(sprintf('%s (n = %i)',stlab1x,n));
% TT2015 coordinates
MTT = [v w kappa sigma h]'; Mlab = {'v','w','kappa','sigma','h'};
%MTT = [v w kappa sigma theta]'; Mlab = {'v','w','kappa','sigma','theta'};
RTT = corrcoef(MTT');
subplot(nr,nc,2); imagesc(RTT); axis equal; axis tight; caxis(clims); colorbar
set(gca,'xtick',[1:6],'xticklabel',Mlab,'ytick',[1:6],'yticklabel',Mlab);
if bprint, orient tall; print(gcf,'-dpng',sprintf('run_uniformMT_iex%i_corr',iex)); end

% source types
figure;
subplot(1,2,1); plot(gamma,delta,'.'); axis equal, axis([-30 30 -90 90])
xlabel('lune longitude'); ylabel('lune latitude');
subplot(1,2,2); plot(v,w,'.'); xlabel('v'); ylabel('w');
axis equal, axis([-1/3 1/3 -3*pi/8 3*pi/8])

error

%% ADDITIONAL HISTOGRAMS

rho = sqrt(2)*M0;
mfac = mean(rho);

% Mij entries
%Medges = [-sqrt(2):0.1:sqrt(2)];
Medges = linspace(-mfac,mfac,21);
plotMT_Mij(M,Medges);

% eigenvalues
[lam,U] = CMTdecom(M);
plotMT_lam(lam);

error

%% ADDITIONAL HISTOGRAMS (note: lam must be evaluated in previous block)

% lune longitude and latitude within latitude bands and longitude bands
plotMT_lune(gamma,delta);

% classical model and CDC
[nu,alpha] = lam2nualpha(lam);
[phi,zeta] = lam2phizeta(lam);
plotMT_phizeta(phi,zeta);
plotMT_nualpha(nu,alpha);
if iex==2
    [alpha,falpha] = pdf_alpha(nu0);
    plot(alpha,falpha,'r--','linewidth',2);
end

%% ADDITIONAL HISTOGRAMS (note: U must be evaluated in prior block)

% plunge-azimuth angles of eigenvectors
Useu = convertv(1,5,U);
Uout = U2pa(Useu,1);
plotMT_eigvec(Uout(:,1),Uout(:,2),Uout(:,3),Uout(:,4),Uout(:,5),Uout(:,6));

%==========================================================================
