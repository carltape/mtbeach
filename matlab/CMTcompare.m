function [omega,dM,omegadc,xi0,dgamma,ddelta,Delta,dthetaDC] = CMTcompare(M1,M2,bfigure)
%CMTCOMPARE compute metrics for two sets of moment tensors
%
% INPUT
%   M1,M2       6 x n moment tensors: M = [M11 M22 M33 M12 M13 M23]
%
% OUTPUT
%   omega       9-angle between moment tensors
%   dM          difference in magnitudes
%   omegadc     9-angle between closest double couples
%   xi0         minimum rotation angle between principal axes (kagan angle)
%   dgamma      difference in lune longitude
%   ddelta      difference in lune latitude
%   Delta       distance between lune points
%   dthetaDC    difference in thetaDC angles
%
% This is a wrapper function.
%
% dgamma, ddelta, and dthetaDC are signed, but there may be cases where
% using the absolute value would be desirable.
%
% calls CMT2omega.m, CMT2omegadc_xi0.m, CMT2TT.m
% 
% Carl Tape, 2026-09-17
%

if nargin==2, bfigure = false; end

omega = CMT2omega(M1,M2);

iorthoU = 0;
[omegadc,xi0] = CMT2omegadc_xi0(M1,M2,iorthoU);

[gamma1,delta1,M01,~,~,~,~,~,~,thetadc1,lam1] = CMT2TT(M1);
[gamma2,delta2,M02,~,~,~,~,~,~,thetadc2,lam2] = CMT2TT(M2);

% magnitude
Mw1 = m02mw(1,M01);
Mw2 = m02mw(2,M02);
dM = Mw2 - Mw1;

% arc distance between lune points
ulam1 = lam1 ./ vecnorm(lam1);
ulam2 = lam2 ./ vecnorm(lam2);
bdot = dot(ulam1,ulam2);
bdot(bdot > 1) = 1;
bdot(bdot <-1) = -1;
Delta = acosd(bdot);

% difference in lune longitude and lune latitude
dgamma = gamma2 - gamma1;
ddelta = delta2 - delta1;

% change in source type angle to DC
% (By the same reasoning, we could calculate the difference in angles to
% closest ISO or closest CLVD.)
dthetaDC = thetadc2 - thetadc1;

% make everything column vectors
omega   = omega(:);
dM      = dM(:);
omegadc = omegadc(:);
xi0     = xi0(:);
dgamma  = dgamma(:);
ddelta  = ddelta(:);
Delta   = Delta(:);
dthetaDC = dthetaDC(:);

if bfigure
    % compare orientation differences only
    figure; plot(omegadc,xi0,'ko','markersize',2,'markerfacecolor','k');
    xlabel('\omega_{DC}, degrees');
    ylabel('\xi_0, degrees');
    grid on;

    % Delta is <= omega (if the frames U1 = U2, then Delta = omega)
    figure; hold on;
    plot([0 180],[0 180],'r--')
    plot(Delta,omega,'ko','markersize',2,'markerfacecolor','k');
    xlabel('\Delta, degrees');
    ylabel('\omega, degrees');
    axis([0 180 0 180])
    grid on;

    % NOTE: Modifications will be needed to generalize the plotting in order
    %       to best display the values for a particular set of data.
    figure; nr=4; nc=2;
    ihist = 1; bmanual_ticks = true; dtick = 30;

    subplot(nr,nc,1); hold on; plot_histo(omega,[0:5:180],ihist);
    if bmanual_ticks, set(gca,'xtick',0:dtick:180); end
    xlabel('\omega : difference between moment tensors');
    title(sprintf('min = %.2f, max = %.2f',min(omega),max(omega)));

    dMX = round(max(abs(dM)));
    subplot(nr,nc,2); hold on; plot_histo(dM,[-dMX:0.1:dMX],ihist);
    xlabel('Mw2 - Mw1');
    title(sprintf('min = %.2f, max = %.2f',min(dM),max(dM)));

    subplot(nr,nc,3); hold on; plot_histo(omegadc,[0:5:180],ihist);
    if bmanual_ticks, set(gca,'xtick',0:dtick:180); end
    xlabel('\omega_{DC} : difference in orientation');
    title(sprintf('min = %.2f, max = %.2f',min(omega),max(omega)));

    subplot(nr,nc,4); hold on; plot_histo(xi0,[0:5:120],1);
    if bmanual_ticks, set(gca,'xtick',0:dtick:120); end
    xlabel('\xi_0 : difference in orientation');
    title(sprintf('min = %.2f, max = %.2f',min(xi0),max(xi0)));

    subplot(nr,nc,5); hold on; plot_histo(dgamma,[-60:5:60],ihist);
    if bmanual_ticks, set(gca,'xtick',-60:dtick:60); end
    xlabel('\gamma_2 - \gamma_1 : change in lune longitude');
    title(sprintf('min = %.2f, max = %.2f',min(dgamma),max(dgamma)));

    subplot(nr,nc,6); hold on; plot_histo(ddelta,[-180:10:180],1);
    if bmanual_ticks, set(gca,'xtick',-180:dtick:180); end
    xlabel('\delta_2 - \delta_1 : change in lune latitude');
    title(sprintf('min = %.2f, max = %.2f',min(ddelta),max(ddelta)));

    subplot(nr,nc,7); hold on; plot_histo(Delta,[0:10:180],ihist);
    if bmanual_ticks, set(gca,'xtick',0:dtick:180); end
    xlabel('\Delta : distance on the lune');
    title(sprintf('min = %.2f, max = %.2f',min(Delta),max(Delta)));

    subplot(nr,nc,8); hold on; plot_histo(dthetaDC,[-90:5:90],1);
    if bmanual_ticks, set(gca,'xtick',-90:dtick:90); end
    xlabel('\theta_{DC2} - \theta_{DC1} : change in DC angle');
    title(sprintf('min = %.2f, max = %.2f',min(dthetaDC),max(dthetaDC)));
end

%==========================================================================
% EXAMPLES

if 0==1
    % specify a set of moment tensors M1 and M2
    % here we choose a set of uniformly distributed moment tensors
    n = 10000;
    M1o = uniformMT(n);
    M2o = uniformMT(n);
    b_samemag = false;
    if b_samemag
        M1 = M1o; M2 = M2o; dMX = 1;
    else
        % generate magnitudes
        Mmin = 4; Mmax = 6;
        dMX = Mmax - Mmin;
        Mw1 = Mmin + dMX*rand(n,1);   % uniform distribution (rand)
        Mw2 = Mmin + dMX*rand(n,1);
        % assign magnitudes - 1
        M01o = CMT2m0(1,M1o);
        M01  = mw2m0(1,Mw1);
        Mrat1 = M01(:)' ./ M01o(:)';
        M1 = M1o .* Mrat1;
        % assign magnitudes - 2
        M02o = CMT2m0(1,M2o);
        M02  = mw2m0(1,Mw2);
        Mrat2 = M02(:)' ./ M02o(:)';
        M2 = M2o .* Mrat2;
    end

    bfigure = true;
    [omega,dM,omegadc,xi0,dgamma,ddelta,Delta,dthetaDC] = CMTcompare(M1,M2,bfigure);

end

%==========================================================================