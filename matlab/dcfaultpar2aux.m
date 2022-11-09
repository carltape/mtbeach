function [kap2,theta2,sigma2] = dcfaultpar2aux(kap0,theta0,sigma0,bdisplay)
%DCFAULTPAR2AUX given strike-dip-rake, compute strike-dip-rake for auxiliary plane
% 
% INPUT (and OUTPUT):
%   kap0    strike (0 to 360)
%   theta0  dip (0 to 90)
%   sigma0  rake (-180 to 180)
%
% Potentially we could instead force a convention of this for sigma
% INPUT:    |rake| <= 90 (-90 to 90)
% OUTPUT:   |rake| > 90 (180 to -90 or 90 to 180)
%
% Very likely the codes could be greatly simplified; for now there are two
% approaches that seem to give identical (and correct) results.
%
% See also run_CMT2sph.m
% See Tape and Tape (2012 GJI) for details (especially Eq 37).
%
% Carl Tape, 2011-08-01
%

% column vectors
kap0 = kap0(:);
theta0 = theta0(:);
sigma0 = sigma0(:);
n = length(kap0);

if nargin==3
    bdisplay = false;
end

deg = 180/pi;

% iversion=1 (DEFAULT) is much faster and more reliable than option 2
% iversion=2 is essentially a check on other algorithms
iversion = 1;

if iversion==1
    costheta = cos(theta0/deg);
    sintheta = sin(theta0/deg);
    cossigma = cos(sigma0/deg);
    sinsigma = sin(sigma0/deg);
    coskappa = cos(kap0/deg);
    sinkappa = sin(kap0/deg);
    
    % DIP: since theta0 is [0,90], so will this theta
    theta2 = deg * acos( abs(sinsigma) .* sintheta);
    
    % RAKE: acos will return angle from 0 to 180, but we want -180 to 180
    sigma2 = deg * acos(-cossigma.*sintheta ./ sqrt(cossigma.^2 + costheta.^2.*sinsigma.^2));
    bneg = sigma0 < 0;
    sigma2(bneg) = sign(sigma0(bneg)) .* sigma2(bneg);
    
    % STRIKE
    Kxyz = zeros(n,3);
    Kxyz(:,1) = -(-sinkappa.*cossigma + coskappa.*sinsigma.*costheta);
    Kxyz(:,2) = -(coskappa.*cossigma + sinkappa.*sinsigma.*costheta);
    ineg = find(sign(sigma0) < 0);
    Kxyz(ineg,:) = -Kxyz(ineg,:);
    [~,kap2] = xyz2tp(Kxyz');
    kap2 = wrapTo360(kap2*deg);
    
elseif iversion==2
    % obtain both fault planes by converting to MT and back
    MDC = dcfaultpar2CMT(kap0,theta0,sigma0);
    %Mx = TT2CMT(0,0,1,kap0,theta0,sigma0);
    [~,kap1,theta1,sigma1,kap2,theta2,sigma2] = CMT2dcfaultpar(MDC,0);

    % WARNING: this is an ad hoc way to determine whether the input plane
    % matches the first or second fault plane returned from CMT2dcfaultpar.m
    DTRSH = 0.1;
    bkap = (abs(kap1 - kap0) > DTRSH);
    btheta = (abs(theta1 - theta0) > DTRSH);
    bsigma = (abs(sigma1 - sigma0) > DTRSH);
    ball = bkap + btheta + bsigma;
    ikeep = find(ball==0);      % default: don't swap
    iswap = find(ball > 0);
    if n ~= length(ikeep) + length(iswap)
        ibad = setdiff([1:n]', [ikeep ; iswap]);
        bkap(ibad), btheta(ibad); bsigma(ibad)
        for jj=1:length(ibad)
           ii = ibad(jj);
           disp(sprintf(['%4i input ' stfmt ', 1 ' stfmt ', 2 ' stfmt],...
               ii,kap0(ii),theta0(ii),sigma0(ii),kap1(ii),theta1(ii),sigma1(ii),kap2(ii),theta2(ii),sigma2(ii)));
        end
        error('non-recovery of the input fault parameters');
    end

    kap2(iswap) = kap1(iswap);
    theta2(iswap) = theta1(iswap);
    sigma2(iswap) = sigma1(iswap);
end

% do not allow -180 rake angles (see readCMT.m)
sigma2(sigma2==-180) = 180;

% display info
stfmt = '(%7.3f, %6.3f, %8.3f)';
if bdisplay
    if iversion==1
        for ii=1:n
           disp(sprintf(['%4i input ' stfmt ' --> ' stfmt],...
               ii,kap0(ii),theta0(ii),sigma0(ii),kap2(ii),theta2(ii),sigma2(ii)));
        end
    elseif iversion==2
        for ii=1:n
            disp(sprintf(['%4i input ' stfmt ', 1 ' stfmt ', 2 ' stfmt],...
            ii,kap0(ii),theta0(ii),sigma0(ii),kap1(ii),theta1(ii),sigma1(ii),kap2(ii),theta2(ii),sigma2(ii)));
        end
    end
end

%==========================================================================
% EXAMPLES

if 0==1
    clc, close all, clear
    
    %% EXAMPLE: given one nodel plane and slip direction,
    %          get the other nodal plane and slip direction
    kap0 = 200; theta0 = 52; sigma0 = -85;
    [kap2,theta2,sigma2] = dcfaultpar2aux(kap0,theta0,sigma0);
    [kap3,theta3,sigma3] = dcfaultpar2aux(kap2,theta2,sigma2);
    disp('input strike, dip, rake:'); kap0, theta0, sigma0
    disp('other strike, dip, rake:'); kap2,theta2,sigma2
    disp('check by converting to original:'); kap3,theta3,sigma3
    disp('check that both triples correspond to the same MT:');
    MDC1 = dcfaultpar2CMT(kap0,theta0,sigma0);
    MDC2 = dcfaultpar2CMT(kap2,theta2,sigma2);
    disp([MDC1 MDC2]);
    disp('only one triple is within the MT orientation space');
    disp('strike = [0,360], dip = [0,90], rake = [-90,90]');
    
    %% EXAMPLE 1: two pairs of nodal planes -- 
    %            perform operation twice to recover the starting angles
    kap0 = [126 349.8 87 290]'; theta0 = [43 56 10 81]'; sigma0 = [-125 -61.9 67 93]';
    [kap2,theta2,sigma2] = dcfaultpar2aux(kap0,theta0,sigma0);
    [kap3,theta3,sigma3] = dcfaultpar2aux(kap2,theta2,sigma2);
    [kap0 kap3 theta0 theta3 sigma0 sigma3]
    [kap0-kap3 theta0-theta3 sigma0-sigma3]
    
    %% EXAMPLE 2: random points in strike-dip-rake space
    n = 1e5;
    kap0 = randomvec(0,360,n);
    h = randomvec(0,1,n);
    theta0 = acos(h)*180/pi;
    sigma0 = randomvec(-180,180,n);
    %sigma0 = randomvec(-90,90,n);
    [kap2,theta2,sigma2] = dcfaultpar2aux(kap0,theta0,sigma0);
    
    % plot input plane and 2nd plane
    plotMT_dc(kap0,theta0,sigma0,kap2,theta2,sigma2);

    % recover the original angles
    [kap3,theta3,sigma3] = dcfaultpar2aux(kap2,theta2,sigma2);
    figure; nr=3; nc=1;
    subplot(nr,nc,1); plot(kap3-kap0,'.'); xlabel('strike angle, degrees');
    subplot(nr,nc,2); plot(theta3-theta0,'.'); xlabel('dip angle, degrees');
    subplot(nr,nc,3); plot(sigma3-sigma0,'.'); xlabel('rake angle, degrees');
    
    %% fault angles listed in the GCMT catalog
    % (here we are ignoring any non-double-couple components of the moment tensors)
    clc, close all, clear
    isub = [1:59681]';
    %isub = [1:4000]';
    %isub = [117 234 366 472 526]';
    %isub = [53 2068 2494]';  % for itype=2 dip discrepancies
    [otime,~,~,lat,lon,dep,M,~,~,eid,~,str1,dip1,rk1,str2,dip2,rk2] = readCMT(isub);
    %[otime,~,~,lat,lon,dep,M,~,~,eid,~,str1,dip1,rk1,str2,dip2,rk2] = readCMT;
    n = length(otime);

    % the GCMT convention is that the first set of fault angles has the shallower dip angle
    figure; plot(dip2-dip1,'.');
    
    % if the first set of fault angles has |rake| < 90, then the second set
    % of fault angles will have |rake| > 90 (and vice versa)
    iUbox = find(abs(rk1) <= 90);
    iaux = setdiff(isub,iUbox);
    figure; nr=2; nc=2; edges = [-180:5:180];
    subplot(nr,nc,1); plot_histo(rk1(iUbox),edges); title('1st rake angle listed');
    subplot(nr,nc,2); plot_histo(rk1(iaux),edges); title('1st rake angle listed')
    subplot(nr,nc,3); plot_histo(rk2(iaux),edges); title('2nd rake angle listed');
    subplot(nr,nc,4); plot_histo(rk2(iUbox),edges); title('2nd rake angle listed');
    %fontsize(16); print(gcf,'-dpng','GCMT_slip_angle');
    
    % the 2nd plane is the one with a larger dip angle
    figure; plot(dip1,dip2,'.'); axis equal, axis([0 90 0 90])
    xlabel('1st dip angle listed');
    ylabel('2nd dip angle listed');
    %fontsize(16); print(gcf,'-dpng','GCMT_dip_angle');
    
    % key command: get the angles of the auxilliary plane
    [str2x,dip2x,rk2x] = dcfaultpar2aux(str1,dip1,rk1);
    [str1x,dip1x,rk1x] = dcfaultpar2aux(str2,dip2,rk2);
    
    % check
    ATRSH = 3;  % need large threshold, since GCMT angles ar rounded to nearest degree
    istr1 = find(abs(str1x-str1) >= ATRSH);
    istr2 = find(abs(str2x-str2) >= ATRSH);
    idip1 = find(abs(dip1x-dip1) >= ATRSH);
    idip2 = find(abs(dip2x-dip2) >= ATRSH);
    irk1 = find(abs(rk1x-rk1) >= ATRSH);
    irk2 = find(abs(rk2x-rk2) >= ATRSH);
    
    figure; nr=3; nc=2;
    subplot(nr,nc,1); plot(isub,str1x-str1,'.',isub(istr1),str1x(istr1)-str1(istr1),'ro');
    title(sprintf('strike 1 diff: %i / %i',length(istr1),n)); xlabel('GCMT index');
    subplot(nr,nc,2); plot(isub,str2x-str2,'.',isub(istr2),str2x(istr2)-str2(istr2),'ro');
    title(sprintf('strike 2 diff: %i / %i',length(istr2),n)); xlabel('GCMT index');
    subplot(nr,nc,3); plot(isub,dip1x-dip1,'.',isub(idip1),dip1x(idip1)-dip1(idip1),'ro');
    title(sprintf('dip 1 diff: %i / %i',length(idip1),n)); xlabel('GCMT index');
    subplot(nr,nc,4); plot(isub,dip2x-dip2,'.',isub(idip2),dip2x(idip2)-dip2(idip2),'ro');
    title(sprintf('dip 2 diff: %i / %i',length(idip2),n)); xlabel('GCMT index');
    subplot(nr,nc,5); plot(isub,rk1x-rk1,'.',isub(irk1),rk1x(irk1)-rk1(irk1),'ro');
    title(sprintf('rake 1 diff: %i / %i',length(irk1),n)); xlabel('GCMT index');
    subplot(nr,nc,6); plot(isub,rk2x-rk2,'.',isub(irk2),rk2x(irk2)-rk2(irk2),'ro');
    title(sprintf('rake 2 diff: %i / %i',length(irk2),n)); xlabel('GCMT index');
    orient tall, wysiwyg
end

%==========================================================================
