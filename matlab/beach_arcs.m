%
% beach_arcs.m
%
% Plot points, arcs, and patches for fundamental lune plots in lune . pl .
% The text files are are used for GMT plotting and may be useful for testing other plotting codes . 
% This file is included for archival purposes and is not intended to work for others .
%
% Carl Tape, 2012-08-22
%

clear, clc, close all

deg = 180/pi;

bwrite = false;

% odir = '/home/carltape/PROJECTS/cmt/figs/bfiles/';
odir = '/home/carltape/REPOSITORIES/REPOSITORIES/mtbeach/gmt/dfiles/';

% min and max limits for each parameter
u1 = 0;         u2 = 3*pi/4;    % similar to lune latitude (ISO)
v1 = -1/3;      v2 = 1/3;       % similar to lune longitude (CLVD)
k1 = 0;         k2 = 360;       % strike
s1 = -90;       s2 = 90;        % slip (rake)
h1 = 0;         h2 = 1;         % cos(dip)

% get GCDC points
nus = [0.25 0.36];
nu1 = nus(1);
nu2 = nus(2);
[gammatop1,deltatop1] = nu2gcdc(nu1);
[gammatop2,deltatop2] = nu2gcdc(nu2);

% all possible values for longitudes and latitudes
lln = [0 30];
llt = [0
    deg*asin(1/sqrt(3))
    deg*asin(sqrt(2/3))
    deltatop1(1)
    deltatop1(2)
    deltatop1(3)
    90
    deltatop2(1)
    deltatop2(2)
    deltatop2(3)
];

% reference points
nextra = 6;
plabs1 = {'ISO','-','CLVD','LVD','ISO','-','CLVD','LVD',...
    'DC',sprintf('C(nu=%.2f)',nu1),sprintf('C(nu=%.2f)',nu1),'-','-','-','-',...
    sprintf('C(nu=%.2f)',nu2),sprintf('C(nu=%.2f)',nu2),'-','-','-','-'};
% locations
lon0 = [lln(1) -lln(2) -lln(2) -lln(2) lln(1) lln(2) lln(2) lln(2) ...
    lln(1) -lln(2) lln(2) lln(2) lln(2) -lln(2) -lln(2) ...
    -lln(2) lln(2) lln(2) lln(2) -lln(2) -lln(2)]';
lat0 = [-llt(7) -llt(3) llt(1) llt(2) llt(7) llt(3) llt(1) -llt(2) ...
    llt(1) llt(4) -llt(4) llt(5) llt(6) -llt(5) -llt(6) ...
    llt(8) -llt(8) llt(9) llt(10) -llt(9) -llt(10)]';
np = length(lon0);
% pstext: L, C, R (for left, center, or right) and T, M, B for top, middle, or bottom
textalign = repmat(cellstr(''),np,1);
textalign(lat0==90) = {'CB'};
textalign(lat0==-90) = {'CT'};
textalign(lon0==-30) = {'RM'};
textalign(lon0==30) = {'LM'};
textalign(and(lon0==0,lat0==0)) = {'LB'};  % DC
% shift 
dinc = 10;  % pts 
dx = dinc*ones(1,np);
dx(lat0==90) = 0;
dx(lat0==-90) = 0;
dx(lon0==-30) = -dinc;
dx(lon0==30) = dinc;
dx(and(lon0==0,lat0==0)) = dinc;  % DC
dy = dinc*ones(1,np);
dy(lat0==90) = dinc;
dy(lat0==-90) = -dinc;
dy(lon0==-30) = 0;
dy(lon0==30) = 0;
dy(and(lon0==0,lat0==0)) = dinc;  % DC

% compute eigenvalues
lam = lune2lam(lon0,lat0,1*ones(np,1));
lnorm = NaN(1,np);
for ii=1:np
    ltemp = lam(:,ii);
    ltemp(abs(ltemp) < 1e-4) = [];
    lnorm(ii) = min(abs(ltemp));
end

% normalize to nice integer values
lamint = lam ./ repmat(lnorm,3,1);
%lamint(:,np-1:np) = lamint(:,np-1:np) * 9;  % (16,9,9) for nu=0.36
[lon0 lat0 lamint']

% eigenvalue labels
plabs2 = repmat(cellstr('-'),1,np);
for ii=1:np-nextra
    % note: rounding is somewhat dangerous, since we are assuming that the
    % eigenvalues have been normalized all to integer values
    plabs2(ii) = cellstr(sprintf('(%i,%i,%i)',round(lamint(:,ii))));
end

disp('REFERENCE POINTS FOR SOURCE TYPE:');
for ii=1:np
   disp(sprintf('%3i %12s %12s (%6.1f,%6.1f)',...
       ii,plabs1{ii},plabs2{ii},lon0(ii),lat0(ii)));
end

if bwrite
    filename = sprintf('%ssourcetype_points_lune.dat',odir);
    fid = fopen(filename,'w');
    for kk = 1:np
        fprintf(fid,'%12.4f%12.4f%12s%12s%6s%12.2f%12.2f\n',...
            lon0(kk),lat0(kk),plabs2{kk},plabs1{kk},textalign{kk},dx(kk),dy(kk));   
    end
    fclose(fid);
    
    [v0,w0] = lune2rect(lon0,lat0)
    filename = sprintf('%ssourcetype_points_rect.dat',odir);
    fid = fopen(filename,'w');
    for kk = 1:np
        fprintf(fid,'%12.4f%12.4f%12s%12s%6s%12.2f%12.2f\n',...
            v0(kk),w0(kk),plabs2{kk},plabs1{kk},textalign{kk},dx(kk),dy(kk));   
    end
    fclose(fid);
    
%     filename = sprintf('%ssourcetype_points_supp.dat',odir);
%     fid = fopen(filename,'w');
%     for kk = 12:np
%         fprintf(fid,'%12.4f%12.4f%12s%12s%12.2f%12.2f\n',...
%             lon0(kk),lat0(kk),plabs2{kk},plabs1{kk},dx(kk),dy(kk));   
%     end
%     fclose(fid);
end

% endpoints for arcs
% dev, DC+iso, top-iso, bot-iso, CDC(nu=1/4), lam2=0 (nu=0), CDC(nu=0.36)
ip = [3 7 ; 1 5 ; 4 6 ; 2 8 ; 10 11 ; 4 8 ; 16 17 ];
[na,~] = size(ip);

figure; hold on;
plot(lon0,lat0,'ko','markersize',6,'markerfacecolor','r');
text(lon0,lat0,plabs2,'fontsize',10);

figure; hold on;
plot(lon0,lat0,'ko','markersize',6,'markerfacecolor','r');
text(lon0,lat0,num2str([1:np]'),'fontsize',10);

n = 50;    % number of points to use to discretize each arc
for ii=1:na
    i1 = ip(ii,1);
    i2 = ip(ii,2);
    [lat,lon] = track2(lat0(i1),lon0(i1),lat0(i2),lon0(i2),[],'degree',n);
    
    plot(lon,lat);
    [v,w] = lune2rect(lon,lat);
    if bwrite
        filename = sprintf('%ssourcetype_arc_%2.2i.dat',odir,ii);
        write_sourcetype(filename,lon,lat,v,w);
    end
end

% patches (solid top, white bottom)
if bwrite
    % indices of the patch -- THREE vertices
    % these zones are separated by lam_i = 0 arcs
    ips = [4 6 5 4 ; 8 2 1 8 ];
    [np,nk] = size(ips);
    for ii=1:np
        xgamma = []; xdelta = []; xv = []; xw = [];
        for kk=1:nk-1
            i1 = ips(ii,kk);
            i2 = ips(ii,kk+1);
            [lat,lon] = track2(lat0(i1),lon0(i1),lat0(i2),lon0(i2),[],'degree',n);
            [v,w] = lune2rect(lon,lat);
            xgamma = [xgamma ; lon];
            xdelta = [xdelta ; lat];
            xv = [xv ; v];
            xw = [xw ; w];
        end
        filename = sprintf('%ssourcetype_patch_%2.2i.dat',odir,ii);
        write_sourcetype(filename,xgamma,xdelta,xv,xw);
    end
    
    % accessible region for nu1 and nu2
    % indices of the patch FOUR vertices
    ips = [9 13 12 10 9 ; 9 15 14 11 9 ; 
        9 19 18 16 9 ; 9 21 20 17 9];
    [np,nk] = size(ips);
    for ii=1:np     % loop over different nu values
        inu = round(ii/2);
        iarc = mod(ii-1,2)+1;
        xgamma = []; xdelta = []; xv = []; xw = [];
        for kk=1:nk-1
            i1 = ips(ii,kk);
            i2 = ips(ii,kk+1);
            [lat,lon] = track2(lat0(i1),lon0(i1),lat0(i2),lon0(i2),[],'degree',n);
            [v,w] = lune2rect(lon,lat);
            xgamma = [xgamma ; lon];
            xdelta = [xdelta ; lat];
            xv = [xv ; v];
            xw = [xw ; w];
        end
        filename = sprintf('%ssourcetype_patch_nu0p%2.2i_%2.2i.dat',odir,round(100*nus(inu)),iarc);
        write_sourcetype(filename,xgamma,xdelta,xv,xw);
    end

end

%--------------------------------------------------------------------------

irefpts = 1;            % choose a reference set of lune points
M0ref = mw2m0(1,5);     % specify a magnitude for all beachballs

switch irefpts
    case 1
        % reference moment tensors on the sphere
        dinc = 15;
        deps = 1e-3;
        %gvec = [(-30+deps):dinc:(30-deps)]; % avoid the exact lune boundary
        gvec = [-30:dinc:30];
        dvec = [(-90+dinc):15:(90-dinc)];   
        [X,Y] = meshgrid(gvec,dvec);
        gam = X(:);
        del = Y(:);
        % add the poles, but avoid the exact poles (by using deps)
        gam = [gam ; 0 ; 0];
        del = [del ; -90+deps ; 90-deps];
        n = length(gam);
        M0 = M0ref * ones(n,1);
        lam = lune2lam(gam,del,M0);
        %Uref = zeros(3,3,n); Uref(1,1,:) = 1; Uref(2,2,:) = 1; Uref(3,3,:) = 1; 
        %M = CMTrecom(lam,Uref);
        % reference basis for center
        %Uref = eye(3);
        
    case 2
        % the arc for lam2 = 0
        dp1 = [1 0 0]';
        rotax = [0 1 0]';
        n = 11;
        angvec = linspace(0,90,n);
        lam = zeros(3,n);
        for ii=1:n
            U = rotmat_gen(rotax,angvec(ii));
            lam(:,ii) = U*dp1;
        end
        % cut the dipoles, which do not seem to plot properly in psmeca
        lam(:,[1 n]) = []; n = n-2;
        lam = lam*M0ref*sqrt(2);
        [gam,del,M0] = lam2lune(lam);
        
    case 3
        
        a = [1 1 1];
        b = [1 0 0];
        c = [2 -1 -1];
        d = [0 -1 -1];
        e = [-1 -1 -1];
        f = [0 0 -1];
        g = [1 1 -2];
        h = [1 1 0];
        i = [1 0 -1];

        % for beachballs near the isotropic regions, we push them JUST
        % outside so that we can see something besides an all-red beachball
        amod = lune2lam(0,89.999);       % kludge for GMT plotting problem
        bmod = lune2lam(-30,33.2644);    % 35.2644
        dmod = lune2lam(-30,-52.7356);   % -54.7356
        emod = lune2lam(0,-89.999);      % kludge for GMT plotting problem
        fmod = lune2lam(30,-33.2644);    % -35.2644
        hmod = lune2lam(30,52.7356);     % 54.7356
 
        Lam_boundary = [
            lamarc(amod,bmod,2) ...
            lamarc(b,c,3) ...
            lamarc(c,d,4) ...
            lamarc(dmod,emod,2) ...
            lamarc(emod,fmod,2) ...
            lamarc(f,g,3) ...
            lamarc(g,h,4) ...
            lamarc(hmod,a,2) ...
            ];
        % remove: b, c, d, e, f, g
        Lam_boundary(:,[3 6 9 12 14 16 20 22]) = [];

        % keep DC once here
        lamnu1 = nualpha2lam(0.25,0);
        lamnu2 = nualpha2lam(0.25,180);
        %Lam_nu = lamarc(lamnu1,lamnu2,7);
        %Lam_nu(:,[1 4 7]) = [];
        Lam_nu = lamarc(lamnu1,lamnu2,2);

        Lam_2z = lamarc(b,f,7);
        Lam_2z(:,[1 4 7]) = [];
        
        Lam_dev = lamarc(c,g,7);
        Lam_dev(:,[1 7]) = [];

        Lam_isop = lamarc(a,i,5); Lam_isop(:,[1 end]) = [];
        Lam_ison = lamarc(i,e,5); Lam_ison(:,[1 end]) = [];
        
        lam = [Lam_nu Lam_isop Lam_ison Lam_dev Lam_2z Lam_boundary];
        [~,n] = size(lam);

        lam = lam*M0ref*sqrt(2);
        [gam,del,M0] = lam2lune(lam);

        %figure; hold on;
        %plot(gamma,delta,'bo'); grid on; axis([-31 31 -91 91]);

        %Uref = Mvec2Mmat(TT2CMT(0,0,1,90,90,0),1)
        %Uref = [0 0 0 ; 0 0 1 ; 0 1 0];
        %M = transform_MT(Uref,[Lam_all ; zeros(3,nlam)]);
        
    case 4
        NPT = 5;

        % dev
        delta0 = 0;
        gamma = linspace(-30,0,NPT);
        gamma([end]) = [];
        lam1 = lune2lam(gamma,delta0*ones(size(gamma)));
        
        % DC + iso
        gamma0 = 0;
        delta = linspace(0,90,NPT+1);
        delta([1 end]) = [];
        lam2 = lune2lam(gamma0*ones(size(delta)),delta);
        
        % NW boundary of lune
        %zeta0 = 90;
        %phi = linspace(0,90,NPT+1); phi(1) = [];
        %lam3 = phizeta2lam(phi,zeta0*ones(size(phi)));
        gamma0 = -30;
        delta = linspace(0,90,NPT+1);
        delta([1 end]) = [];
        lam3 = lune2lam(gamma0*ones(size(delta)),delta);

        % nu=0.25 in NW quadrant of lune
        nu0 = 0.25;
        lam0 = nualpha2lam(nu0,0);
        [phi0,zeta0] = lam2phizeta(lam0);
        zeta = linspace(0,90,NPT);
        lam4 = phizeta2lam(phi0*ones(size(zeta)),zeta);

        lam = [lam1 lam2 lam3 lam4];
        [~,n] = size(lam);

        lam = lam*M0ref*sqrt(2);
        [gam,del,M0] = lam2lune(lam);
        
    case 5
        gamma0 = -30;
        %n = 10; delta = linspace(0,90,n+1); delta(end) = [];
        n = 21; delta = linspace(-90,90,n+2); delta([1 end]) = [];
        lam = lune2lam(gamma0*ones(size(delta)),delta);
        lam = lam*M0ref*sqrt(2);
        [gam,del,M0] = lam2lune(lam);
        
    case 6
        gamma0 = 30;
        %n = 10; delta = linspace(0,-90,n+1); delta(end) = [];
        n = 21; delta = linspace(-90,90,n+2); delta([1 end]) = [];
        lam = lune2lam(gamma0*ones(size(delta)),delta);
        lam = lam*M0ref*sqrt(2);
        [gam,del,M0] = lam2lune(lam);

end

[v,w] = lune2rect(gam,del);
[nu,alpha] = lam2nualpha(lam);
[phi,zeta] = lam2phizeta(lam);

% remove numerical artifacts
ulam = normc(lam);
ulam(find(abs(ulam) < 0.001)) = 0;
del(find(abs(del) < 0.001)) = 0;
gam(find(abs(gam) < 0.001)) = 0;
nu(find(abs(nu) < 0.001)) = 0;
v(find(abs(v) < 0.001)) = 0;
nu(find(and(gam==0,del==0))) = 0;   % nu is undefined at DC
phi(abs(abs(phi)-180) < 0.01) = 180;

figure; plot(gam,del,'.'); axis equal; axis tight;

% loop over reference orientations
for kk=1:5
    % CHOICE OF REFERENCE ORIENTATION 
    if kk==1, str = 0; dip = 45; rk = 90; end
    if kk==2, str = 0; dip = 45; rk = -90; end
    if kk==3, str = 45; dip = 90; rk = 0; end
    if kk==4, str = -45; dip = 90; rk = 0; end
    if kk==5, str = 40; dip = 50; rk = 60; end
    MDC = dcfaultpar2CMT(str,dip,rk,0);
    [~,Uref] = CMTdecom(MDC);
    Uref
    Mdiag = [lam ; zeros(3,n)];
    M = transform_MT(Uref,Mdiag)

    % psmeca does not handle tiny numbers well
    % -- it can even flip the sign of the mechanism!
    izero = find(abs(M) < 1e-12*M0ref);
    M(izero) = 0;
    
    if bwrite
        % FUTURE WORK: add set of files for rect (vw) coordinates
        ftag = sprintf('%sbeachballs_ipts%i_iref%i',odir,irefpts,kk);
        %write_psmeca(filename,datenum(0,0,0)*ones(n,1),del,gam,0*ones(n,1),M);
        
        for jj=1:5
            % custom text labels for beachballs
            slabel = cellstr(repmat(' ',n,1));
            switch jj
                case 1, slabeltag = 'lam';
                    for ii=1:n, slabel(ii) = cellstr(sprintf('%.2f,%.2f,%.2f',ulam(:,ii))); end
                    % remove text labels for deviatoric arc
                    %if n > 9, slabel([9 10 12 13]) = cellstr('.'); end
                    
                case 2, slabeltag = 'gammadelta';
                    for ii=1:n, slabel(ii) = cellstr(sprintf('%.0f, %.0f',gam(ii),del(ii))); end
                    
                case 3, slabeltag = 'alphanu';
                    for ii=1:n
                        if abs(nu(ii)-0.25) < 0.01
                            stfmt = '%.0f, %.2f';
                        elseif or(abs(abs(nu(ii))-1) < 0.01, abs(nu(ii)) < 0.1)
                            stfmt = '%.0f, %.0f';
                        else
                            stfmt = '%.0f, %.1f';
                        end
                        slabel(ii) = cellstr(sprintf(stfmt,alpha(ii),nu(ii)));
                    end
                    slabel(find(and(gam==0,del==0))) = {'0,-'};
                    slabel(find(abs(abs(del)-90) < 0.01)) = cellstr('.');
                    
                case 4, slabeltag = 'zetaphi';
                    for ii=1:n, slabel(ii) = cellstr(sprintf('%.0f, %.0f',zeta(ii),phi(ii))); end
                    %slabel(find(abs(abs(del)-90) < 0.01)) = cellstr('.');
                    slabel(find(and(gam==0,del==0))) = {'0,-'};
                    
                case 5, slabeltag = 'vw';
                    for ii=1:n
                        if and(abs(v(ii)) < 0.01, abs(w(ii)) < 0.01)
                            stfmt = '%.0f, %.0f';
                        elseif abs(v(ii)) < 0.01
                            stfmt = '%.0f, %.2f';
                        elseif abs(w(ii)) < 0.01
                            stfmt = '%.2f, %.0f';
                        else
                            stfmt = '%.2f, %.2f';
                        end
                        slabel(ii) = cellstr(sprintf(stfmt,v(ii),w(ii)));
                    end
                    %slabel(find(abs(abs(del)-90) < 0.01)) = cellstr('.');
            end
            
            write_psmeca([ftag '_lune'],datenum(0,0,0)*ones(n,1),del,gam,0*ones(n,1),M,[],slabel,slabeltag);
            write_psmeca([ftag '_vw'],datenum(0,0,0)*ones(n,1),w,v,0*ones(n,1),M,[],slabel,slabeltag);
            if 1==1
                % write to text file for other plotting codes besides psmeca
                ofile = [ftag '.txt'];
                fid1 = fopen(ofile,'w');
                for xx=1:length(del)
                    fprintf(fid1,'%14.6f%14.6f%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n',...
                        gam(xx),del(xx),M(1,xx),M(2,xx),M(3,xx),M(4,xx),M(5,xx),M(6,xx));
                end
                fclose(fid1);
            end
        end
        
    end
end

error

%==========================================================================

a = [1 1 1];
b = [1 0 0];
c = [2 -1 -1];
d = [0 -1 -1];
e = [-1 -1 -1];
f = [0 0 -1];
g = [1 1 -2];
h = [1 1 0];

Lam_boundary = [lamarc(a,b,2) lamarc(b,c,3) lamarc(c,d,4) lamarc(e,f,2) ...
    lamarc(f,g,3) lamarc(g,h,4) ];
Lam_boundary(:,[3 6 12 15]) = []

Lam_nu = lamarc(b,f,7);
Lam_nu(:,[1 end]) = [];

Lam_dev = lamarc(c,g,7);
Lam_dev(:,[1 end]) = [];

Lam_all = [Lam_boundary Lam_nu Lam_dev];
[~,nlam] = size(Lam_all);

[gamma,delta] = lam2lune(Lam_all);
figure; hold on;
plot(gamma,delta,'bo'); grid on; axis([-31 31 -91 91]);

Uref = Mvec2Mmat(TT2CMT(0,0,1,90,90,0),1)
Uref = [0 0 0 ; 0 0 1 ; 0 1 0];
M = transform_MT(Uref,[Lam_all ; zeros(3,nlam)]);

if bwrite
    % future: add set of files for rect (vw) coordinates
    filename = sprintf('%sbeachballs_new_ipts%i_iref%i_lune',odir,irefpts,kk);
    write_psmeca(filename,datenum(0,0,0)*ones(nlam,1),delta,gamma,0*ones(nlam,1),M);
end

%==========================================================================
