%
% TT2013AppA.m
%
% This script reproduces the calculations in Appendix A of
% Tape and Tape (2013), "The Classical Model of Moment Tensors"
%
% It may be helpful in understanding the equations within the paper.
%
% It uses several functions stored within the github project surfacevel2strain:
% https://github.com/carltape/surfacevel2strain/
%
% Carl Tape, 2013-10-22
%

clear, clc, close all

deg = 180/pi;

% EXAMPLE DATA: full moment tensors from Minson et al. (2007)
sdir = './';
%sdir = '/home/carltape/carltape/REPOSITORIES/mtbeach/matlab/';
ifile = [sdir 'Minson2007_Table8.txt'];
if ~exist(ifile,'file'), error('set sdir correctly'); end
[eid3,Tlam,Ttr,Tpl,Nlam,Ntr,Npl,Plam,Ptr,Ppl] = ...
textread(ifile,'%s%f%f%f%f%f%f%f%f%f','headerlines',5);
n3 = length(eid3);

% calculate the basis
Uin = [Tpl Ttr Npl Ntr Ppl Ptr]
U = U2pa(Uin,0);    % south-east-up basis
detUall = detU(U);
for ii=1:n3, disp(sprintf('%6s det(U) = %20.16f',eid3{ii},detUall(ii))); end

% check
Xcheck = U2pa(U,1)
for ii=1:n3
   disp(sprintf('%s',eid3{ii}));
   disp(sprintf('  Tpl %8.4f %8.4f',Uin(ii,1),Xcheck(ii,1)));
   disp(sprintf('  Taz %8.4f %8.4f',wrap360(Uin(ii,2)),wrap360(Xcheck(ii,2))));
   disp(sprintf('  Npl %8.4f %8.4f',Uin(ii,3),Xcheck(ii,3)));
   disp(sprintf('  Naz %8.4f %8.4f',wrap360(Uin(ii,4)),wrap360(Xcheck(ii,4))));
   disp(sprintf('  Ppl %8.4f %8.4f',Uin(ii,5),Xcheck(ii,5)));
   disp(sprintf('  Paz %8.4f %8.4f',wrap360(Uin(ii,6)),wrap360(Xcheck(ii,6))));
end

% calculate quantities from the eigenvalues
lambdaall = [Tlam Nlam Plam]' * 1e15;
MM = CMTrecom(lambdaall,U);    % south-east-up (same as input U)
[gamma3,delta3,M0,mu,lamdev,lamiso] = lam2lune(lambdaall);
nuall = lambdaall(2,:) ./ (lambdaall(1,:) + lambdaall(3,:));

% target event
ipick = 8;
fac = 1e17;
lams = lambdaall(:,ipick) / fac;
Mvec = MM(:,ipick) / fac;
Mmat = Mvec2Mmat(Mvec,1);
Umat = U(:,:,ipick);                % det(Umat) = 1
% for every U, there are three equivalent versions.
% flip sign of two columns to match Eq A3
Umat(:,[2 3]) = -Umat(:,[2 3]);     % det(Umat) = 1

if 0==1
   % shortest set of equations
   Mmat = [ 3.108 -4.855 -1.949 ;
           -4.855  3.044  1.110 ; 
           -1.949  1.110  3.382 ];
   [Umat,D] = eig(Mmat);
   lams = diag(D);
   lams = lams([3 2 1]);    % swap 1 and 3
   Umat = Umat(:,[3 2 1]);  % swap 1 and 3
   det(Umat)
   Umat(:,1) = -Umat(:,1);
   det(Umat)
   
   % check with CMTdecom.m
   Mx = Mvec2Mmat(Mmat,0);
   [lams,U] = CMTdecom(Mx,1)
   det(U)
  
   %Mvec = [3.108 3.044 3.382 -4.855 -1.949 1.110]';
   %[lams,Umat] = CMTdecom(Mvec)
   %Mmat = Mvec2Mmat(Mvec,1);
end

lam1 = lams(1);
lam2 = lams(2);
lam3 = lams(3);

disp(sprintf('======= Appendix A of TapeTape2013 (checked %s) ===========',datestr(now,29)));

% Eq A1
Mmat

% Eq A2
lams

% Eq A3
Umat

Mcheck = Umat*diag(lams)*Umat'

detU = det(Umat)

%------------
% Table A1

disp('======= Table A1 of TapeTape2013 ===========');

rho = sqrt(sum(lams.^2))

% normalized eigenvalues
lamh = lams / rho
lam1h = lamh(1);
lam2h = lamh(2);
lam3h = lamh(3);

u = (lam1-lam3) / sqrt(2)
v = (-lam1 + 2*lam2 - lam3)/sqrt(6)
w = (lam1 + lam2 + lam3)/sqrt(3)

delta = asin(w/rho)*deg         % delta3(ipick)
beta  = acos(w/rho)*deg         % 90 - delta3(ipick)
gamma = atan(v/u)*deg           % gamma3(ipick)

alpha = acos( (lam1-2*lam2+lam3) / (lam1-lam3) )*deg
nu = lam2 / (lam1 + lam3)       % nu3(ipick)

zeta = acos( sqrt(2*(lam1-lam2)*(lam2-lam3)) / rho )*deg
cosz = cos(zeta/deg);
sinz = sin(zeta/deg);

% Eqs 24
%phi = atan( (lam1-2*lam2+lam3)/(sqrt(2)*(lam1+lam2+lam3)) )*deg
phi = atan2( lam1-2*lam2+lam3, sqrt(2)*(lam1+lam2+lam3) )*deg;

% Eq S1 (not listed in table)
thetadc = acos( (lam1-lam3)/(sqrt(2)*rho) )*deg

% unit double couple tensors
D = 1/sqrt(2) * [0 0 1 ; 0 0 0 ; 1 0 0]

% unit crack tensor
K0 = [lam2 0 0 ; 0 lam2 0 ; 0 0 lam1-lam2+lam3 ];
Knorm = sqrt(sum(diag(K0).^2));     % norm(K0(:))
K = K0 / Knorm
% check with Eq 47, which uses normalized eigenvalues
% Kcheck = 1/sinz * [lam2h 0 0 ; 0 lam2h 0 ; 0 0 lam1h-lam2h+lam3h ];

% Eq 48, which uses normalized eigenvalues
Mh = [lam2h 0 sqrt((lam1h-lam2h)*(lam2h-lam3h)) ; ...
    0 lam2h 0 ; ...
    sqrt((lam1h-lam2h)*(lam2h-lam3h)) 0 lam1h-lam2h+lam3h]

t1 = sqrt((lam2-lam3)/(lam1-lam3));
t2 = sqrt((lam1-lam2)/(lam1-lam3));
Uh = [t1 0 -t2 ; 0 1 0 ; t2 0 t1]

% various checks
% norm(D(:)), norm(K(:)), norm(Mh(:)), Uh'*Uh

%------------
% Table A2

disp('======= Table A2 of TapeTape2013 ===========');

% note: Eqs 36 and 56 are the same
N1 = Umat*rotmat(-alpha/2,2)*[1 0 0]'
N2 = Umat*rotmat(alpha/2,2)*[1 0 0]'

% checks using other functions
if 0
    % check using a different function
    [nucheck,alphacheck,N1check,N2check,lamcheck] = CMT2faultpar(Mvec,true)
    % N1check = -N1 and N2check = -N2, but we can flip the sign of BOTH N1 and N2 without consequence
    % get fault angles from N and S
    [kappaX,thetaX,sigmaX,KX,NX,SX] = NS2sdr(N1,N2,true)
    
    % convert south-east-up to up-south-east
    [Muse,T] = convert_MT(5,1,Mvec);
    % calculate TT2012 parameters
    [gamma,delta,M0,kappa,theta,sigma,K,N,S,thetadc,lam,Useu] = CMT2TT(Muse);
    beta = 90 - delta
    
    % build U from the (non-orthogonal) fault vectors
    Useu
    T = (N1 + N2)/norm(N1+N2);
    P = (N1 - N2)/norm(N1-N2);
    B = cross(P,T);
    Useu_check = [T B P]
    
    % get eigenvalues, then create M
    lam_check = nualpha2lam(nu,alpha) * norm(lam);
    Mseu_check = CMTrecom(lam_check,Useu_check);
    [Mvec Mseu_check]
end

% unit double couple tensors
chi = (180 - alpha)/2;       % Eq 51
U1 = Umat*rotmat(chi,2)
U2 = Umat*rotmat(180,1)*rotmat(chi,2)
DU1 = U1*D*U1'
DU2 = U2*D*U2'

% unit crack tensors
KU1 = U1*K*U1'
KU2 = U2*K*U2'

% various checks
% norm(DU1(:)), norm(DU2(:)), norm(KU1(:)), norm(KU2(:)), U1'*U1, U2'*U2

% ultimate check
% Mmat, rho*(cosz*DU1 + sinz*KU1), rho*(cosz*DU2 + sinz*KU2)

if 1==1
    % checking with a separate function
    [K1,D1,K2,D2,KU1x,DU1x,KU2x,DU2x,rho,zeta,Ux] = CMT2CDC(Mmat);

    % Due to the different U
    Umat,Ux
    % the order of the output from CMT2CDC.m is swapped with what is shown above.
    DU1,DU2x,DU2,DU1x,KU1,KU2x,KU2,KU1x
end

% ======= Appendix A of TapeTape2013 (checked 2022-08-23) ===========
% Mmat =
%     3.1083   -4.8550   -1.9493
%    -4.8550    3.0444    1.1105
%    -1.9493    1.1105    3.3823
% lams =
%     8.8020
%     2.5840
%    -1.8510
% Umat =
%     0.6727   -0.1772    0.7184
%    -0.6391    0.3502    0.6848
%    -0.3729   -0.9198    0.1223
% Mcheck =
%     3.1083   -4.8550   -1.9493
%    -4.8550    3.0444    1.1105
%    -1.9493    1.1105    3.3823
% detU =
%     1.0000
% ======= Table A1 of TapeTape2013 ===========
% rho =
%     9.3583
% lamh =
%     0.9406
%     0.2761
%    -0.1978
% u =
%     7.5328
% v =
%    -0.7279
% w =
%     5.5050
% delta =
%    36.0329
% beta =
%    53.9671
% gamma =
%    -5.5194
% alpha =
%    80.3650
% nu =
%     0.3717
% zeta =
%    37.4790
% thetadc =
%    36.3965
% D =
%          0         0    0.7071
%          0         0         0
%     0.7071         0         0
% K =
%     0.4538         0         0
%          0    0.4538         0
%          0         0    0.7669
% Mh =
%     0.2761         0    0.5611
%          0    0.2761         0
%     0.5611         0    0.4666
% Uh =
%     0.6452         0   -0.7640
%          0    1.0000         0
%     0.7640         0    0.6452
% ======= Table A2 of TapeTape2013 ===========
% N1 =
%     0.9775
%    -0.0465
%    -0.2060
% N2 =
%     0.0504
%    -0.9301
%    -0.3638
% U1 =
%    -0.1149   -0.1772    0.9775
%    -0.9355    0.3502   -0.0465
%    -0.3340   -0.9198   -0.2060
% U2 =
%     0.9829    0.1772    0.0504
%     0.1108   -0.3502   -0.9301
%    -0.1472    0.9198   -0.3638
% DU1 =
%    -0.1588   -0.6428   -0.2141
%    -0.6428    0.0615    0.1472
%    -0.2141    0.1472    0.0973
% DU2 =
%     0.0700   -0.6425   -0.2581
%    -0.6425   -0.1457    0.0683
%    -0.2581    0.0683    0.0757
% KU1 =
%     0.7529   -0.0142   -0.0630
%    -0.0142    0.4545    0.0030
%    -0.0630    0.0030    0.4671
% KU2 =
%     0.4546   -0.0147   -0.0057
%    -0.0147    0.7247    0.1060
%    -0.0057    0.1060    0.4952

%==========================================================================
