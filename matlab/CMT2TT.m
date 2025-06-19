function [gamma,delta,M0,kappa,theta,sigma,K,N,S,thetadc,lam,U] = CMT2TT(M,bdisplay)
%CMT2TT converts a moment tensor to six parameters of TapeTape2012beach
%
% INPUT
%   M           6 x n moment tensors in CMT convention (UP-SOUTH-EAST)
%               M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%   bdisplay    OPTIONAL (if present, will display details)
%               
% OUTPUT
%   gamma       angle from DC meridian to MT point (-30 <= gamma <= 30)
%   delta       angle from deviatoric plane to MT point (-90 <= delta <= 90)
%   M0          seismic moment, N-m
%   kappa       strike angle, degrees: [0,360]
%   theta       dip angle, degrees: [0,90]
%   sigma       rake (or slip) angle, degrees: [-90,90]
% optional:
%   K           strike vector (SOUTH-EAST-UP)
%   N           normal vector (SOUTH-EAST-UP)
%   S           slip vector (SOUTH-EAST-UP)
%   thetadc     angle to DC
%   lam         eigenvalues
%   U           basis (SOUTH-EAST-UP)
%
% Reverse program for TT2CMT.m
% See TapeTape2012beach "A geometric setting for moment tensors"
%
% Carl Tape, 2012-12-01
%

if nargin==1, bdisplay=false; end

BIGN = 1e5;

% make sure M is 6 x n
[M,n] = Mdim(M);

disp(sprintf('CMT2TT.m: %i moment tensors to convert into lune + strike/dip/rake',n));

% KEY: convert M into another basis
% YOU MUST ALSO CHANGE north AND zenith IN faultvec2ang BELOW
% --> U will be with respect to this basis (from CMTdecom.m)
%M = convert_MT(1,3,M);
M = convert_MT(1,5,M);  % moment tensor in south-east-up basis

%---------------------
% PART 1: moment tensor source type (or pattern)

% decomposition of moment tensor into eigenvalues + basis (M = U*lam*U')
% NOTE: ordering of eigenvalues is important
if n >= BIGN, disp('CMT2TT: M to lam + U...'); end
isort = 1;
[lam,U] = CMTdecom(M,isort);

% compute lune coordinates and magnitude from eigenvalues
if n >= BIGN, disp('CMT2TT: lam to lune...'); end
[gamma,delta,M0,thetadc] = lam2lune(lam);

%---------------------
% PART 2: moment tensor orientation
% TT2012beach Section 6.3

[kappa,theta,sigma,K,N,S] = U2sdr(U,bdisplay);

%==========================================================================
% EXAMPLE

if 0==1
    clear, close all, clc
    n = 10;
    M0 = 1*ones(n,1);
    % random source types
    gamma = randomvec(-30,30,n); b = randomvec(-1,1,n); delta = asin(b)*180/pi;
    % random orientations
    kappa = randomvec(0,360,n); h = randomvec(0,1,n); sigma = randomvec(-90,90,n);
    theta = acos(h)*180/pi;
    %kappa = -10; theta = 30; sigma = 45; gamma = -20; delta = 70; M0 = 1;
    %kappa = -10; theta = 30; sigma = 45; gamma = -20; delta = 70; M0 = 1;
    M = TT2CMT(gamma,delta,M0,kappa,theta,sigma);
    
    [gammac,deltac,M0c,kappac,thetac,sigmac,K,N,S,thetadc,lam,U] = CMT2TT(M,1);
    disp([gamma gammac delta deltac M0 M0c kappa kappac theta thetac sigma sigmac]);
    
    % horizontal fault
    kappa = 30; theta = 0; sigma = 310;
    M = TT2CMT(0,0,1,kappa,theta,sigma)
    M = [0 0 0 -sqrt(3)/2 1/2 0]'
    [gammac,deltac,M0c,kappac,thetac,sigmac,K,N,S,thetadc,lam,U] = CMT2TT(M,1);
    disp([kappa kappac theta thetac sigma sigmac]);
    
    % GCMT catalog
    [otime,tshift,hdur,lat,lon,dep,M,M0,Mw,eid] = readCMT;
    display_eq_summary(otime,lon,lat,dep,Mw);
    [gammac,deltac,M0c,kappac,thetac,sigmac] = CMT2TT(M);
    
    % Mij listed in a USGS xml file (Mrr, Mtt, Mpp, Mrt, Mrp, Mtp; units in N-m)
    % https://earthquake.usgs.gov/product/moment-tensor/us_7000ar1r_mwr/us/1601139517040/quakeml.xml
    % 2020-07-19 00:35:27 Alaska event
    M = [-159500000000000 828490000000000 -668990000000000 244900000000000 -314060000000000 -161080000000000];
    [gamma,delta,M0,kappa,theta,sigma] = CMT2TT(M,1)
    %Mw = m02mw(1,M0)
    %plot_beachballs(M)
end

%==========================================================================
