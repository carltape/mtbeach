function [nu,fnu] = pdf_nu(n)
%PDF_NU probability density function for nu (Poisson parameter) for uniform moment tensors
%
% INPUT
% optional
%   n       number of points
%
% OUTPUT
%   nu      mathematical Poisson parameter (-infinity, infinity)
%   fnu     f(nu)
%
% See TapeTape2013 "The classical model for moment tensors"
% called by plotMT_nualpha.m
%
% Carl Tape, 2022-09-09
%

if nargin==0, n = 500; end

% note: the plotting range for nu is tricky, since nu can range from -infinity to infinity
% note: for real materials, nu varies from -1 to 0.5
NUMAX = 5;
NUMIN = -NUMAX;
nu = linspace(NUMIN,NUMAX,n);

if 1
    % Untitled-3.nb -- this is for an arbitrary interval nu = [NUMIN,NUMAX]
    qplot = 1 - 2*nu + 3*nu.^2;
    K = sqrt(2) / (atan((1-3*NUMIN)/sqrt(2)) +  atan((-1+3*NUMAX)/sqrt(2)));
    fnu = K ./ qplot;
else
    % NuForCarloTOSS.pdf
    qplot = 1./ ((nu - 1/3).^2 + 2/9);
    % crude integration (should be analytical)
    dx = nu(2) - nu(1);
    C = 1 / (sum(qplot) * dx);
    fnu = C*qplot;
    % maximum of distribution
    %umax = C*9/2;
end