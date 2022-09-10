function zeta = uniform_zeta(n,brandom)
%UNIFORM_ZETA get a set of zeta angles that correspond to a uniform distribution of moment tensors
%
% INPUT 
%   n       number of zeta values
% optional:
%   brand   =true for random points; =false for regular spacing
%
% OUTPUT
%   zeta    crack fraction within CDC model, degrees [0,90]
%
% EXAMPLE: zeta = uniform_zeta(1e6); [set bfigure = true]
%
% The gist of this code is to generate samples of a given probability
% density function. Presumably there are built-in functions in Matlab to do
% this kind of thing.
%
% See TapeTape2013 "The classical model for moment tensors".
%
% Carl Tape, 2022-08-30
%

if nargin==1, brandom = false; end

bfigure = false;

% inverse transform mapping
% (probably there is a more elegant way to do this)
zetai = linspace(0,pi/2,1000);
% probability density function for random moment tensors (fixed magnitude)
p = 4*cos(zetai).^3.*sin(zetai);
% integrate to cumulative distribution
pcum = cumsum(p) / max(cumsum(p));
pcum(end) = [];
zetai(end) = [];
p(end) = [];

if brandom
    u = rand(n,1);
else
    u = linspace(0,1,n); 
end
zeta = interp1(pcum,zetai,u);
zeta = zeta(:);

if bfigure
    figure; plot(zetai,pcum);
    figure; plot(pcum,zetai);
    figure; hold on;
    plot_histo(zeta,linspace(0,pi/2,40),3,true);
    plot(zetai,p,'r','linewidth',2);
end

% convert to degrees
zeta = zeta*180/pi;
